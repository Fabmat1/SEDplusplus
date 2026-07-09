#!/usr/bin/env python3
"""bulk_fit.py -- bulk SED fitting of the local parquet SED tables with
sedfit, trying several grid / initial-value candidates per star and keeping
the best solution.

Everything is declared in one campaign JSON (see configs/campaign_example.json):
paths (sedfit, refdata, SED table), a star selection, shared fit settings, and
a list of *candidates* -- each a named combination of model grid(s) and
initial parameter values (`par` / `par_full`, same syntax as a sedfit
config). Every selected star is fitted once per candidate and the winner is
the candidate with the smallest (norm_chi_red, chisqr_reduced): sedfit
normalises chi2_red > 1 fits to 1 by adding excess noise `norm_chi_red`, so
that pair orders solutions by actual fit quality.

Pipeline per HEALPix pixel file (fully offline):
  1. Optional per-pixel joins: the local Gaia store (`gaia_store`) attaches
     corrected parallaxes (fallback: uncorrected) + GSP-Phot priors; the
     local 3D dust map (`reddening: "3d"` + `extinction_maps`) attaches the
     Wang et al. 2025 E(44-55) at the star's distance and the line-of-sight
     total. Band columns are classified via `catalogs` (catalogs_use.txt;
     catname vs filter system) when configured.
  2. sed_local.prepare_star writes photometry.dat + astrometry header per
     star.
  3. Fits: with `ms_first`, the named (main-sequence) candidate runs first
     for every star, Teff initialised from Gaia GSP-Phot and logg from a
     ZAMS estimate; stars whose fit is good (excess noise below
     accept_norm_chi_red) AND plausible against the MS mass-radius relation
     skip the remaining candidates. Everything else (or without `ms_first`:
     everything) is fitted with every candidate. Config lists are grouped by
     candidate so sedfit --multi reloads each model grid once per process,
     split over N parallel sedfit processes.
  4. Results are collected into one parquet (or CSV) per pixel:
     identity, per-candidate chi2 quality, ms_accepted, best candidate, its
     fitted parameters (+ confidence bounds) and MC radius/mass/luminosity.

Usage:
  python scripts/bulk_fit.py configs/campaign_example.json
  python scripts/bulk_fit.py campaign.json --pixels 300 301 --max-stars 50
"""
import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sed_local  # noqa: E402

NAN = float("nan")


# ---------------------------------------------------------------------------
# campaign config
# ---------------------------------------------------------------------------

DEFAULT_FIT = {
    "conf_level": 0,
    "apply_ZPO_corr": 1,
    "remove_outliers": 5,
    # 20k MC samples keep the R/M/L mode/HDI stable to ~1% of the interval
    # width -- bulk-appropriate, unlike the publication-grade 2e6 default.
    "nMC": 20000,
    "mc_seed": 42,
    # degenerate bulk fits can "improve" by micro-amounts forever in the
    # conf-limit loop; cap the restart cycles (sedfit default 1000 = ISIS)
    "max_conf_restarts": 10,
}


def load_campaign(path):
    c = json.loads(Path(path).read_text())
    for key in ("sedfit", "refdata", "sed_dir", "work_root", "out_root",
                "candidates"):
        if key not in c:
            sys.exit(f"campaign config: missing '{key}'")
    names = [cand["name"] for cand in c["candidates"]]
    if len(set(names)) != len(names):
        sys.exit("campaign config: candidate names must be unique")
    c.setdefault("jobs", 4)
    c.setdefault("distance", "parallax")     # parallax | bj_geo | bj_pgeo
    c.setdefault("reddening", "sfd")         # sfd | none | 3d | <E(B-V)>
    c.setdefault("keep_blends", False)
    c.setdefault("selection", {})
    # optional local refdata / catalogue paths (fully offline operation):
    c.setdefault("catalogs", None)          # catalogs_use.txt
    c.setdefault("gaia_store", None)        # big_surveys/store/gaia_edr3
    c.setdefault("extinction_maps", None)   # 3D dust maps (reddening: "3d")
    if c["reddening"] == "3d" and not c["extinction_maps"]:
        sys.exit("campaign config: reddening '3d' needs 'extinction_maps'")
    if c["catalogs"]:
        sed_local.load_catalogs(c["catalogs"])
    # optional MS-first mode: fit the named candidate for every star first
    # (Teff init from Gaia GSP-Phot, ZAMS logg); stars whose MS fit is good
    # AND plausible against the MS mass-radius relation skip the remaining
    # candidates.
    ms = c.get("ms_first")
    if ms:
        if ms.get("candidate") not in {n for n in names}:
            sys.exit("campaign config: ms_first.candidate must name a "
                     "candidate")
        # thresholds tuned for faint survey stars (weak parallaxes make the
        # MC mass/radius noisy and photometric logg is barely constrained)
        ms.setdefault("accept_norm_chi_red", 0.1)  # mag of excess noise
        ms.setdefault("mr_tol_dex", 0.3)   # |log10 R - log10 R_MS(M)| limit
        ms.setdefault("logg_range", [3.4, 5.6])  # fallback when R/M are NaN
    c["fit"] = {**DEFAULT_FIT, **c.get("fit", {})}
    return c


# ---------------------------------------------------------------------------
# main-sequence helpers (ms_first mode)
# ---------------------------------------------------------------------------

# crude ZAMS logg(Teff) (piecewise-linear in log10 Teff; from ZAMS M-R-Teff
# tables) -- only an initial value for the fit, accuracy ~0.1 dex suffices
_ZAMS_LOGG = [(3000.0, 4.85), (4000.0, 4.65), (5800.0, 4.44), (7300.0, 4.30),
              (9800.0, 4.15), (20000.0, 4.05), (40000.0, 3.95)]


def zams_logg(teff):
    pts = _ZAMS_LOGG
    lt = math.log10(max(teff, 1.0))
    if lt <= math.log10(pts[0][0]):
        return pts[0][1]
    for (t0, g0), (t1, g1) in zip(pts, pts[1:]):
        if teff <= t1:
            f = (lt - math.log10(t0)) / (math.log10(t1) - math.log10(t0))
            return g0 + f * (g1 - g0)
    return pts[-1][1]


def ms_radius(mass):
    """MS mass-radius relation (Demircan & Kahraman 1991), R/Rsun."""
    return 1.06 * mass ** 0.945 if mass < 1.66 else 1.33 * mass ** 0.555


def ms_plausible(res, ms):
    """Is a fit result consistent with a main-sequence star?  Mass-radius
    check when the MC gave both; else the fitted logg must look dwarf-like."""
    m, r = res.get("p_c1_M", NAN), res.get("p_c1_R", NAN)
    if math.isfinite(m) and math.isfinite(r) and m > 0 and r > 0:
        return abs(math.log10(r) - math.log10(ms_radius(m))) \
            <= ms["mr_tol_dex"]
    lo, hi = ms["logg_range"]
    logg = res.get("p_c1_logg", NAN)
    return math.isfinite(logg) and lo <= logg <= hi


def ms_accept(res, ms):
    """Full acceptance: fit quality (excess noise below the threshold) and
    MS plausibility."""
    if res is None:
        return False
    n = res["norm_chi_red"]
    if not (math.isfinite(n) and n <= ms["accept_norm_chi_red"]):
        return False
    return ms_plausible(res, ms)


_grid_trange_cache = {}


def grid_ranges(campaign, cand):
    """(teff_min, teff_max, logg_min, logg_max) of the candidate's first
    model grid (grid.fits), for clamping per-star initial values."""
    key = cand["name"]
    if key in _grid_trange_cache:
        return _grid_trange_cache[key]
    gdir = cand["griddirectories"][0]
    for bp in cand.get("bpaths", ["./"]):
        p = Path(bp) / gdir / "grid.fits"
        if p.exists():
            break
    else:
        sys.exit(f"candidate {cand['name']}: grid.fits not found under "
                 f"bpaths")
    from astropy.io import fits as pf
    import numpy as np
    with pf.open(p) as h:
        d = h[1].data
        t = np.unique(np.concatenate([np.atleast_1d(x) for x in d["t"]]))
        g = np.unique(np.concatenate([np.atleast_1d(x) for x in d["g"]]))
    r = (float(t.min()), float(t.max()), float(g.min()), float(g.max()))
    _grid_trange_cache[key] = r
    return r


def ms_candidate_for_star(campaign, cand, meta):
    """Per-star copy of the MS candidate with Teff/logg initial values from
    Gaia GSP-Phot (clamped into the grid) and a ZAMS logg estimate."""
    tmin, tmax, gmin, gmax = grid_ranges(campaign, cand)
    teff = meta.get("teff_gspphot")
    if teff is None or not math.isfinite(teff):
        return cand  # no prior: keep the configured initial values
    teff = min(max(teff, tmin), tmax)
    logg = min(max(zams_logg(teff), gmin), gmax)
    out = dict(cand)
    par = {k: list(v) for k, v in cand["par"].items()}
    for name, val in (("c*_teff", round(teff, 1)), ("c*_logg", round(logg, 3))):
        if name in par["name"]:
            par["value"][par["name"].index(name)] = val
        else:
            par["name"].append(name)
            par["value"].append(val)
            par["freeze"].append(0)
    out["par"] = par
    return out


def select_rows(campaign, args):
    """Yield (pixel, row) for the campaign's star selection."""
    sel = dict(campaign["selection"])
    if args.pixels:
        sel["pixels"] = args.pixels
        sel.pop("ids", None)
    if args.max_stars:
        sel["max_stars"] = args.max_stars
    sed = campaign["sed_dir"]

    def passes(row):
        if sel.get("max_ruwe") is not None:
            ruwe = row.get("ruwe")
            if ruwe is not None and math.isfinite(ruwe) \
                    and ruwe > sel["max_ruwe"]:
                return False
        if sel.get("min_surveys") is not None:
            ns = row.get("n_surveys") or 0
            if ns < sel["min_surveys"]:
                return False
        return True

    n = 0
    limit = sel.get("max_stars")
    if sel.get("ids"):
        for uid in sel["ids"]:
            row = sed_local.row_by_id(sed, uid)
            if row is None:
                print(f"WARNING: u_obj_id {uid} not in {sed}", file=sys.stderr)
            elif passes(row):
                yield None, row
                n += 1
                if limit and n >= limit:
                    return
        return
    pixels = sel.get("pixels") or sed_local.list_pixels(sed)
    for pix in pixels:
        for row in sed_local.pixel_rows(sed, pix):
            if passes(row):
                yield pix, row
                n += 1
                if limit and n >= limit:
                    return


# ---------------------------------------------------------------------------
# sedfit configs + execution
# ---------------------------------------------------------------------------

def write_config(campaign, cand, stardir, meta):
    """One sedfit config for (star, candidate); returns its path."""
    fit = {**campaign["fit"], **cand.get("fit", {})}
    cfg = {
        "star": f"uid{meta['u_obj_id']}",
        "coordinates": {"ra": meta["ra"], "dec": meta["dec"]},
        "griddirectories": cand["griddirectories"],
        "bpaths": cand.get("bpaths", ["./"]),
        "refdata": campaign["refdata"],
        "workdir": str(stardir),
        "outdir": str(stardir / f"cand_{cand['name']}"),
        **fit,
    }
    if "par" in cand:
        cfg["par"] = cand["par"]
    if "par_full" in cand:
        cfg["par_full"] = cand["par_full"]
    if "fix_distance" in meta:
        cfg["fix_distance"] = meta["fix_distance"]
        cfg["fix_distance_err"] = meta["fix_distance_err"]
    path = stardir / f"cand_{cand['name']}.json"
    path.write_text(json.dumps(cfg, indent=1))
    return path


def run_multi(campaign, config_paths, log_dir):
    """Run the configs through `sedfit --multi`, `jobs` processes in
    parallel. Config order is preserved within each job list (grouped by
    candidate by the caller -> each process reloads grids rarely)."""
    jobs = max(1, int(campaign["jobs"]))
    log_dir.mkdir(parents=True, exist_ok=True)
    chunks = [config_paths[i::jobs] for i in range(jobs)]
    procs = []
    for i, chunk in enumerate(c for c in chunks if c):
        lst = log_dir / f"configs_{i}.txt"
        lst.write_text("".join(str(p) + "\n" for p in chunk))
        log = open(log_dir / f"job_{i}.log", "w")
        procs.append((subprocess.Popen(
            [campaign["sedfit"], "--multi", str(lst)],
            stdout=log, stderr=subprocess.STDOUT), log))
    fails = 0
    for p, log in procs:
        fails += p.wait() != 0
        log.close()
    if fails:
        print(f"WARNING: {fails} sedfit process(es) reported failed configs "
              f"(see {log_dir}/job_*.log)", file=sys.stderr)


# ---------------------------------------------------------------------------
# result parsing
# ---------------------------------------------------------------------------

def _header_kv(lines):
    kv = {}
    for line in lines:
        if not line.startswith("#"):
            break
        for part in line[1:].split(";"):
            if "=" in part:
                k, _, v = part.partition("=")
                kv[k.strip()] = v.strip()
    return kv


def parse_result(outdir):
    """photometry_results.txt (+ stellar MC files) of one candidate run ->
    flat dict, or None if the fit did not produce results."""
    f = Path(outdir) / "photometry_results.txt"
    if not f.exists():
        return None
    lines = f.read_text().splitlines()
    kv = _header_kv(lines)
    out = {
        "norm_chi_red": float(kv.get("norm_chi_red", "nan")),
        "chisqr_reduced": float(kv.get("chisqr_reduced", "nan")),
        "nmag_good": int(float(kv.get("nmag_good", "0"))),
        "grid": kv.get("grid", ""),
    }
    in_table = False
    for line in lines:
        if line.startswith("index"):
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        tok = line.split()
        if len(tok) < 10:
            continue
        name = tok[1]
        if name.startswith("dummy_"):
            continue
        out[f"p_{name}"] = float(tok[2])
        out[f"p_{name}_lo"] = float(tok[6])
        out[f"p_{name}_hi"] = float(tok[7])
    for st in sorted(Path(outdir).glob("photometry_results_stellar_c*.txt")):
        for line in st.read_text().splitlines()[1:]:
            tok = line.split()
            if len(tok) >= 4:
                out[f"p_{tok[0]}"] = float(tok[1])
                out[f"p_{tok[0]}_lo"] = float(tok[2])
                out[f"p_{tok[0]}_hi"] = float(tok[3])
    return out


def quality_key(res):
    """Ordering key: smaller is better. Failed/NaN fits sort last."""
    n, c = res["norm_chi_red"], res["chisqr_reduced"]
    if not (math.isfinite(n) and math.isfinite(c)):
        return (math.inf, math.inf)
    return (n, c)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def write_table(rows, out_path):
    if not rows:
        print(f"no results for {out_path}")
        return
    cols = []
    for r in rows:  # union of keys, first-seen order
        for k in r:
            if k not in cols:
                cols.append(k)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        table = pa.table({c: [r.get(c) for r in rows] for c in cols})
        pq.write_table(table, out_path.with_suffix(".parquet"))
        print(f"wrote {out_path.with_suffix('.parquet')} ({len(rows)} stars)")
    except Exception as e:  # pyarrow missing/typing trouble: CSV fallback
        import csv
        with open(out_path.with_suffix(".csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out_path.with_suffix('.csv')} ({len(rows)} stars) "
              f"[parquet failed: {e}]")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("campaign", help="campaign JSON")
    ap.add_argument("--pixels", type=int, nargs="+",
                    help="override the campaign's pixel selection")
    ap.add_argument("--max-stars", type=int,
                    help="override the campaign's max_stars")
    ap.add_argument("--clean", action="store_true",
                    help="delete per-star workdirs after collecting results")
    args = ap.parse_args()
    campaign = load_campaign(args.campaign)

    work_root = Path(campaign["work_root"])
    out_root = Path(campaign["out_root"])
    out_root.mkdir(parents=True, exist_ok=True)

    # ---- 1. prepare star workdirs (grouped by pixel for the output tables
    #         and for the batched local Gaia-store / 3D-dust-map joins)
    t0 = time.time()
    rows_by_pixel = {}
    for pix, row in select_rows(campaign, args):
        rows_by_pixel.setdefault(pix, []).append(row)
    by_pixel = {}   # pixel -> [(stardir, meta)]
    n_skip = 0
    for pix, rows in rows_by_pixel.items():
        if campaign["gaia_store"]:
            if pix is None:  # by-id selection: rows carry their own pixel
                for row in rows:
                    sed_local.attach_gaia([row], campaign["gaia_store"],
                                          row["_pixel"])
            else:
                sed_local.attach_gaia(rows, campaign["gaia_store"], pix)
        if campaign["reddening"] == "3d":
            sed_local.attach_reddening_3d(rows, campaign["extinction_maps"])
        pixdir = work_root / (f"hpx_{pix:07d}" if pix is not None else "byid")
        for row in rows:
            stardir = pixdir / str(int(row["u_obj_id"]))
            meta = sed_local.prepare_star(
                row, stardir, campaign["distance"], campaign["reddening"],
                campaign["keep_blends"])
            if meta is None:
                n_skip += 1
                continue
            by_pixel.setdefault(pix, []).append((stardir, meta))
    n_stars = sum(len(v) for v in by_pixel.values())
    print(f"prepared {n_stars} stars ({n_skip} skipped: no usable "
          f"photometry) in {time.time() - t0:.1f}s")
    if n_stars == 0:
        return

    # ---- 2. fits: either MS-first with early exit, or all candidates
    ms = campaign.get("ms_first")
    accepted = {}   # stardir -> True if the MS fit was accepted
    nfits = 0
    t0 = time.time()
    if ms:
        mscand = next(c for c in campaign["candidates"]
                      if c["name"] == ms["candidate"])
        cfgs = []
        for stars in by_pixel.values():
            for stardir, meta in stars:
                cfgs.append(write_config(
                    campaign, ms_candidate_for_star(campaign, mscand, meta),
                    stardir, meta))
        run_multi(campaign, cfgs, work_root / "logs_ms")
        nfits += len(cfgs)
        for stars in by_pixel.values():
            for stardir, meta in stars:
                r = parse_result(stardir / f"cand_{ms['candidate']}")
                accepted[stardir] = ms_accept(r, ms)
        n_acc = sum(accepted.values())
        print(f"MS-first: {n_acc}/{n_stars} stars accepted by the "
              f"'{ms['candidate']}' fit; trying the remaining candidates "
              f"for the other {n_stars - n_acc}")
        rest = [c for c in campaign["candidates"]
                if c["name"] != ms["candidate"]]
        cfgs = []
        for cand in rest:
            for stars in by_pixel.values():
                for stardir, meta in stars:
                    if not accepted[stardir]:
                        cfgs.append(write_config(campaign, cand, stardir,
                                                 meta))
        if cfgs:
            run_multi(campaign, cfgs, work_root / "logs")
            nfits += len(cfgs)
    else:
        cfgs = []
        for cand in campaign["candidates"]:
            for stars in by_pixel.values():
                for stardir, meta in stars:
                    cfgs.append(write_config(campaign, cand, stardir, meta))
        run_multi(campaign, cfgs, work_root / "logs")
        nfits += len(cfgs)
    dt = time.time() - t0
    print(f"{nfits} fits for {n_stars} stars in {dt:.1f}s  "
          f"[{dt / max(nfits, 1):.2f}s/fit]")

    # ---- 3. collect, pick best per star, write per-pixel tables
    for pix, stars in by_pixel.items():
        rows = []
        for stardir, meta in stars:
            results = {}
            for cand in campaign["candidates"]:
                r = parse_result(stardir / f"cand_{cand['name']}")
                if r is not None:
                    results[cand["name"]] = r
            row = {"u_obj_id": meta["u_obj_id"],
                   "source_id": meta["source_id"],
                   "ra": meta["ra"], "dec": meta["dec"],
                   "ruwe": meta["ruwe"],
                   "n_good_bands": meta["n_good_bands"]}
            if ms:
                row["ms_accepted"] = bool(accepted.get(stardir, False))
            for name in (c["name"] for c in campaign["candidates"]):
                r = results.get(name)
                row[f"cand_{name}_norm_chi_red"] = \
                    r["norm_chi_red"] if r else NAN
                row[f"cand_{name}_chisqr_reduced"] = \
                    r["chisqr_reduced"] if r else NAN
            if ms and accepted.get(stardir):
                row["best_candidate"] = ms["candidate"]
                row.update(results[ms["candidate"]])
            elif results:
                best = min(results, key=lambda k: quality_key(results[k]))
                row["best_candidate"] = best
                row.update(results[best])
            else:
                row["best_candidate"] = ""
            rows.append(row)
            if args.clean:
                shutil.rmtree(stardir, ignore_errors=True)
        name = f"hpx_{pix:07d}" if pix is not None else "byid"
        write_table(rows, out_root / name)


if __name__ == "__main__":
    main()
