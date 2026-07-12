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

The campaign is processed ONE HEALPix pixel at a time (bounded memory; a
100M-star table never lives in RAM), and every pixel's result table is
written before the next pixel is fitted, which makes runs resumable: a
pixel whose out_root table already exists is skipped, so interrupting and
rerunning the same command continues where the run stopped (--overwrite
refits). While one pixel is being fitted by the sedfit processes, the NEXT
pixel is loaded, joined and prepared in a background thread, so the
per-star preparation stays off the critical path.

The work tree is scratch space: every file in it (photometry.dat,
star_meta.json, the per-candidate configs and result text files) is either
regenerated from the SED table + campaign or fully parsed into the output
parquet. The script therefore cleans up after itself BY DEFAULT: once a
pixel's results are collected, its per-star workdirs, the pixel dir and its
sedfit logs are deleted, bounding the work tree to ~one pixel in flight
(essential at 100M-star scale, where keeping the workdirs costs terabytes
and hundreds of millions of inodes). Logs of a unit whose sedfit runs
failed or produced missing results are kept for debugging. Use --keep-work
(or campaign "clean": false) to keep everything instead.

Pipeline per HEALPix pixel file (fully offline):
  0. Mag-grid preprocessing (once per campaign): every candidate grid gets
     a precomputed per-band magnitude grid (<griddir>/mag/,
     `sedfit --premag`) that replaces the per-fit passband integration;
     rebuilt only when the passbands (filter archive / GaiaXP boxes)
     change. Fits default to the bulk fast paths (error_mode "covar",
     nMC 0, use_mag_grid, fast_ext) -- deliberately NOT ISIS-parity; see
     DEFAULT_FIT.
  1. Optional per-pixel joins: the local Gaia store (`gaia_store`) attaches
     corrected parallaxes (fallback: uncorrected) + GSP-Phot priors; the
     reddening join (`reddening: "table"` + `reddening_tables`, the
     photometry pipeline's precomputed per-pixel shards or the E4455 columns
     merged into the SED table itself as in sed_qmost) attaches E(44-55) at
     the star's distance and the line-of-sight total. Band columns are
     classified via `catalogs` (catalogs_use.txt; catname vs filter system)
     when configured.
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
import importlib.util
import json
import math
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# interpreter guard: the machine's default `python` is a foreign conda env
# without the campaign dependencies. When required modules are missing,
# re-exec this script with the known-good interpreter (override:
# SEDPP_PYTHON) instead of crashing mid-pipeline with a ModuleNotFoundError.
# ---------------------------------------------------------------------------

_PREFERRED_PYTHON = os.environ.get(
    "SEDPP_PYTHON", "/home/taurus/userdata/fmattig/miniforge3/bin/python3")
_REQUIRED_MODULES = ("numpy", "pyarrow", "astropy")
_OPTIONAL_MODULES = ("dustmaps", "rich", "psutil")


def _ensure_interpreter():
    missing = [m for m in _REQUIRED_MODULES + _OPTIONAL_MODULES
               if importlib.util.find_spec(m) is None]
    if not missing:
        return
    exe = Path(_PREFERRED_PYTHON)
    if (not os.environ.get("_SEDPP_REEXEC") and exe.is_file()
            and os.access(exe, os.X_OK)
            and Path(sys.executable).resolve() != exe.resolve()):
        print(f"bulk_fit: module(s) {', '.join(missing)} missing in "
              f"{sys.executable} -- re-executing with {exe}",
              file=sys.stderr)
        os.environ["_SEDPP_REEXEC"] = "1"
        os.execv(str(exe), [str(exe), str(Path(__file__).resolve()),
                            *sys.argv[1:]])
    hard = [m for m in missing if m in _REQUIRED_MODULES]
    if hard:
        sys.exit(f"bulk_fit: missing required module(s) "
                 f"{', '.join(hard)} in {sys.executable}.\n"
                 f"Run with {_PREFERRED_PYTHON} (or set SEDPP_PYTHON to an "
                 f"interpreter that has pyarrow/astropy).")
    print(f"bulk_fit: optional module(s) {', '.join(missing)} unavailable "
          f"in {sys.executable} -- continuing with reduced functionality",
          file=sys.stderr)


_ensure_interpreter()

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sed_local  # noqa: E402
from bulk_tui import Dashboard  # noqa: E402

NAN = float("nan")


# ---------------------------------------------------------------------------
# campaign config
# ---------------------------------------------------------------------------

DEFAULT_FIT = {
    "conf_level": 0,
    "apply_ZPO_corr": 1,
    "remove_outliers": 5,
    # bulk fast paths (not ISIS-parity, validated against the parity path):
    # covariance errors instead of the conf-limit loop, analytic R/M/L error
    # propagation instead of the MC, precomputed per-band mag grids instead
    # of flux interpolation + passband integration, fast 10^x extinction.
    "error_mode": "covar",
    "nMC": 0,
    "use_mag_grid": 1,
    "fast_ext": 1,
    "mc_seed": 42,
    # only relevant when a campaign switches error_mode back to "conf":
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
    # jobs: "auto" (or 0/None) = all available cores
    c.setdefault("jobs", "auto")
    if c["jobs"] in (None, 0, "auto"):
        c["jobs"] = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") \
            else (os.cpu_count() or 4)
    c.setdefault("premag", True)  # build/refresh mag grids before fitting
    c.setdefault("clean", True)   # delete per-star workdirs after collection
    c.setdefault("distance", "parallax")     # parallax | bj_geo | bj_pgeo
    c.setdefault("reddening", "sfd")     # sfd | none | table | <E(B-V)>
    c.setdefault("keep_blends", False)
    # stars with fewer good photometry points are skipped outright: the SED
    # is unconstrained and outlier removal can leave fewer points than free
    # parameters
    c.setdefault("min_points", 5)
    c.setdefault("selection", {})
    # optional local refdata / catalogue paths (fully offline operation):
    c.setdefault("catalogs", None)          # catalogs_use.txt
    c.setdefault("gaia_store", None)        # big_surveys/store/gaia_edr3
    c.setdefault("reddening_tables", None)  # precomputed reddening shards
    # reddening "table" without 'reddening_tables' is allowed when the SED
    # table itself carries the merged E4455_bj/E4455_los columns (sed_qmost);
    # sed_local errors out per pixel if neither source is available.
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
        # use the mass-radius check only when the R/M relative 1-sigma is
        # below this; noisier estimates (weak parallax) use the logg range
        ms.setdefault("mr_max_rel_err", 1.0)
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
    check when R and M are actually constrained (relative 1-sigma below
    ms_mr_max_rel_err), with the tolerance widened by the propagated
    uncertainty of log(R/R_MS(M)); else the fitted logg must look dwarf-like.
    (The analytic nMC=0 error propagation returns R/M for weak parallaxes
    where the MC used to yield NaN -- without the gate those meaningless
    point estimates would drive the acceptance.)"""
    m, r = res.get("p_c1_M", NAN), res.get("p_c1_R", NAN)
    if math.isfinite(m) and math.isfinite(r) and m > 0 and r > 0:
        def rel_halfwidth(name, value):
            lo, hi = res.get(f"{name}_lo", NAN), res.get(f"{name}_hi", NAN)
            if not (math.isfinite(lo) and math.isfinite(hi) and hi >= lo):
                return None
            return (hi - lo) / (2 * value)
        sr = rel_halfwidth("p_c1_R", r)
        sm = rel_halfwidth("p_c1_M", m)
        if sr is not None and sm is not None \
                and max(sr, sm) <= ms["mr_max_rel_err"]:
            slope = 0.945 if m < 1.66 else 0.555  # dlogR_MS/dlogM
            sig = math.hypot(sr, slope * sm) / math.log(10)
            return abs(math.log10(r) - math.log10(ms_radius(m))) \
                <= ms["mr_tol_dex"] + sig
        # R/M unconstrained (weak parallax): fall through to the logg check
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
    Gaia GSP-Phot (clamped into the grid) and a ZAMS logg estimate.
    Parameters frozen in the campaign keep their configured value (e.g. a
    campaign that freezes c*_logg is not overridden by the ZAMS estimate)."""
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
            i = par["name"].index(name)
            if par["freeze"][i]:  # frozen in the campaign: keep its value
                continue
            par["value"][i] = val
        else:
            par["name"].append(name)
            par["value"].append(val)
            par["freeze"].append(0)
    out["par"] = par
    return out


def _passes(sel, row):
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


def _result_exists(base):
    return (base.with_suffix(".parquet").exists()
            or base.with_suffix(".csv").exists())


def plan_batches(campaign, args, out_root, dash):
    """The campaign's pending work units as (selection, [(pixel, out_base)],
    star estimate): one unit per HEALPix pixel (or a single by-id unit),
    resume-skipping every unit whose result table already exists (--overwrite
    refits). The star count is a cheap pre-selection estimate from the
    parquet metadata, for the dashboard only."""
    sel = dict(campaign["selection"])
    if args.pixels:
        sel["pixels"] = args.pixels
        sel.pop("ids", None)
    if args.max_stars:
        sel["max_stars"] = args.max_stars
    sed = campaign["sed_dir"]
    if sel.get("ids"):
        base = out_root / "byid"
        if _result_exists(base) and not args.overwrite:
            dash.log("resume: the by-id results exist -- nothing to do "
                     "(--overwrite refits)")
            return sel, [], 0
        est = len(sel["ids"])
        return sel, [(None, base)], min(est, sel.get("max_stars") or est)
    pixels = sel.get("pixels") or sed_local.list_pixels(sed)
    pending, n_done = [], 0
    for pix in pixels:
        base = out_root / f"hpx_{pix:07d}"
        if _result_exists(base) and not args.overwrite:
            n_done += 1
        else:
            pending.append((pix, base))
    if n_done:
        dash.log(f"resume: {n_done}/{len(pixels)} pixel(s) already have "
                 f"results -- skipped (--overwrite refits)")
    est = sum(sed_local.pixel_nrows(sed, pix) for pix, _ in pending)
    if sel.get("max_stars"):
        est = min(est, sel["max_stars"])
    return sel, pending, est


def load_batch(campaign, sel, work_root, pix, limit, log):
    """Load, filter, join and prepare ONE work unit (a pixel, or the by-id
    selection when pix is None). Returns ([(stardir, meta)], n_selected,
    n_skipped). Runs in the prefetch thread; only log() touches the UI."""
    sed = campaign["sed_dir"]
    t0 = time.time()
    if pix is None:
        rows = []
        for uid in sel["ids"]:
            row = sed_local.row_by_id(sed, uid)
            if row is None:
                log(f"WARNING: u_obj_id {uid} not in {sed}")
            elif _passes(sel, row):
                rows.append(row)
    else:
        rows = [r for r in sed_local.pixel_rows(sed, pix) if _passes(sel, r)]
    if limit is not None:
        rows = rows[:max(limit, 0)]
    label = f"pixel {pix}" if pix is not None else "by-id selection"
    if campaign["gaia_store"]:
        if pix is None:  # by-id selection: rows carry their own pixel
            for row in rows:
                sed_local.attach_gaia([row], campaign["gaia_store"],
                                      row["_pixel"])
        else:
            sed_local.attach_gaia(rows, campaign["gaia_store"], pix)
    if campaign["reddening"] == "table":
        sed_local.attach_reddening_table(
            rows, campaign["reddening_tables"], pix)
    pixdir = work_root / (f"hpx_{pix:07d}" if pix is not None else "byid")
    stars, n_skip = [], 0
    for row in rows:
        stardir = pixdir / str(int(row["u_obj_id"]))
        meta = sed_local.prepare_star(
            row, stardir, campaign["distance"], campaign["reddening"],
            campaign["keep_blends"], campaign["min_points"])
        if meta is None:
            n_skip += 1
        else:
            stars.append((stardir, meta))
    log(f"{label}: prepared {len(stars)} stars ({n_skip} skipped: fewer "
        f"than {campaign['min_points']} usable points) in "
        f"{time.time() - t0:.1f}s")
    return stars, len(rows), n_skip


def prepared_batches(campaign, sel, pending, work_root, log):
    """Generator over the prepared work units: (pix, out_base, stars,
    n_skip). Enforces the global max_stars budget across units."""
    left = sel.get("max_stars")
    for pix, base in pending:
        stars, n_sel, n_skip = load_batch(campaign, sel, work_root, pix,
                                          left, log)
        if left is not None:
            left -= n_sel
        yield pix, base, stars, n_skip
        if left is not None and left <= 0:
            return


_SENTINEL = object()


def prefetched(gen, depth=1):
    """Iterate `gen` through a background thread that keeps up to `depth`
    items ready: the next pixel is loaded and prepared WHILE the current one
    is being fitted, so the (single-threaded) prep never blocks the fit
    processes between pixels."""
    q = queue.Queue(maxsize=depth)
    err = []

    def worker():
        try:
            for item in gen:
                q.put(item)
        except BaseException as e:   # re-raised on the consumer side
            err.append(e)
        finally:
            q.put(_SENTINEL)

    threading.Thread(target=worker, daemon=True).start()
    while True:
        item = q.get()
        if item is _SENTINEL:
            if err:
                raise err[0]
            return
        yield item


# ---------------------------------------------------------------------------
# mag-grid preprocessing (sedfit --premag)
# ---------------------------------------------------------------------------

def resolve_griddir(cand, gdir):
    """Resolve a candidate grid directory against its bpaths."""
    for bp in cand.get("bpaths", ["./"]):
        p = Path(bp) / gdir
        if (p / "grid.fits").exists():
            return p
    sys.exit(f"candidate {cand['name']}: grid {gdir} not found under bpaths")


def gaiaxp_boxes():
    """Box passband names ('l0_l1') of all GaiaXP bands known to the loaded
    catalogs table -- the fixed universe of synthetic tophat boxes that the
    mag grids must cover."""
    cats = getattr(sed_local, "_catalogs", None) or {}
    return sorted({band for (_, system, band, _) in cats.values()
                   if system == "GaiaXP"})


def ensure_mag_grids(campaign, dash):
    """Build/refresh the precomputed mag grid (<griddir>/mag/) of every
    candidate grid via `sedfit --premag`. Cheap when up to date: sedfit
    compares the manifest (filter-archive hash + box coverage) and exits
    immediately unless the passbands changed."""
    if not campaign.get("premag", True):
        return
    boxes = gaiaxp_boxes()
    dirs = sorted({str(resolve_griddir(cand, g))
                   for cand in campaign["candidates"]
                   for g in cand["griddirectories"]})
    for gdir in dirs:
        cmd = [campaign["sedfit"], "--premag", gdir, campaign["refdata"],
               *boxes]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        phase = False
        done = 0
        for line in proc.stdout:
            tok = line.split()
            if not tok or tok[0] != "PREMAG":
                continue
            if tok[1] == "UPTODATE":
                dash.log(f"mag grid up-to-date: {gdir}")
            elif tok[1] == "DONE":
                dash.log(f"mag grid ready: {gdir} ({tok[2]} nodes)")
            elif len(tok) == 3:
                if not phase:
                    dash.set_phase(f"premag {Path(gdir).parent.name}",
                                   int(tok[2]), count_stars=False, fits=False)
                    phase = True
                dash.step(int(tok[1]) - done)
                done = int(tok[1])
                dash.tick()
        if proc.wait() != 0:
            sys.exit(f"sedfit --premag failed for {gdir} "
                     f"(rerun: {' '.join(cmd)})")


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
    # a result file left over from a previous run would be mistaken for this
    # run's output (by the progress polling AND by parse_result)
    stale = stardir / f"cand_{cand['name']}" / "photometry_results.txt"
    stale.unlink(missing_ok=True)
    return path


def run_multi(campaign, config_paths, log_dir, dash=None, clean=False):
    """Run the configs through `sedfit --multi`, `jobs` processes in
    parallel. Config order is preserved within each job list (grouped by
    candidate by the caller -> each process reloads grids rarely).

    While the processes run, completion is tracked by polling each config's
    result file (<stardir>/cand_X/photometry_results.txt) and reported to
    the dashboard: fits done, plus stars done once ALL of a star's configs
    in this batch have results.

    With `clean`, log_dir is deleted when every fit succeeded; the logs of
    units with failures/missing results are kept for debugging."""
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

    # cand_X.json -> cand_X/photometry_results.txt written on fit completion
    pending = {p: p.parent / p.stem / "photometry_results.txt"
               for p in config_paths}
    cfgs_left = {}   # stardir -> configs of this batch still running
    for p in config_paths:
        cfgs_left[p.parent] = cfgs_left.get(p.parent, 0) + 1
    poll = min(5.0, max(1.0, len(pending) / 2e5))
    while True:
        alive = any(p.poll() is None for p, _ in procs)
        done = [c for c, m in pending.items() if m.exists()]
        n_stars = 0
        for c in done:
            del pending[c]
            cfgs_left[c.parent] -= 1
            if cfgs_left[c.parent] == 0:
                n_stars += 1
        if dash:
            dash.batch_advance(len(done), n_stars if dash.count_stars else 0)
            dash.tick(force=not alive)
        if not alive:
            break
        time.sleep(poll)

    fails = 0
    for p, log in procs:
        fails += p.wait() != 0
        log.close()
    warn = print if dash is None else dash.log
    if pending:
        warn(f"WARNING: {len(pending)} fit(s) produced no result file "
             f"(see {log_dir}/job_*.log)")
    if fails:
        warn(f"WARNING: {fails} sedfit process(es) reported failed configs "
             f"(see {log_dir}/job_*.log)")
    if clean and not pending and not fails:
        shutil.rmtree(log_dir, ignore_errors=True)


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

def write_table(rows, out_path, log=print):
    if not rows:
        # a 0-row marker table: the unit is DONE (e.g. every star skipped),
        # so the resume logic must not refit it forever
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
            pq.write_table(pa.table({"u_obj_id": pa.array([], pa.int64())}),
                           out_path.with_suffix(".parquet"))
        except Exception:
            out_path.with_suffix(".csv").write_text("u_obj_id\n")
        log(f"no results for {out_path} -- wrote an empty marker table")
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
        log(f"wrote {out_path.with_suffix('.parquet')} ({len(rows)} stars)")
    except Exception as e:  # pyarrow missing/typing trouble: CSV fallback
        import csv
        with open(out_path.with_suffix(".csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        log(f"wrote {out_path.with_suffix('.csv')} ({len(rows)} stars) "
            f"[parquet failed: {e}]")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("campaign", help="campaign JSON")
    ap.add_argument("--pixels", type=int, nargs="+",
                    help="override the campaign's pixel selection")
    ap.add_argument("--max-stars", type=int,
                    help="override the campaign's max_stars (a global budget "
                         "across pixels; NB a pixel fitted under a budget "
                         "writes a partial table and counts as done for "
                         "the resume logic)")
    ap.add_argument("--jobs", type=int,
                    help="override the campaign's parallel sedfit processes "
                         "(campaign default: all available cores)")
    ap.add_argument("--no-premag", action="store_true",
                    help="skip the mag-grid preprocessing check")
    ap.add_argument("--clean", action="store_true",
                    help="(default; kept for compatibility) delete per-star "
                         "workdirs + logs after collecting a unit's results")
    ap.add_argument("--keep-work", action="store_true",
                    help="keep the per-star workdirs and sedfit logs instead "
                         "of cleaning them up (debugging; also a campaign "
                         "key 'clean': false). NB a full-catalogue run "
                         "leaves terabytes of small files this way")
    ap.add_argument("--overwrite", action="store_true",
                    help="refit pixels whose result table already exists "
                         "(default: skip them = resume an interrupted run)")
    ap.add_argument("--no-tui", action="store_true",
                    help="plain status lines instead of the live dashboard")
    args = ap.parse_args()
    campaign = load_campaign(args.campaign)
    if args.jobs:
        campaign["jobs"] = args.jobs
    if args.no_premag:
        campaign["premag"] = False

    work_root = Path(campaign["work_root"])
    out_root = Path(campaign["out_root"])
    out_root.mkdir(parents=True, exist_ok=True)

    with Dashboard(enabled=not args.no_tui) as dash:
        _run_campaign(campaign, args, work_root, out_root, dash)


def _fit_batch(campaign, pix, stars, work_root, dash, clean):
    """Fit the prepared stars of one work unit; returns ({stardir:
    ms_accepted}, n_fits). MS-first semantics as before, scoped to the
    unit."""
    label = f"pixel {pix}" if pix is not None else "by-id selection"
    name = f"hpx_{pix:07d}" if pix is not None else "byid"
    ms = campaign.get("ms_first")
    accepted = {}
    nfits = 0
    if ms:
        mscand = next(c for c in campaign["candidates"]
                      if c["name"] == ms["candidate"])
        cfgs = [write_config(campaign,
                             ms_candidate_for_star(campaign, mscand, meta),
                             stardir, meta)
                for stardir, meta in stars]
        dash.set_phase(f"{label}: MS-first fits ('{ms['candidate']}')",
                       len(cfgs), count_stars=False)
        run_multi(campaign, cfgs, work_root / "logs" / f"{name}_ms", dash,
                  clean)
        nfits += len(cfgs)
        for stardir, meta in stars:
            r = parse_result(stardir / f"cand_{ms['candidate']}")
            accepted[stardir] = ms_accept(r, ms)
        n_acc = sum(accepted.values())
        dash.add_stars(n_acc)
        dash.log(f"{label}: MS-first accepted {n_acc}/{len(stars)} stars")
        rest = [c for c in campaign["candidates"]
                if c["name"] != ms["candidate"]]
        cfgs = [write_config(campaign, cand, stardir, meta)
                for cand in rest
                for stardir, meta in stars if not accepted[stardir]]
        if cfgs:
            dash.set_phase(f"{label}: remaining candidates", len(cfgs))
            run_multi(campaign, cfgs, work_root / "logs" / f"{name}_rest",
                      dash, clean)
            nfits += len(cfgs)
    else:
        cfgs = [write_config(campaign, cand, stardir, meta)
                for cand in campaign["candidates"]
                for stardir, meta in stars]
        dash.set_phase(f"{label}: fitting all candidates", len(cfgs))
        run_multi(campaign, cfgs, work_root / "logs" / name, dash, clean)
        nfits += len(cfgs)
    return accepted, nfits


def _collect_batch(campaign, pix, stars, accepted, out_base, clean, dash):
    """Collect one work unit's results, pick the best candidate per star and
    write the unit's output table (this is what marks the unit done for the
    resume logic)."""
    label = f"pixel {pix}" if pix is not None else "by-id selection"
    ms = campaign.get("ms_first")
    dash.set_phase(f"{label}: collecting results", len(stars),
                   count_stars=False, fits=False)
    rows = []
    for stardir, meta in stars:
        dash.step()
        dash.tick()
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
        if clean:
            shutil.rmtree(stardir, ignore_errors=True)
    if clean and stars:
        # every stardir of this unit is gone -- drop the (empty) pixel dir
        # too so a full-catalogue run doesn't leave one dir per pixel behind
        shutil.rmtree(stars[0][0].parent, ignore_errors=True)
    write_table(rows, out_base, dash.log)


def _run_campaign(campaign, args, work_root, out_root, dash):
    # ---- 0. mag-grid preprocessing (no-op when the passbands are unchanged)
    ensure_mag_grids(campaign, dash)

    # ---- 1. plan: one work unit per pixel; units with an existing result
    #         table are skipped (resume)
    sel, pending, est = plan_batches(campaign, args, out_root, dash)
    if not pending:
        dash.log("nothing to do")
        return
    dash.set_total_stars(est)
    dash.log(f"{len(pending)} work unit(s) pending (~{est} stars, "
             f"pre-selection estimate)")
    clean = (args.clean or campaign["clean"]) and not args.keep_work
    if not clean:
        dash.log("cleanup disabled -- the work tree will keep every "
                 "per-star dir (fine for debugging, not for large runs)")

    # ---- 2. per-unit pipeline: load+join+prepare in a background thread
    #         (one unit ahead), fit + collect + write in the main thread.
    #         Each unit's table is written before the next unit's fits, so
    #         an interrupted campaign resumes at the unit that was running.
    t0 = time.time()
    tot_stars = tot_fits = tot_skip = n_units = 0
    dash.set_phase("preparing the first work unit", None, count_stars=False,
                   fits=False)
    for pix, base, stars, n_skip in prefetched(
            prepared_batches(campaign, sel, pending, work_root, dash.log)):
        n_units += 1
        tot_skip += n_skip
        accepted = {}
        if stars:
            accepted, nfits = _fit_batch(campaign, pix, stars, work_root,
                                         dash, clean)
            tot_fits += nfits
        _collect_batch(campaign, pix, stars, accepted, base, clean, dash)
        tot_stars += len(stars)
        dash.tick(force=True)
    dt = time.time() - t0
    dash.log(f"campaign: {tot_stars} stars fitted ({tot_fits} fits, "
             f"{tot_skip} stars skipped) in {n_units} work unit(s), "
             f"{_fmt_dur(dt)}  [{dt / max(tot_fits, 1):.2f}s/fit]")
    if n_units < len(pending):
        dash.log("hit the max_stars budget -- rerun the same command to "
                 "continue with the remaining pixels")


def _fmt_dur(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


if __name__ == "__main__":
    main()
