#!/usr/bin/env python3
"""sed_local.py -- feed sedfit from the local SED parquet tables (no remote
queries).

sedfit already runs fully offline when its workdir contains a cached
`photometry.dat` and a `photometry_results.txt` astrometry header (that is
exactly the cache the live-query path writes). This module generates those two
files per star from mdorsch's per-HEALPix SED tables
(sed_qmost/ / sed_xgaia/, see their README.md), so the existing, validated
fit engine is used unchanged and the remote-query mode keeps working.

Per star it writes into a workdir:
  photometry.dat          photometric_table.sl format; one row per band with
                          finite _mag; parquet _flag kept as-is (a bitmask,
                          0 = pristine: 1 survey_bad, 2 superseded,
                          4 internal_bad, 8 sed_outlier; sedfit fits only
                          flag == 0); blends (survey nWithin > 1) get
                          bit 4 (internal_bad) unless --keep-blends.
  photometry_results.txt  "# key = value;" astrometry header (RA/DEC,
                          parallax [zero-point-corrected where available],
                          ruwe) -- omitted in the BJ-distance modes, where
                          the distance goes into the config as fix_distance.
  star_meta.json          everything the bulk driver needs (ids, coords,
                          distances, band counts).

Reddening (the E(44-55) first guess sedfit takes from the photometry.dat
header): 'table' = the photometry pipeline's precomputed 3D-dust reddening
(reddening_enrich.py shards keyed by u_obj_id, e.g. photometry/
reddening_qmost/, or the same E4455_bj/E4455_los columns merged into the SED
table itself as in sed_qmost; E(44-55) at the star's distance as meanStilism,
line-of-sight total as meanSandF -- sedfit picks the LOS value itself for
stars without significant parallax), 'sfd' = local 2D SFD map
via `dustmaps` (meanSandF = 0.86 * SFD), a fixed value, or none.

Fully offline extras (see the bulk driver / CLI options): load_catalogs()
reads the photometry pipeline's catalogs_use.txt so band columns are
classified by (catname, system, band) instead of the builtin map;
attach_gaia() joins a pixel with the local Gaia store for corrected
parallaxes (parallax_corr; callers fall back to the uncorrected value where
NaN) and the GSP-Phot priors; attach_reddening_table() joins the precomputed
reddening shards.

CLI (one star / one pixel; the bulk driver imports the functions instead):
  python scripts/sed_local.py --sed <sed_dir> --id <u_obj_id> --workdir W
  python scripts/sed_local.py --sed <sed_dir> --pixel 300 --workroot ROOT
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# parquet band column -> sedfit (system, passband, VizieR_catalog)
# ---------------------------------------------------------------------------

# Preferred source: the photometry pipeline's config/catalogs_use.txt
# (one row per (catalogue, catname, system, band); catname can differ from
# system, e.g. APASS9 measures the SDSS and Johnson systems). Loaded with
# load_catalogs(); the hardcoded map below is the fallback when no
# catalogs_use.txt is configured.

_catalogs = None      # prefix -> (catname, system, passband, catalogue)


def _sanitize(s):
    """build_sed.py column-name rule: non-alphanumerics -> '_'."""
    return re.sub(r"[^0-9A-Za-z]+", "_", s)


def load_catalogs(path):
    """Parse catalogs_use.txt into the band-column lookup used by
    classify_band. SED-table column prefixes follow build_sed.py:
    '<system>_<band>' when catname == system, else
    '<catname>_<system>_<band>'; a legacy '<catname>_<band>' spelling is
    accepted too (tables built before a filter-system rename, e.g. DES_g
    from when the DES system was not yet called DECam)."""
    global _catalogs
    _catalogs = {}
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            tok = line.split()
            if header is None:
                header = tok  # "catalogue catname ... system ... filter ..."
                icat = header.index("catalogue")
                icatname = header.index("catname")
                isystem = header.index("system")
                iband = header.index("filter")
                continue
            catalogue, catname = tok[icat], tok[icatname]
            system, band = tok[isystem], tok[iband]
            entry = (catname, system, band, catalogue)
            if catname == system:
                keys = [f"{system}_{band}"]
            else:
                keys = [f"{catname}_{system}_{band}", f"{catname}_{band}"]
            for key in keys:
                _catalogs.setdefault(_sanitize(key), entry)


def _classify_from_catalogs(prefix):
    hit = _catalogs.get(prefix)
    if hit is None:
        return None
    catname, system, band, catalogue = hit
    if system == "GaiaXP":  # synthetic tophat boxes -> sedfit "box" system
        system = "box"
    return (catname, system, band, catalogue)


# hardcoded fallback map (no catalogs_use.txt configured); systems must exist
# in src/zp_table.inc (or be "box"); catalogs are provenance strings in
# photometry.dat (any value is accepted by sedfit).
_SIMPLE_SURVEYS = {
    # survey prefix: (system, catalog)
    "Gaia": ("Gaia", "I/350/gaiaedr3"),
    "GALEX": ("GALEX", "II/335/galex_ais"),
    "SkyMapper": ("SkyMapper", "Skymapper_DR4"),
    "SPLUS": ("SPLUS", "II/380/splusdr4"),
    "PS1": ("PS1", "PS1_DR2"),
    "DES": ("DES_DR2", "DES_DR2"),
    "VST": ("VST", "II/385/vst_atlas4"),
    "2MASS": ("2MASS", "II/246/out"),
    "DENIS": ("DENIS", "B/denis/denis"),
    "WISE": ("WISE", "II/363/unwise"),
    "Spitzer": ("Spitzer", "II/368/sstsl2"),
    "SDSS": ("SDSS", "sdss_dr17"),
}

_warned = set()


def classify_band(prefix):
    """Map a parquet band prefix (column name minus '_mag') to
    (survey, system, passband, catalog); None if unknown. survey is the
    catname (the '<catname>_angDist'/'_nWithin' column prefix)."""
    if _catalogs is not None:
        return _classify_from_catalogs(prefix)
    if prefix.startswith("GaiaXP_"):
        return ("GaiaXP", "box", prefix[len("GaiaXP_"):], "GaiaXP")
    if prefix.startswith("APASS9_SDSS_"):
        return ("APASS9", "SDSS", prefix.rsplit("_", 1)[1], "II/336/apass9")
    if prefix.startswith("APASS9_Johnson_"):
        return ("APASS9", "Johnson", prefix.rsplit("_", 1)[1], "II/336/apass9")
    if prefix.startswith("SAGE_Spitzer_"):
        return ("SAGE", "Spitzer", prefix[len("SAGE_Spitzer_"):],
                "II/305/catalog")
    if prefix.startswith("Legacy_DECam_"):
        # DECam AB mags; the generic zero-ZP "DES" rows + DES filter curves
        return ("Legacy", "DES", prefix[len("Legacy_DECam_"):], "legacy_dr10")
    survey, _, band = prefix.partition("_")
    if survey in _SIMPLE_SURVEYS and band:
        system, catalog = _SIMPLE_SURVEYS[survey]
        return (survey, system, band, catalog)
    if prefix not in _warned:
        _warned.add(prefix)
        print(f"sed_local: unknown band column '{prefix}_mag' -- skipped",
              file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# reddening
# ---------------------------------------------------------------------------

_DUSTMAPS_DATA = Path(__file__).resolve().parent.parent / "data" / "dustmaps"
_sfd_query = None


def sfd_ebv(ra, dec):
    """E(B-V) on the SFD scale from the local dustmaps SFD map."""
    global _sfd_query
    if _sfd_query is None:
        from dustmaps.config import config as _dm_config
        _dm_config["data_dir"] = str(_DUSTMAPS_DATA)
        from dustmaps.sfd import SFDQuery
        _sfd_query = SFDQuery()
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    return float(_sfd_query(c))


def reddening_header(row, source):
    """The reddening header lines for photometry.dat, or ''.
    source: 'sfd' | 'none' | 'table' | a float (E(B-V), SFD scale).

    'table' uses the fields attached from the precomputed reddening shards
    (or the SED table's own merged columns) by attach_reddening_table.
    sedfit turns meanStilism into the E(44-55) initial guess via
    E = E3d/0.91 (and falls back to meanSandF for stars beyond parallax
    significance), so E(44-55) is written on the 0.91 scale here."""
    if source == "none":
        return ""
    if source == "table":
        e_star = row.get("E4455_3d")
        e_los = row.get("E4455_los_3d")
        lines = ""
        if e_star is not None and math.isfinite(e_star):
            lines += (f"# meanStilism = {0.91 * e_star:.4f} "
                      f"stdStilism = 0.0000\n")
        if e_los is not None and math.isfinite(e_los):
            lines += (f"# meanSandF = {0.91 * e_los:.4f} "
                      f"stdSandF = 0.0000\n")
        return lines
    ra, dec = row["ra"], row["dec"]
    ebv = sfd_ebv(ra, dec) if source == "sfd" else float(source)
    if not math.isfinite(ebv):
        return ""
    return (f"# meanSFD = {ebv:.4f} stdSFD = 0.0000\n"
            f"# meanSandF = {0.86 * ebv:.4f} stdSandF = 0.0000\n")


def attach_reddening_table(rows, tables_dir, pixel=None):
    """Attach precomputed 3D-dust reddening (reddening_enrich.py shards, e.g.
    photometry/reddening_qmost/: one hpx_NNNNNNN.parquet per Norder=4 pixel,
    keyed by u_obj_id, E(44-55) values -- no scale conversion needed) to each
    row: E4455_3d = E4455_bj (at the star's distance) and E4455_los_3d =
    E4455_los (total sightline). Rows that already carry E4455_bj/E4455_los
    (SED tables enriched via reddening_enrich.py --to-sed, e.g. sed_qmost)
    are used directly without touching `tables_dir`.

    `pixel` is the shared HEALPix pixel of all rows; None means each row
    carries its own in row['_pixel'] (by-id selection). Stars absent from
    the shard (or in-row columns) keep NaN -> no reddening header lines
    (sedfit then treats the reddening as unconstrained)."""
    def _num(v):
        return float(v) if v is not None else float("nan")

    need = []                                # (row, pixel) still to join
    for r in rows:
        if "E4455_bj" in r or "E4455_los" in r:   # merged sed table
            r["E4455_3d"] = _num(r.get("E4455_bj"))
            r["E4455_los_3d"] = _num(r.get("E4455_los"))
        else:
            r["E4455_3d"] = r["E4455_los_3d"] = float("nan")
            need.append((r, pixel if pixel is not None else r.get("_pixel")))
    if not need:
        return
    if not tables_dir:
        sys.exit("sed_local: reddening 'table' needs the reddening-tables "
                 "dir (the SED rows carry no E4455_bj/E4455_los columns)")
    by_pix = {}
    for r, pix in need:
        by_pix.setdefault(pix, []).append(r)
    n_miss = 0
    for pix, prows in by_pix.items():
        shard = Path(tables_dir) / f"hpx_{pix:07d}.parquet"
        lut = {}
        if shard.exists():
            t = pq.read_table(shard,
                              columns=["u_obj_id", "E4455_bj", "E4455_los"])
            lut = {u: (bj, los) for u, bj, los in
                   zip(t.column("u_obj_id").to_pylist(),
                       t.column("E4455_bj").to_pylist(),
                       t.column("E4455_los").to_pylist())}
        else:
            print(f"sed_local: reddening shard missing: {shard}",
                  file=sys.stderr)
        for r in prows:
            hit = lut.get(int(r["u_obj_id"]))
            if hit is None:
                n_miss += 1
                continue
            r["E4455_3d"] = _num(hit[0])
            r["E4455_los_3d"] = _num(hit[1])
    if n_miss:
        print(f"sed_local: {n_miss}/{len(rows)} stars not in the reddening "
              f"tables ({tables_dir}) -- no reddening prior for those",
              file=sys.stderr)


# ---------------------------------------------------------------------------
# local Gaia store (astrometry + GSP-Phot priors; no remote queries)
# ---------------------------------------------------------------------------

def attach_gaia(rows, store_dir, pixel):
    """Join the rows of one SED-table HEALPix pixel with the local Gaia store
    (big_surveys/store/gaia_edr3, hats layout, same Norder=4 nested pixels):
    attaches parallax_corr (zero-point-corrected parallax; NaN where Gaia
    could not solve the correction -> callers fall back to `parallax`) and
    the GSP-Phot priors teff_gspphot / logg_gspphot / mh_gspphot."""
    if not rows:
        return
    f = (Path(store_dir) / "dataset" / "Norder=4"
         / f"Dir={(pixel // 10000) * 10000}" / f"Npix={pixel}.parquet")
    cols = ["Source", "PlxCorr", "teff_gspphot", "logg_gspphot", "mh_gspphot"]
    if not f.exists():
        print(f"sed_local: no Gaia store pixel {f} -- astrometry priors "
              f"limited to the SED table", file=sys.stderr)
        return
    t = pq.read_table(f, columns=cols)
    src = t.column("Source").to_numpy(zero_copy_only=False)
    order = np.argsort(src)
    data = {c: t.column(c).to_numpy(zero_copy_only=False)[order]
            for c in cols[1:]}
    src = src[order]
    for r in rows:
        sid = r.get("source_id")
        if sid is None:
            continue
        i = np.searchsorted(src, int(sid))
        if i >= len(src) or src[i] != int(sid):
            continue
        r["parallax_corr"] = float(data["PlxCorr"][i])
        r["teff_gspphot"] = float(data["teff_gspphot"][i])
        r["logg_gspphot"] = float(data["logg_gspphot"][i])
        r["mh_gspphot"] = float(data["mh_gspphot"][i])


# ---------------------------------------------------------------------------
# per-star files
# ---------------------------------------------------------------------------

def _f(row, key, default=float("nan")):
    v = row.get(key)
    if v is None:
        return default
    v = float(v)
    return v if math.isfinite(v) else default


def photometry_dat_text(row, redd_lines, keep_blends=False):
    """photometry.dat content for one SED-table row (dict).
    Returns (text, n_good): n_good = rows sedfit will actually fit."""
    lines = [f"# RA = {row['ra']:.10f} DEC = {row['dec']:.10f}\n", redd_lines,
             "flag system passband magnitude uncertainty type"
             " angu_dist_arcsec VizieR_catalog\n"]
    n_good = 0
    for col, val in row.items():
        if not col.endswith("_mag") or val is None:
            continue
        mag = float(val)
        if not math.isfinite(mag):
            continue
        info = classify_band(col[:-4])
        if info is None:
            continue
        survey, system, passband, catalog = info
        err = _f(row, col[:-4] + "_err")
        # flag is a bitmask (0 = pristine): 1 survey_bad, 2 superseded,
        # 4 internal_bad, 8 sed_outlier (reserved for the post-SED pass)
        flag = int(row.get(col[:-4] + "_flag") or 0)
        nwithin = _f(row, survey + "_nWithin", 1.0)
        if not keep_blends and nwithin > 1:
            flag |= 4  # blended: identification not unique -> internal_bad
        angdist = _f(row, survey + "_angDist", 0.0)
        if not math.isfinite(angdist):
            angdist = 0.0
        if flag == 0:
            n_good += 1
        lines.append(f"{flag:4d} {system} {passband} {mag:.6f} "
                     f"{err:.6f} magnitude {angdist:.6f} {catalog}\n")
    return "".join(lines), n_good


def astrometry_header_text(row):
    """photometry_results.txt astrometry-cache header (parallax prior).
    Uses the Lindegren+2021-corrected parallax where present."""
    plx = _f(row, "parallax_corr")
    if not math.isfinite(plx):
        plx = _f(row, "parallax")
    plx_err = _f(row, "parallax_error")
    ruwe = _f(row, "ruwe")
    h = (f"# RA = {row['ra']:.8f}; DEC = {row['dec']:.8f};\n"
         f"# parallax = {plx:.6g}; parallax_error = {plx_err:.6g}; ")
    if math.isfinite(ruwe):
        h += f"ruwe = {ruwe:.4g}; "
    return h + "\n"


def bj_distance(row, kind):
    """(fix_distance, fix_distance_err) in kpc from the Bailer-Jones columns,
    or None. kind: 'bj_geo' | 'bj_pgeo'."""
    base = "BJ2021_r_geo" if kind == "bj_geo" else "BJ2021_r_pgeo"
    r = _f(row, base)
    if not math.isfinite(r) or r <= 0:
        return None
    lo, hi = _f(row, base + "_lo"), _f(row, base + "_hi")
    err = (hi - lo) / 2.0 if math.isfinite(lo) and math.isfinite(hi) else 0.1 * r
    return r / 1000.0, max(err, 1e-4 * r) / 1000.0


def prepare_star(row, workdir, distance="parallax", reddening="sfd",
                 keep_blends=False, min_points=5):
    """Write photometry.dat (+ astrometry header) + star_meta.json for one
    SED-table row into `workdir`. Returns the meta dict (None if the star has
    fewer than `min_points` usable photometry points -- too few to constrain
    the SED, and outlier removal can push the fit below the number of free
    parameters)."""
    redd = reddening_header(row, reddening)
    text, n_good = photometry_dat_text(row, redd, keep_blends)
    if n_good < min_points:
        return None
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "photometry.dat").write_text(text)

    meta = {
        "u_obj_id": int(row["u_obj_id"]),
        "source_id": (int(row["source_id"])
                      if row.get("source_id") is not None else None),
        "ra": float(row["ra"]), "dec": float(row["dec"]),
        "n_good_bands": n_good,
        "ruwe": _f(row, "ruwe"),
        "distance_mode": distance,
        # local Gaia-store extras (NaN when absent): GSP-Phot fit priors and
        # the 3D-dust reddening, for the bulk driver's initial values
        "teff_gspphot": _f(row, "teff_gspphot"),
        "logg_gspphot": _f(row, "logg_gspphot"),
        "mh_gspphot": _f(row, "mh_gspphot"),
        "E4455_3d": _f(row, "E4455_3d"),
    }
    if distance == "parallax":
        (workdir / "photometry_results.txt").write_text(
            astrometry_header_text(row))
    elif distance in ("bj_geo", "bj_pgeo"):
        d = bj_distance(row, distance)
        if d:
            meta["fix_distance"], meta["fix_distance_err"] = d
        else:  # no BJ distance: fall back to the parallax prior
            meta["distance_mode"] = "parallax"
            (workdir / "photometry_results.txt").write_text(
                astrometry_header_text(row))
    # distance == "none": no header, no fix_distance -> sedfit fits without a
    # distance prior only if coordinates come from the config; the driver
    # always sets them.
    (workdir / "star_meta.json").write_text(json.dumps(meta, indent=1))
    return meta


# ---------------------------------------------------------------------------
# SED-table access (mirrors query_sed.py, without its repo imports)
# ---------------------------------------------------------------------------

def pixel_rows(sed_dir, pixel):
    """All rows of one HEALPix pixel file as dicts."""
    f = Path(sed_dir) / f"hpx_{pixel:07d}.parquet"
    return pq.read_table(f).to_pylist()


def pixel_nrows(sed_dir, pixel):
    """Row count of one pixel file (parquet metadata only; cheap)."""
    f = Path(sed_dir) / f"hpx_{pixel:07d}.parquet"
    return pq.ParquetFile(f).metadata.num_rows


def pixel_of_id(sed_dir, uid):
    """HEALPix pixel of one u_obj_id via the _index.parquet sidecar."""
    hit = pq.read_table(Path(sed_dir) / "_index.parquet",
                        filters=[("u_obj_id", "==", int(uid))])
    return int(hit.column("pixel")[0].as_py()) if hit.num_rows else None


def row_by_id(sed_dir, uid):
    """One row by u_obj_id via the _index.parquet sidecar; also returns the
    pixel as row['_pixel'] (for the local Gaia-store join)."""
    pix = pixel_of_id(sed_dir, uid)
    if pix is None:
        return None
    t = pq.read_table(Path(sed_dir) / f"hpx_{pix:07d}.parquet",
                      filters=[("u_obj_id", "==", int(uid))])
    if not t.num_rows:
        return None
    row = t.to_pylist()[0]
    row["_pixel"] = pix
    return row


def list_pixels(sed_dir):
    return sorted(int(p.stem[4:])
                  for p in Path(sed_dir).glob("hpx_*.parquet"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sed", required=True, help="sed_qmost/ or sed_xgaia/")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--id", type=int, help="one star by u_obj_id")
    g.add_argument("--pixel", type=int, help="all stars of one HEALPix pixel")
    ap.add_argument("--workdir", help="output dir (--id mode)")
    ap.add_argument("--workroot", help="root dir; one subdir per star "
                                       "(--pixel mode)")
    ap.add_argument("--distance", default="parallax",
                    choices=["parallax", "bj_geo", "bj_pgeo", "none"])
    ap.add_argument("--reddening", default="sfd",
                    help="'sfd' (local dustmaps SFD), 'table' "
                         "(precomputed reddening shards / merged E4455_bj "
                         "columns; see --reddening-tables), 'none', or a "
                         "fixed E(B-V) value on the SFD scale")
    ap.add_argument("--keep-blends", action="store_true",
                    help="keep flag 0 on bands whose survey has nWithin > 1")
    ap.add_argument("--catalogs", help="catalogs_use.txt (band-column "
                    "classification; default: builtin map)")
    ap.add_argument("--gaia-store", help="local Gaia store "
                    "(big_surveys/store/gaia_edr3): corrected parallaxes + "
                    "GSP-Phot priors")
    ap.add_argument("--reddening-tables", help="precomputed reddening shard "
                    "dir (reddening_enrich.py output, e.g. photometry/"
                    "reddening_qmost/, for --reddening table; not needed "
                    "when the SED table already carries E4455_bj/E4455_los)")
    a = ap.parse_args()

    if a.catalogs:
        load_catalogs(a.catalogs)

    if a.id is not None:
        if not a.workdir:
            ap.error("--id needs --workdir")
        row = row_by_id(a.sed, a.id)
        if row is None:
            sys.exit(f"u_obj_id {a.id} not found in {a.sed}")
        if a.gaia_store:
            attach_gaia([row], a.gaia_store, row["_pixel"])
        if a.reddening == "table":
            attach_reddening_table([row], a.reddening_tables)
        meta = prepare_star(row, a.workdir, a.distance, a.reddening,
                            a.keep_blends)
        if meta is None:
            sys.exit("star has no usable (finite, flag-0) band")
        print(json.dumps(meta, indent=1))
    else:
        if not a.workroot:
            ap.error("--pixel needs --workroot")
        rows = pixel_rows(a.sed, a.pixel)
        if a.gaia_store:
            attach_gaia(rows, a.gaia_store, a.pixel)
        if a.reddening == "table":
            attach_reddening_table(rows, a.reddening_tables, a.pixel)
        n = n_skip = 0
        for row in rows:
            uid = int(row["u_obj_id"])
            meta = prepare_star(row, Path(a.workroot) / str(uid), a.distance,
                                a.reddening, a.keep_blends)
            if meta is None:
                n_skip += 1
            else:
                n += 1
        print(f"pixel {a.pixel}: {n} stars prepared, {n_skip} without "
              f"usable photometry")


if __name__ == "__main__":
    main()
