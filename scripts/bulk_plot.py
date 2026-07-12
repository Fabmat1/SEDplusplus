#!/usr/bin/env python3
"""bulk_plot.py -- plot any bulk_fit.py output parameter against any other.

Reads the per-pixel result parquets a campaign wrote into its out_root
(out/hpx_*.parquet) and draws an x-vs-y plot of any two columns, with
optional confidence-interval error bars, log axes, a third parameter as
colour (sequential colourbar for numbers, fixed hues for categories like
best_candidate) and pandas-query row filters.

Small samples (<= 5000 stars) are drawn as per-star markers; larger ones
switch to an mpl-scatter-density raster (density colourbar, or per-pixel
mean colour / per-category density overlays with --color), which stays
fast and readable up to hundreds of millions of stars. --scatter /
--density force either mode. Note: in density mode --log-x/--log-y plot
log10 of the data on linear axes (the density artist cannot bin on log
axes), so the tick values are exponents.

Column names may be given in shorthand: 'teff' resolves to p_c1_teff, 'R' to
p_c1_R, 'E_44m55' to p_E_44m55 (tried in the order: exact, p_<name>,
p_c1_<name>). `--list` shows every available column.

Examples:
  python scripts/bulk_plot.py out/ --hrd --color best_candidate
  python scripts/bulk_plot.py out/ --hrd --debug configs/bulk_campaign.json
  python scripts/bulk_plot.py out/ -x R -y M --log-x --log-y --errors
  python scripts/bulk_plot.py out/ -x teff -y HE --color best_candidate
  python scripts/bulk_plot.py out/ -x teff -y logg --color norm_chi_red \\
      --where "best_candidate == 'sdB'" --flip-y --save teff_logg.png
  python scripts/bulk_plot.py out/ --list
"""
import argparse
import glob
import math
import os
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# result loading
# ---------------------------------------------------------------------------

def load_results(paths):
    """All result rows from parquet/CSV files, dirs of hpx_*.parquet, or a
    mix, as one pandas DataFrame with a `_file` provenance column."""
    import pandas as pd
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            hits = sorted(p.glob("hpx_*.parquet")) or sorted(
                p.glob("hpx_*.csv"))
            if not hits:
                sys.exit(f"bulk_plot: no hpx_*.parquet/csv in {p}")
            files += hits
        elif p.exists():
            files.append(p)
        else:
            hits = sorted(glob.glob(str(p)))
            if not hits:
                sys.exit(f"bulk_plot: no such file or dir: {p}")
            files += [Path(h) for h in hits]
    frames = []
    for f in files:
        df = (pd.read_csv(f) if f.suffix == ".csv" else pd.read_parquet(f))
        df["_file"] = f.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def resolve_column(df, name):
    """Map a (shorthand) name to a DataFrame column: exact, p_<name>,
    p_c1_<name>."""
    for cand in (name, f"p_{name}", f"p_c1_{name}"):
        if cand in df.columns:
            return cand
    sys.exit(f"bulk_plot: no column for '{name}' (tried {name}, p_{name}, "
             f"p_c1_{name}); see --list")


# axis labels for the common fit parameters (fallback: the column name)
_LABELS = {
    "p_c1_teff": r"$T_\mathrm{eff}$ [K]",
    "p_bb_teff": r"$T_\mathrm{eff,bb}$ [K]",
    "p_bb2_teff": r"$T_\mathrm{eff,bb2}$ [K]",
    "p_c1_logg": r"$\log g$ [cgs]",
    "p_c1_HE": r"$\log n(\mathrm{He})/n(\mathrm{H})$",
    "p_c1_z": r"$z$ [dex]",
    "p_c1_xi": r"$\xi$ [km/s]",
    "p_c1_vsini": r"$v \sin i$ [km/s]",
    "p_c1_R": r"$R$ [$R_\odot$]",
    "p_c1_M": r"$M$ [$M_\odot$]",
    "p_c1_L": r"$L$ [$L_\odot$]",
    "p_logtheta": r"$\log \Theta$",
    "p_E_44m55": r"$E(44-55)$ [mag]",
    "p_R_55": r"$R(55)$",
    "norm_chi_red": "excess noise (norm_chi_red) [mag]",
    "chisqr_reduced": r"$\chi^2_\mathrm{red}$",
    "n_good_bands": "good bands",
    "nmag_good": "fitted magnitudes",
    "ruwe": "RUWE",
}


def axis_label(col):
    return _LABELS.get(col, col)


def error_arrays(df, col, values):
    """Asymmetric (lower, upper) error-bar arrays from the <col>_lo/_hi
    confidence bounds, or None if the bounds are missing."""
    lo, hi = f"{col}_lo", f"{col}_hi"
    if lo not in df.columns or hi not in df.columns:
        return None
    el = np.clip(values - df[lo].to_numpy(float), 0, None)
    eu = np.clip(df[hi].to_numpy(float) - values, 0, None)
    el = np.where(np.isfinite(el), el, 0.0)
    eu = np.where(np.isfinite(eu), eu, 0.0)
    return el, eu


# fixed categorical hue order (validated palette; never cycled)
_CAT_COLORS = ["#2a78d6", "#1baf7a", "#eda100", "#008300",
               "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"]

# sequential blue ramp (same palette) for the plain density raster
_SEQ_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
             "#256abf", "#184f95", "#0d366b"]


def density_cmap(colors, base_alpha=1.0):
    """Colormap for a scatter_density raster: `colors` light->dark (or a
    single hue faded in from base_alpha), empty pixels transparent."""
    from matplotlib.colors import LinearSegmentedColormap, to_rgba
    if isinstance(colors, str):
        colors = [to_rgba(colors, base_alpha), to_rgba(colors, 1.0)]
    cm = LinearSegmentedColormap.from_list("bulk_density", colors)
    cm.set_bad((0, 0, 0, 0))   # log(0) pixels: show the surface
    return cm


# ---------------------------------------------------------------------------
# --debug: click a star -> single-star sedfit re-fit with full outputs
# ---------------------------------------------------------------------------

def debug_refit(campaign, row):
    """Re-fit one plotted star with the single-star sedfit and full output
    products (photometry_results.tex/.pdf, SED_results.fits, SED plot).
    The star is re-prepared from the SED table (the bulk workdir may be
    cleaned or stale), fitted with its winning candidate and the campaign's
    fit settings, into <work_root>/debug/<u_obj_id>/. Returns the path of
    the rendered SED plot (png), or None."""
    import subprocess
    os.environ.setdefault("_SEDPP_REEXEC", "1")  # no bulk_fit self re-exec
    import bulk_fit
    import sed_local

    uid = int(row["u_obj_id"])
    cand_name = row.get("best_candidate") or ""
    cands = {c["name"]: c for c in campaign["candidates"]}
    cand = cands.get(cand_name) or campaign["candidates"][0]
    print(f"\n=== re-fitting u_obj_id {uid} "
          f"(candidate '{cand['name']}') ===")

    srow = sed_local.row_by_id(campaign["sed_dir"], uid)
    if srow is None:
        print(f"debug: u_obj_id {uid} not found in {campaign['sed_dir']}")
        return None
    if campaign["gaia_store"]:
        sed_local.attach_gaia([srow], campaign["gaia_store"],
                              srow["_pixel"])
    if campaign["reddening"] == "table":
        sed_local.attach_reddening_table([srow],
                                         campaign["reddening_tables"])
    stardir = Path(campaign["work_root"]) / "debug" / str(uid)
    meta = sed_local.prepare_star(
        srow, stardir, campaign["distance"], campaign["reddening"],
        campaign["keep_blends"], campaign["min_points"])
    if meta is None:
        print(f"debug: star {uid} has too few usable photometry points")
        return None
    ms = campaign.get("ms_first")
    if ms and cand["name"] == ms["candidate"]:
        cand = bulk_fit.ms_candidate_for_star(campaign, cand, meta)
    cand = dict(cand)
    # full single-star outputs; the SED plot is rendered below with THIS
    # interpreter (sedfit's own `plot` toggle shells out to a bare python3)
    cand["fit"] = {**cand.get("fit", {}),
                   "write_tex": 1, "write_fits": 1, "write_model": 1}
    cfg = bulk_fit.write_config(campaign, cand, stardir, meta)

    r = subprocess.run([campaign["sedfit"], str(cfg)])
    outdir = stardir / f"cand_{cand['name']}"
    results = outdir / "photometry_results.txt"
    if r.returncode != 0 or not results.exists():
        print(f"debug: sedfit failed (exit {r.returncode}); config: {cfg}")
        return None
    print(results.read_text())

    png = outdir / "photometry_SED.png"
    plot_sed = Path(__file__).resolve().parent / "plot_sed.py"
    p = subprocess.run([sys.executable, str(plot_sed),
                        str(outdir / "SED_results.fits"), str(png)])
    if p.returncode != 0:
        print("debug: SED plot rendering failed")
        png = None
    print(f"debug: outputs in {outdir}/ (photometry_results.txt/.tex, "
          f"SED_results.fits{', photometry_SED.png' if png else ''})")
    return png


def install_debug(fig, ax, df, x, y, campaign_path):
    """Bind left-clicks on the scatter to debug_refit; the star's SED plot
    opens in a new figure window."""
    os.environ.setdefault("_SEDPP_REEXEC", "1")  # no bulk_fit self re-exec
    import bulk_fit
    import matplotlib.pyplot as plt
    campaign = bulk_fit.load_campaign(campaign_path)
    pick_r2 = 12.0 ** 2   # max click distance in screen pixels

    def onclick(ev):
        if ev.inaxes is not ax or ev.button != 1:
            return
        if getattr(fig.canvas.toolbar, "mode", ""):
            return   # zoom/pan drag, not a pick
        # cheap data-space prefilter around the click (transforming every
        # point to pixels is infeasible for 1e8 stars), then exact pixel
        # distance on the few candidates
        r = math.sqrt(pick_r2)
        inv = ax.transData.inverted()
        box = inv.transform([(ev.x - r, ev.y - r), (ev.x + r, ev.y + r)])
        (xlo, xhi), (ylo, yhi) = np.sort(box[:, 0]), np.sort(box[:, 1])
        idx = np.flatnonzero((x >= xlo) & (x <= xhi) &
                             (y >= ylo) & (y <= yhi))
        if idx.size == 0:
            return
        pts = ax.transData.transform(np.column_stack([x[idx], y[idx]]))
        d2 = (pts[:, 0] - ev.x) ** 2 + (pts[:, 1] - ev.y) ** 2
        i = int(np.argmin(d2))
        if d2[i] > pick_r2:
            return
        row = df.iloc[int(idx[i])]
        png = debug_refit(campaign, row)
        if png:
            import matplotlib.image as mpimg
            f2 = plt.figure(figsize=(11, 7))
            f2.canvas.manager.set_window_title(
                f"u_obj_id {int(row['u_obj_id'])}")
            a2 = f2.add_axes([0, 0, 1, 1])
            a2.imshow(mpimg.imread(png))
            a2.set_axis_off()
            f2.show()

    fig.canvas.mpl_connect("button_press_event", onclick)
    print("debug: click a star to re-fit it with full outputs "
          "(TeX + SED plot)")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples:")[1])
    ap.add_argument("results", nargs="+",
                    help="campaign out_root dir(s), hpx_*.parquet/csv "
                         "file(s), or globs")
    ap.add_argument("-x", help="x column (shorthand ok, e.g. teff, R, M)")
    ap.add_argument("-y", help="y column")
    ap.add_argument("--hrd", action="store_true",
                    help="HRD preset: L (log scale) over Teff (log scale, "
                         "decreasing rightward); -x/-y still override the "
                         "columns")
    ap.add_argument("--color", "-c",
                    help="third column as colour: numeric -> viridis "
                         "colourbar, categorical (best_candidate, grid, "
                         "ms_accepted) -> fixed hues + legend")
    ap.add_argument("--errors", action="store_true",
                    help="error bars from the <col>_lo/_hi confidence bounds")
    ap.add_argument("--where", help="pandas query filter, e.g. "
                    "\"best_candidate == 'MS' and norm_chi_red < 0.05\"")
    ap.add_argument("--log-x", action="store_true")
    ap.add_argument("--log-y", action="store_true")
    ap.add_argument("--flip-x", action="store_true",
                    help="invert the x axis (e.g. Teff in a Kiel diagram)")
    ap.add_argument("--flip-y", action="store_true",
                    help="invert the y axis (e.g. logg in a Kiel diagram)")
    ap.add_argument("--save", "-o", metavar="FILE",
                    help="write the figure (png/pdf/...) instead of showing "
                         "it; default when no display is available: "
                         "bulk_plot.png")
    ap.add_argument("--density", action="store_true",
                    help="force the scatter-density raster (default for "
                         "> 5000 stars)")
    ap.add_argument("--scatter", action="store_true",
                    help="force per-star markers even for large samples")
    ap.add_argument("--dpi", type=float, default=72.0,
                    help="density raster resolution in dots/inch "
                         "(default 72; higher = finer bins)")
    ap.add_argument("--size", type=float, default=14.0,
                    help="marker area in scatter mode (default 14)")
    ap.add_argument("--alpha", type=float,
                    help="marker opacity in scatter mode (default: auto "
                         "from star count)")
    ap.add_argument("--list", action="store_true",
                    help="list the available columns and exit")
    ap.add_argument("--debug", metavar="CAMPAIGN",
                    help="interactive outlier inspection: click a marker to "
                         "re-fit that star with the single-star sedfit "
                         "(TeX + SED plot outputs into <work_root>/debug/"
                         "<u_obj_id>/) and pop up its SED plot; needs the "
                         "campaign JSON and an interactive window (no "
                         "--save)")
    a = ap.parse_args()
    if a.hrd:
        a.x = a.x or "teff"
        a.y = a.y or "L"
        a.log_x = True
        a.log_y = True
        a.flip_x = True

    df = load_results(a.results)

    if a.list:
        print(f"{len(df)} stars from {df['_file'].nunique()} file(s); "
              f"columns (shorthand = name without p_ / p_c1_):")
        for col in df.columns:
            if col == "_file" or col.endswith(("_lo", "_hi")):
                continue
            err = " (+errors)" if f"{col}_lo" in df.columns else ""
            print(f"  {col}{err}")
        return
    if not a.x or not a.y:
        ap.error("-x and -y are required (or use --list)")

    if a.where:
        try:
            df = df.query(a.where)
        except Exception as e:
            sys.exit(f"bulk_plot: bad --where expression: {e}")
        if df.empty:
            sys.exit("bulk_plot: --where filter matches no stars")

    xcol, ycol = resolve_column(df, a.x), resolve_column(df, a.y)
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y)
    ccol = resolve_column(df, a.color) if a.color else None
    if ccol is not None and df[ccol].dtype.kind in "fc":
        good &= np.isfinite(df[ccol].to_numpy(float))
    df, x, y = df[good], x[good], y[good]
    n = len(df)
    if n == 0:
        sys.exit(f"bulk_plot: no star has finite {xcol} and {ycol}")

    if a.density and a.scatter:
        ap.error("--density and --scatter are mutually exclusive")
    use_density = a.density or (n > 5000 and not a.scatter)
    dlog_x = dlog_y = False
    if use_density and (a.log_x or a.log_y):
        # the density artist bins linearly and cannot live on log axes,
        # so plot log10 of the data instead (tick values are exponents)
        pos = np.ones(n, bool)
        if a.log_x:
            pos &= x > 0
        if a.log_y:
            pos &= y > 0
        if not pos.all():
            print(f"bulk_plot: dropping {n - int(pos.sum())} stars with "
                  f"non-positive values on a log axis", file=sys.stderr)
            df, x, y = df[pos], x[pos], y[pos]
            n = len(df)
            if n == 0:
                sys.exit("bulk_plot: no stars left after the log-axis cut")
        if a.log_x:
            x, dlog_x, a.log_x = np.log10(x), True, False
        if a.log_y:
            y, dlog_y, a.log_y = np.log10(y), True, False

    if a.debug:
        if a.save:
            sys.exit("bulk_plot: --debug is interactive; drop --save")
        if not os.environ.get("DISPLAY"):
            sys.exit("bulk_plot: --debug needs an interactive display "
                     "(X forwarding or a local session)")
        if not Path(a.debug).is_file():
            sys.exit(f"bulk_plot: no such campaign JSON: {a.debug}")
        if "u_obj_id" not in df.columns:
            sys.exit("bulk_plot: --debug needs u_obj_id in the results")
    if a.save or not os.environ.get("DISPLAY"):
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if use_density:
        import warnings
        import mpl_scatter_density  # noqa: F401 -- registers the projection
        # the artist's pre-data first draw warns about all-NaN rasters
        warnings.filterwarnings("ignore", message="All-NaN slice",
                                module="mpl_scatter_density.*")

    alpha = a.alpha if a.alpha is not None else (
        0.9 if n <= 200 else 0.5 if n <= 5000 else 0.2)
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True,
                           subplot_kw=({"projection": "scatter_density"}
                                       if use_density else {}))

    if a.errors and use_density:
        print("bulk_plot: --errors is not drawn in density mode "
              "(use --scatter with a --where filter)", file=sys.stderr)
    elif a.errors:
        for col, vals, axis in ((xcol, x, "x"), (ycol, y, "y")):
            err = error_arrays(df, col, vals)
            if err is None:
                print(f"bulk_plot: no {col}_lo/_hi bounds -- no {axis} "
                      f"error bars", file=sys.stderr)
                continue
            ax.errorbar(x, y, **{f"{axis}err": np.vstack(err)},
                        fmt="none", ecolor="0.55", elinewidth=0.8,
                        alpha=min(1.0, alpha + 0.1), zorder=1)

    mappable = None
    if use_density:
        from matplotlib.colors import LogNorm
        if ccol is None:
            d = ax.scatter_density(x, y, cmap=density_cmap(_SEQ_RAMP),
                                   norm=LogNorm(vmin=1), dpi=a.dpi,
                                   zorder=2)
            fig.canvas.draw()   # rasterize now so the norm has data
            fig.colorbar(d, ax=ax, label="stars per pixel", pad=0.01)
        elif df[ccol].dtype.kind in "fc" and df[ccol].nunique() > 8:
            cm = plt.get_cmap("viridis").copy()
            cm.set_bad((0, 0, 0, 0))    # empty pixels: show the surface
            mappable = ax.scatter_density(x, y,
                                          c=df[ccol].to_numpy(float),
                                          cmap=cm, dpi=a.dpi, zorder=2)
            fig.canvas.draw()   # rasterize now so the norm has data
            fig.colorbar(mappable, ax=ax, pad=0.01,
                         label=f"mean {axis_label(ccol)}")
        else:                       # categorical: one density per hue
            from matplotlib.lines import Line2D
            cats = sorted(df[ccol].dropna().unique(), key=str)
            if len(cats) > len(_CAT_COLORS):
                sys.exit(f"bulk_plot: {ccol} has {len(cats)} categories; "
                         f"use a numeric column or filter with --where")
            handles = []
            for cat, color in zip(cats, _CAT_COLORS):
                m = (df[ccol] == cat).to_numpy()
                if m.any():
                    ax.scatter_density(x[m], y[m],
                                       cmap=density_cmap(color, 0.6),
                                       norm=LogNorm(vmin=1), dpi=a.dpi,
                                       zorder=2)
                handles.append(Line2D([], [], marker="s", ls="none",
                                      color=color, label=str(cat)))
            # every scatter_density call reset the limits to ITS category's
            # range; restore the full data range
            ax.set_xlim(np.min(x), np.max(x))
            ax.set_ylim(np.min(y), np.max(y))
            ax.legend(handles=handles, title=axis_label(ccol),
                      frameon=False, handletextpad=0.2)
    elif ccol is None:
        ax.scatter(x, y, s=a.size, c=_CAT_COLORS[0], alpha=alpha,
                   linewidths=0, zorder=2)
    elif df[ccol].dtype.kind in "fc" and df[ccol].nunique() > 8:
        mappable = ax.scatter(x, y, s=a.size, c=df[ccol].to_numpy(float),
                              cmap="viridis", alpha=alpha, linewidths=0,
                              zorder=2)
        fig.colorbar(mappable, ax=ax, label=axis_label(ccol), pad=0.01)
    else:                           # categorical / few discrete values
        cats = sorted(df[ccol].dropna().unique(), key=str)
        if len(cats) > len(_CAT_COLORS):
            sys.exit(f"bulk_plot: {ccol} has {len(cats)} categories; use a "
                     f"numeric column or filter with --where")
        for cat, color in zip(cats, _CAT_COLORS):
            m = (df[ccol] == cat).to_numpy()
            ax.scatter(x[m], y[m], s=a.size, c=color, alpha=alpha,
                       linewidths=0, zorder=2, label=str(cat))
        ax.legend(title=axis_label(ccol), frameon=False, markerscale=1.6,
                  handletextpad=0.2)

    if a.log_x:
        ax.set_xscale("log")
    if a.log_y:
        ax.set_yscale("log")
    if a.flip_x:
        ax.invert_xaxis()
    if a.flip_y:
        ax.invert_yaxis()
    ax.set_xlabel((r"$\log_{10}$ " if dlog_x else "") + axis_label(xcol))
    ax.set_ylabel((r"$\log_{10}$ " if dlog_y else "") + axis_label(ycol))
    ax.set_title(f"{n:,} stars", fontsize=10, color="0.35", loc="right")
    ax.grid(True, color="0.9", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("0.6")
    ax.tick_params(color="0.6", labelcolor="0.25")

    out = a.save
    if not out and not os.environ.get("DISPLAY"):
        out = "bulk_plot.png"
        print("bulk_plot: no display -- writing bulk_plot.png",
              file=sys.stderr)
    if out:
        fig.savefig(out, dpi=150)
        print(f"bulk_plot: wrote {out} ({n} stars: {xcol} vs {ycol})")
    else:
        if a.debug:
            install_debug(fig, ax, df, x, y, a.debug)
        plt.show()


if __name__ == "__main__":
    main()
