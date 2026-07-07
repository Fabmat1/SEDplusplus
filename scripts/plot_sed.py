#!/usr/bin/env python3
"""SED plot from SED_results.fits — Python replacement for the xfig figure
rendered by ISIS photometry.sl:1380-1525 via xfig_photometry.sl.

Usage: plot_sed.py <SED_results.fits> <output.pdf>

Reads the extensions written by sedfit (spectrum_fit, filters, colours, IUE)
and reproduces the three-pane layout: main SED pane (model spectrum, observed
filter-averaged fluxes with error bars, filter widths and passband labels,
per-component curves), a chi_magnitude pane below, and a chi_colour side pane
when colours exist. Semantics follow xfig_photometry.sl (chi + errbar +
filter_width + colored qualifiers, exponent=3, xmin=1100, xmax=58999); exact
typography parity with xfig/LaTeX is not a goal.
"""

import sys

import numpy as np
from astropy.io import fits

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EXPONENT = 3
XMIN, XMAX = 1100.0, 58999.0  # photometry.sl:1311 (xmax-1 passed to the plot)

# xfig colour names used by filter_color_xypos -> hex (X11 values).
XFIG_COLORS = {
    "darkviolet": "#9400D3", "blue": "#0000FF", "green4": "#008B00",
    "crimson": "#DC143C", "goldenrod": "#DAA520", "royalblue": "#4169E1",
    "maroon": "#B03060", "cornflowerblue": "#6495ED", "cyan4": "#008B8B",
    "brown": "#A52A2A", "gold": "#FFD700", "steelblue": "#4682B4",
    "orange": "#FFA500", "red": "#FF0000", "pink4": "#8B636C",
    "red4": "#8B0000", "red3": "#CD0000", "blue2": "#0000EE",
    "magenta": "#FF00FF", "deeppink": "#FF1493", "purple": "#A020F0",
    "black": "#000000", "gray": "#808080", "skyblue": "#87CEEB",
    "pink3": "#CD919E", "violet": "#EE82EE",
}

# system -> (xfig colour, rel_y_pos default); rel_y <= 0 places the label
# above the error bar, > 0 below (filter_color_xypos defaults).
SYSTEM_STYLE = {
    "GALEX": ("darkviolet", -1), "TD1": ("darkviolet", 1),
    "FAUST": ("darkviolet", -1), "Johnson": ("blue", -1),
    "Stroemgren": ("green4", 1), "Geneva": ("crimson", 1),
    "SDSS": ("goldenrod", -1), "VST": ("goldenrod", -1),
    "SkyMapper": ("green4", -1.5), "PS1": ("crimson", -1),
    "ZTF": ("royalblue", -2), "INT": ("maroon", 1),
    "Gaia": ("cornflowerblue", 1), "Hipparcos": ("cyan4", 1),
    "Tycho": ("brown", -1), "BATC": ("gold", 1), "JPLUS": ("steelblue", -1),
    "SPLUS": ("steelblue", -1), "DES": ("gold", 1), "DES_DR2": ("gold", 1),
    "DELVE": ("gold", 1), "DECaPS": ("gold", 1), "SMASH": ("gold", 1),
    "DENIS": ("orange", -1), "2MASS": ("red", -1), "UKIDSS": ("pink4", -1),
    "WFCAM": ("pink4", -1), "NSFCam": ("red4", 1), "ERIS_NIX": ("red3", -1),
    "AstroSat": ("blue2", -1), "WISE": ("magenta", -1), "HST": ("black", -1),
    "SWIFT": ("deeppink", 1), "XMM": ("deeppink", 1), "Spitzer": ("purple", 1),
    "VISTA": ("red4", -1), "box": ("black", -1),
}


def pretty_passband(system, passband):
    """Label text per filter_color_xypos (subset relevant to sedfit output)."""
    p = passband
    if system == "box":
        return "$.$"
    if system in ("Johnson", "Geneva"):
        return "$" + p.replace("m", "-") + "$"
    if system == "Stroemgren":
        if "beta" in p:
            return r"$H\beta$"
        p = p.replace("bmy", "b-y").replace("umb", "u-b").replace("1", "_1")
        return "$" + p + "$"
    if system == "Gaia":
        return {"G": "$G$", "BP": r"$G_{\mathrm{BP}}$",
                "RP": r"$G_{\mathrm{RP}}$"}.get(p, "$" + p + "$")
    if system == "Hipparcos":
        return "$" + p.replace("p", "_p") + "$"
    if system == "Tycho":
        return "$" + p.replace("t", "_T") + "$"
    if system in ("JPLUS", "SPLUS"):
        return "$" + p.replace("J0", "") + "$"
    if system in ("NSFCam", "ERIS_NIX"):
        return "$" + p.replace("Ks", "K_s") + "$"
    if system == "AstroSat":
        for a, b in (("F", ""), ("W", ""), ("M", ""), ("N2", "2"),
                     ("VIS1", "V_1"), ("VIS2", "V_2"), ("VIS3", "V_3")):
            p = p.replace(a, b)
        return "$" + p + "$"
    if system in ("SWIFT", "XMM"):
        return p.replace("UV", "") if "UV" in p else "$" + p + "$"
    if system == "Spitzer":
        return p.replace("_", ".")
    if system == "HST":
        return "".join(c for c in p if c.isdigit())
    return "$" + p + "$"


def style_for(system):
    for key, val in SYSTEM_STYLE.items():
        if system == key or (key == "DELVE" and "DELVE" in system):
            return val
    return ("black", -1)


def decode(col):
    return np.array([v.decode() if isinstance(v, bytes) else v for v in col])


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    infile, outfile = sys.argv[1], sys.argv[2]

    with fits.open(infile) as hdul:
        names = [h.name.upper() for h in hdul]
        if "SPECTRUM_FIT" not in names or "FILTERS" not in names:
            sys.exit("plot_sed.py: %s lacks spectrum_fit/filters extensions "
                     "(run sedfit with write_fits and write_model)" % infile)
        spec = hdul["SPECTRUM_FIT"].data
        filt = hdul["FILTERS"].data
        comp_cols = [c for c in spec.columns.names if c.startswith("f_c")]
        cols = hdul["COLOURS"].data if "COLOURS" in names else None
        iue = hdul["IUE"].data if "IUE" in names else None

    l = np.asarray(spec["l"], float)
    f = np.asarray(spec["f"], float)
    win = (l >= XMIN) & (l <= XMAX)
    lm, fm = l[win], f[win] * l[win] ** EXPONENT

    good = filt["flag"] == 0
    lam = np.asarray(filt["lambda"], float)[good]
    inwin = (lam >= XMIN) & (lam <= XMAX)
    lam = lam[inwin]
    lam_lo = np.asarray(filt["lambda_min"], float)[good][inwin]
    lam_hi = np.asarray(filt["lambda_max"], float)[good][inwin]
    flux = np.asarray(filt["flux"], float)[good][inwin] * lam ** EXPONENT
    flux_lo = np.asarray(filt["flux_min"], float)[good][inwin] * lam ** EXPONENT
    flux_hi = np.asarray(filt["flux_max"], float)[good][inwin] * lam ** EXPONENT
    chi = (np.asarray(filt["diff"], float) /
           np.asarray(filt["diff_err"], float))[good][inwin]
    systems = decode(filt["system"])[good][inwin]
    passbands = decode(filt["passband"])[good][inwin]

    # decimal power factored out of f*l^3 (xfig_photometry.sl:901-916)
    ymax_data = max(fm.max(), flux_hi.max())
    factor = np.floor(np.log10(ymax_data))
    scale = 10.0 ** factor
    fm, flux, flux_lo, flux_hi = (a / scale for a in (fm, flux, flux_lo,
                                                      flux_hi))
    ymin = min(fm.min(), flux_lo.min())
    ymax = ymax_data / scale
    pad = 0.1 * (ymax - ymin)
    ymin, ymax = ymin - pad, ymax + pad

    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix",
        "font.size": 9, "pdf.fonttype": 42,
    })
    have_cols = cols is not None and len(cols) > 0
    fig = plt.figure(figsize=(15 / 2.54 * (1.18 if have_cols else 1.0),
                              10 / 2.54))
    # main 0.85w x 0.75h, residual 0.85w x 0.25h, side 0.14w x 1.0h
    gs = fig.add_gridspec(2, 2 if have_cols else 1, height_ratios=[3, 1],
                          width_ratios=[0.85, 0.15] if have_cols else [1],
                          hspace=0.0, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)

    ax.set_xscale("log")
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(ymin, ymax)
    ax.plot(lm, fm, color=XFIG_COLORS["gray"], lw=0.6, zorder=2)
    # per-component curves (photometry.sl:1429-1455): c1 skyblue, c2 pink3
    for i, c in enumerate(comp_cols):
        fc = np.asarray(spec[c], float)[win] * lm ** EXPONENT / scale
        ax.plot(lm, fc, color=XFIG_COLORS["skyblue" if i == 0 else "pink3"],
                lw=0.6, zorder=1.9)
    if iue is not None:
        ax.plot(np.asarray(iue["wavelength"], float),
                np.asarray(iue["flux"], float) *
                np.asarray(iue["wavelength"], float) ** EXPONENT / scale,
                color=XFIG_COLORS["violet"], lw=0.5, zorder=2.1)

    for i in range(len(lam)):
        cname, rel_y = style_for(systems[i])
        c = XFIG_COLORS[cname]
        ax.plot([lam[i]], [flux[i]], marker="x", ms=3.5, mew=0.9, color=c,
                ls="none", zorder=4)
        ax.plot([lam[i], lam[i]], [flux_lo[i], flux_hi[i]], color=c, lw=0.8,
                zorder=3)
        # filter width at tenth maximum, dashed (filter_width qualifier)
        ax.plot([lam_lo[i], lam_hi[i]], [flux[i], flux[i]], color=c, lw=1.0,
                ls=(0, (5, 3)), zorder=3)
        label = pretty_passband(systems[i], passbands[i])
        if rel_y <= 0 and flux_hi[i] <= ymax:
            ax.annotate(label, (lam[i], flux_hi[i]), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom",
                        color=c, fontsize=7.5, zorder=5)
        elif flux_lo[i] >= ymin:
            ax.annotate(label, (lam[i], flux_lo[i]), xytext=(0, -2),
                        textcoords="offset points", ha="center", va="top",
                        color=c, fontsize=7.5, zorder=5)
    exp_txt = "" if factor == 0 else r"10^{%d}\," % int(factor)
    ax.set_ylabel(r"$f\lambda^{%d}$ ($%s$erg cm$^{-2}$ s$^{-1}$ Å$^{%d}$)"
                  % (EXPONENT, exp_txt, EXPONENT - 1))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.tick_params(which="both", direction="in", top=True, right=True)

    # chi_magnitude pane (xfig_photometry.sl:976-1046)
    cmax = max(chi.max() + 0.5, 2.25)
    cmin = min(chi.min() - 0.5, -2.25)
    cpad = 0.04 * (cmax - cmin)
    axr.set_ylim(cmin - cpad, cmax + cpad)
    axr.axhline(0, color="#808080", lw=0.7, ls=(0, (5, 3)), zorder=1)
    axr.axhline(1, color="#808080", lw=0.7, ls=(0, (1, 2)), zorder=1)
    axr.axhline(-1, color="#808080", lw=0.7, ls=(0, (1, 2)), zorder=1)
    for i in range(len(lam)):
        cname, _ = style_for(systems[i])
        axr.plot([lam[i]], [chi[i]], marker="x", ms=3.5, mew=0.9,
                 color=XFIG_COLORS[cname], ls="none", zorder=3)
    axr.set_xlabel(r"$\lambda$ (Å)")
    axr.set_ylabel(r"$\chi_{\mathrm{magnitude}}$")
    axr.tick_params(which="both", direction="in", top=True, right=True)
    ticks = [t for t in (1500, 2000, 3000, 5000, 8000, 12000, 20000, 30000,
                         50000) if XMIN <= t <= XMAX]
    axr.set_xticks(ticks)
    axr.set_xticklabels(["%d" % t for t in ticks])
    axr.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    # chi_colour side pane (xfig_photometry.sl:1050-1122)
    if have_cols:
        axs = fig.add_subplot(gs[:, 1])
        cchi = np.asarray(cols["diff"], float) / np.asarray(cols["diff_err"],
                                                            float)
        csys = decode(cols["system"])
        cpass = decode(cols["passband"])
        n = len(cchi)
        xhi = max(cchi.max() + 0.5, 2.25)
        xlo = min(cchi.min() - 0.5, -2.25)
        xpad = 0.04 * (xhi - xlo)
        axs.set_xlim(xlo - xpad, xhi + xpad)
        axs.set_ylim(0, n + 1)
        axs.axvline(0, color="#808080", lw=0.7, ls=(0, (5, 3)), zorder=1)
        axs.axvline(1, color="#808080", lw=0.7, ls=(0, (1, 2)), zorder=1)
        axs.axvline(-1, color="#808080", lw=0.7, ls=(0, (1, 2)), zorder=1)
        for i in range(n):
            cname, _ = style_for(csys[i])
            axs.plot([cchi[i]], [i + 1], marker="x", ms=3.5, mew=0.9,
                     color=XFIG_COLORS[cname], ls="none", zorder=3)
        axs.set_yticks(range(1, n + 1))
        axs.set_yticklabels([pretty_passband(csys[i], cpass[i])
                             for i in range(n)])
        axs.yaxis.tick_right()
        axs.set_xlabel(r"$\chi_{\mathrm{color}}$")
        axs.tick_params(which="both", direction="in", top=True)

    fig.savefig(outfile, bbox_inches="tight")
    print("- wrote %s" % outfile)


if __name__ == "__main__":
    main()
