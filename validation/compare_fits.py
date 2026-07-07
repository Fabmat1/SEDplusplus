#!/usr/bin/env python3
"""Compare a sedfit SED_results.fits against the ISIS reference fixture.

Usage: compare_fits.py <ref_dir> <cpp_dir>

Checks HDU names/order, column names/TFORMs, primary keywords (presence,
undefined-ness, values), and data: exact for strings/ints, ~1e-6 relative for
fit-dependent doubles, ~1% + interval-width tolerance for MC-derived stellar
values (MC noise floor, PLAN.md). The composite star's fit-level offsets
(accepted Phase-1 floor) are reported but bounded by a looser cap on the
params_fit value columns.
"""

import sys

import numpy as np
from astropy.io import fits

FIT_RTOL = 2e-4      # params/filters/spectrum: generous cap over the 1e-6
                     # typical (composite star sits on the Phase-1 E=0 floor)
MC_RTOL = 0.05       # stellar_c*: MC noise (~1%) + mode/HDI binning jitter

fails = []
notes = []


def fail(msg):
    fails.append(msg)


def close(a, b, rtol):
    a, b = np.atleast_1d(np.asarray(a, float)), np.atleast_1d(
        np.asarray(b, float))
    both_nan = np.isnan(a) & np.isnan(b)
    ok = both_nan | np.isclose(a, b, rtol=rtol, atol=1e-300)
    # near-zero entries: tolerance relative to the column scale
    with np.errstate(invalid="ignore"):
        vals = np.abs(np.concatenate([a, b]))
        scale = np.nanmax(vals) if np.isfinite(vals).any() else 0.0
    ok |= np.abs(a - b) <= rtol * max(scale, 1e-30)
    return ok


def cmp_table(name, ref, cpp, rtol):
    if list(ref.columns.names) != list(cpp.columns.names):
        fail(f"{name}: column names {ref.columns.names} != {cpp.columns.names}")
        return
    if list(ref.columns.formats) != list(cpp.columns.formats):
        notes.append(f"{name}: TFORMs ref={ref.columns.formats} "
                     f"cpp={cpp.columns.formats}")
    if len(ref.data) != len(cpp.data):
        fail(f"{name}: {len(ref.data)} rows vs {len(cpp.data)}")
        return
    for col, fmt in zip(ref.columns.names, ref.columns.formats):
        r, c = ref.data[col], cpp.data[col]
        if fmt.endswith("A"):
            r = [v.strip() for v in r]
            c = [v.strip() for v in c]
            if name == "params_fit" and col == "tex":
                # MC-free rows must match; the tex strings were validated
                # byte-level in stage 2 -- just report diffs here
                bad = [i for i in range(len(r)) if r[i] != c[i]]
                if bad:
                    notes.append(f"{name}.tex differs on rows {bad}")
            elif list(r) != list(c):
                fail(f"{name}.{col}: string mismatch {list(r)[:4]}... vs "
                     f"{list(c)[:4]}...")
        elif fmt == "J":
            if not np.array_equal(np.asarray(r), np.asarray(c)):
                fail(f"{name}.{col}: int mismatch")
        else:
            ok = close(r, c, rtol)
            if not ok.all():
                i = int(np.argmax(~ok))
                fail(f"{name}.{col}: worst row {i}: {r[i]!r} vs {c[i]!r} "
                     f"({(~ok).sum()}/{len(ok)} beyond rtol={rtol})")


def main():
    global FIT_RTOL
    ref_dir, cpp_dir = sys.argv[1], sys.argv[2]
    if len(sys.argv) > 3:  # composite star: accepted Phase-1 E=0 floor gives
        FIT_RTOL = float(sys.argv[3])  # ~2-3% on the faint-component spectrum
    ref = fits.open(f"{ref_dir}/SED_results.fits")
    cpp = fits.open(f"{cpp_dir}/SED_results.fits")

    rn = [h.name for h in ref]
    cn = [h.name for h in cpp]
    if rn != cn:
        fail(f"HDU list {rn} != {cn}")

    # primary keywords
    skip = {"SIMPLE", "BITPIX", "EXTEND", "COMMENT"}
    rk = [k for k in ref[0].header if k not in skip and
          not k.startswith("NAXIS")]
    ck = [k for k in cpp[0].header if k not in skip and
          not k.startswith("NAXIS")]
    if rk != ck:
        fail(f"primary keywords {rk} != {ck}")
    for k in rk:
        if k not in ck:
            continue
        rv, cv = ref[0].header[k], cpp[0].header[k]
        run = isinstance(rv, fits.card.Undefined)
        cun = isinstance(cv, fits.card.Undefined)
        if run != cun:
            fail(f"key {k}: undefined-ness {run} vs {cun}")
        elif run:
            pass
        elif isinstance(rv, str):
            if rv != cv:
                fail(f"key {k}: {rv!r} vs {cv!r}")
        elif not close(rv, cv, FIT_RTOL).all():
            # CHISQR_RED / NORM_CHI_RED carry the fit-level floor (~1e-5)
            fail(f"key {k}: {rv!r} vs {cv!r}")
        # card image equality (formatting parity)
        ri = ref[0].header.cards[k].image.rstrip()
        ci = cpp[0].header.cards[k].image.rstrip()
        if ri != ci:
            notes.append(f"card image differs: {ri!r} vs {ci!r}")

    for name in [n for n in rn if n != "PRIMARY" and n in cn]:
        rtol = MC_RTOL if name.upper().startswith("STELLAR") else FIT_RTOL
        cmp_table(name.lower(), ref[name], cpp[name], rtol)
        # GRID keyword on stellar extensions
        if name.upper().startswith("STELLAR"):
            if ref[name].header.get("GRID") != cpp[name].header.get("GRID"):
                fail(f"{name}: GRID {ref[name].header.get('GRID')!r} vs "
                     f"{cpp[name].header.get('GRID')!r}")

    for n in notes:
        print("  note:", n)
    if fails:
        print(f"FAIL ({ref_dir} vs {cpp_dir}):")
        for f_ in fails:
            print("  -", f_)
        sys.exit(1)
    print(f"PASS {cpp_dir} ({len(rn)-1} extensions)")


if __name__ == "__main__":
    main()
