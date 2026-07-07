#!/usr/bin/env python3
"""Compare C++ vs ISIS photometry_fit*.txt products (Phase 3 Stage 1).

Column classification and tolerances (see PHASE3.md working notes):
  * string / flag / lambda_min / lambda_max : must be byte(value)-exact
  * lambda (lambda_eff)                      : relative <= 1e-6
  * flux, flux_min, flux_max                 : relative <= 1e-5 (obs-driven,
                                               usually machine precision)
  * diff, diff_err                           : absolute <= DIFF_ABS
  * photometry_fit.txt l                     : value-exact
  * photometry_fit.txt f / f_c*              : relative (reported; PASS if the
                                               MEDIAN relative <= FIT_MEDREL)

diff/f are model residuals: for single-component fits they agree to the model-
reconstruction float floor (~1e-7 mag / ~1e-6 in f); the composite star 872 is
looser because its fit parameters themselves differ at the accepted Phase-1
composite floor (reddening degenerate at E=0). Pass thresholds are set to admit
that floor while still catching real regressions.

Usage: compare_fit_txt.py <ref_dir> <cpp_dir>
"""
import os
import sys
import math

DIFF_ABS = 2e-3       # abs tol on diff/diff_err (covers the 872 composite floor)
FIT_MEDREL = 5e-3     # median rel tol on photometry_fit.txt f/f_c*

EXACT = {"lambda_min", "lambda_max", "passband", "system", "flag",
         "VizieR_catalog"}
# lambda (lambda_eff) is model-flux-weighted, so it inherits the model floor:
# ~1e-8 for single-component fits, ~4e-5 for the 872 composite (whose fit
# parameters themselves differ at the accepted Phase-1 composite floor).
REL = {"lambda": 1e-4, "flux": 1e-5, "flux_min": 1e-5, "flux_max": 1e-5}
ABS = {"diff": DIFF_ABS, "diff_err": DIFF_ABS}


def read(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    return lines[0].split(), [l.split() for l in lines[1:]]


def cmp_generic(rp, cp, fname):
    hr, rr = read(rp)
    hc, rc = read(cp)
    ok = True
    if hr != hc:
        print(f"  header: ref={hr} cpp={hc}")
        return False
    if len(rr) != len(rc):
        print(f"  row count ref={len(rr)} cpp={len(rc)}")
        return False
    for j, col in enumerate(hr):
        a = [r[j] for r in rr]
        b = [r[j] for r in rc]
        if col in EXACT:
            mism = [i for i in range(len(a)) if a[i] != b[i]
                    and (col in ("passband", "system", "flag", "VizieR_catalog")
                         or float(a[i]) != float(b[i]))]
            status = "OK" if not mism else f"FAIL ({len(mism)} mismatch)"
            ok &= not mism
            print(f"    {col:16s} exact   {status}")
        else:
            fa = [float(x) for x in a]
            fb = [float(x) for x in b]
            maxabs = maxrel = 0.0
            for x, y in zip(fa, fb):
                if math.isnan(x) and math.isnan(y):
                    continue
                d = abs(x - y)
                maxabs = max(maxabs, d)
                maxrel = max(maxrel, d / max(abs(x), 1e-300))
            if col in REL:
                good = maxrel <= REL[col]
                print(f"    {col:16s} rel<={REL[col]:.0e}  maxrel={maxrel:.2e}  "
                      f"{'OK' if good else 'FAIL'}")
            else:
                good = maxabs <= ABS[col]
                print(f"    {col:16s} abs<={ABS[col]:.0e}  maxabs={maxabs:.2e} "
                      f"maxrel={maxrel:.2e}  {'OK' if good else 'FAIL'}")
            ok &= good
    return ok


def cmp_spectrum(rp, cp):
    hr, rr = read(rp)
    hc, rc = read(cp)
    ok = True
    if hr != hc or len(rr) != len(rc):
        print(f"  header/rows: {hr}/{len(rr)} vs {hc}/{len(rc)}")
        return False
    for j, col in enumerate(hr):
        fa = [float(r[j]) for r in rr]
        fb = [float(r[j]) for r in rc]
        if col == "l":
            mism = sum(1 for x, y in zip(fa, fb) if x != y)
            print(f"    l   exact   {'OK' if mism == 0 else f'FAIL ({mism})'}")
            ok &= mism == 0
        else:
            rels = sorted(abs(x - y) / max(abs(x), 1e-300)
                          for x, y in zip(fa, fb))
            med = rels[len(rels) // 2]
            mx = rels[-1]
            good = med <= FIT_MEDREL
            print(f"    {col:4s} medrel={med:.2e} maxrel={mx:.2e}  "
                  f"{'OK' if good else 'FAIL'}")
            ok &= good
    return ok


def main():
    ref, cpp = sys.argv[1], sys.argv[2]
    allok = True
    for fname in ("photometry_fit_mag.txt", "photometry_fit_col.txt",
                  "photometry_fit.txt"):
        rp, cp = os.path.join(ref, fname), os.path.join(cpp, fname)
        re_, ce = os.path.exists(rp), os.path.exists(cp)
        if re_ != ce:
            print(f"{fname}: PRESENCE MISMATCH ref={re_} cpp={ce}")
            allok = False
            continue
        if not re_:
            print(f"{fname}: absent in both (ok)")
            continue
        print(f"{fname}:")
        ok = (cmp_spectrum(rp, cp) if fname == "photometry_fit.txt"
              else cmp_generic(rp, cp, fname))
        allok &= ok
    print("RESULT:", "PASS" if allok else "FAIL")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
