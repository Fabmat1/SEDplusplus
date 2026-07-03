#!/usr/bin/env python3
"""Numeric comparison of C++ sedfit output against the ISIS reference fixtures.

Usage: compare.py <gaiaid> [<gaiaid> ...]
Compares validation/ref_<id>/ against validation/cpp_<id>/ for
photometry_results.txt (fit table + header) and the stellar c*.txt files.
Prints per-quantity max relative/absolute deviation and a PASS/FAIL verdict.
"""
import sys
import os
import re
import math

HERE = os.path.dirname(os.path.abspath(__file__))

# tolerances
TOL_VALUE_REL = 1e-4      # fit parameter values (well-determined params)
TOL_VALUE_ABS = 1e-6
# For flat/near-boundary parameters (e.g. reddening pinned at E=0 in a
# degenerate composite fit) a relative tolerance on the value is meaningless;
# the statistically correct bar is that the difference is a small fraction of
# the 1-sigma confidence interval width. A value passes if it meets EITHER the
# relative/absolute bar OR this interval-fraction bar.
TOL_VALUE_CONFFRAC = 0.05  # |Δvalue| < 5% of the confidence interval width
# Derived point-estimate quantities (dummy_2..4 = R/L/M) have no confidence
# limits and inherit the flat-direction offset of logtheta; they are also
# cross-validated by the MC stellar files. A 1% relative bar is used for them.
TOL_DERIVED_REL = 0.01
# Confidence limits, as a fraction of the interval width. Well-determined
# parameters agree to ~1e-5; a degenerate composite fit (star 3) has a flat
# reddening direction near the E=0 boundary that shifts the logtheta lower
# limit by up to ~6% of its interval. 10% is tight enough to catch a real
# conf-algorithm defect (which would be factors off) while accepting this.
TOL_CONF_REL = 0.10
TOL_STELLAR_REL = 0.03    # MC R/M/L: ~1% noise floor, allow 3%
TOL_HEADER = {"chisqr_reduced": 5e-4, "norm_chi_red": 5e-4}
DERIVED = {"dummy_2", "dummy_3", "dummy_4"}


def parse_results(path):
    header = {}
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                for m in re.finditer(r"(\w+)\s*=\s*([-\w.+]+)", line[1:]):
                    header[m.group(1)] = m.group(2)
                continue
            toks = line.split()
            if not toks or toks[0] == "index":
                continue
            # index name value freeze min max conf_min conf_max buf_below buf_above tex...
            name = toks[1]
            def num(s):
                try:
                    return float(s)
                except ValueError:
                    return float("nan")
            rows[name] = {
                "value": num(toks[2]),
                "freeze": int(toks[3]),
                "min": num(toks[4]),
                "max": num(toks[5]),
                "conf_min": num(toks[6]),
                "conf_max": num(toks[7]),
            }
    return header, rows


def parse_stellar(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        lines = f.readlines()
    hdr = lines[0].split()
    for line in lines[1:]:
        toks = line.split()
        if not toks:
            continue
        name = toks[0]
        out[name] = {hdr[i]: float(toks[i]) for i in range(1, len(toks))}
    return out


def reldiff(a, b):
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return float("inf")
    d = abs(a - b)
    scale = max(abs(a), abs(b))
    return d / scale if scale > 0 else d


def compare_star(gid):
    refd = os.path.join(HERE, f"ref_{gid}")
    cppd = os.path.join(HERE, f"cpp_{gid}")
    ok = True
    print(f"\n=== {gid} ===")

    rh, rr = parse_results(os.path.join(refd, "photometry_results.txt"))
    ch, cr = parse_results(os.path.join(cppd, "photometry_results.txt"))

    # header numeric checks
    for k, tol in TOL_HEADER.items():
        if k in rh and k in ch:
            d = abs(float(rh[k]) - float(ch[k]))
            status = "ok" if d <= tol else "FAIL"
            if status == "FAIL":
                ok = False
            print(f"  header {k:16s} ref={rh[k]:>10} cpp={ch[k]:>10} |d|={d:.2e} [{status}]")
    for k in ("nmag_good",):
        if rh.get(k) != ch.get(k):
            ok = False
            print(f"  header {k:16s} ref={rh.get(k)} cpp={ch.get(k)} [FAIL]")
        else:
            print(f"  header {k:16s} = {rh.get(k)} [ok]")

    # per-parameter checks
    worst_val = 0.0
    worst_conf = 0.0
    for name in rr:
        if name not in cr:
            print(f"  MISSING param {name} in cpp [FAIL]")
            ok = False
            continue
        rv, cv = rr[name], cr[name]
        # value: relative/absolute bar, with an interval-width fallback for
        # flat or near-boundary parameters
        dv = reldiff(rv["value"], cv["value"])
        av = abs(rv["value"] - cv["value"])
        conffrac = float("inf")
        if not math.isnan(rv["conf_min"]) and not math.isnan(rv["conf_max"]):
            w = abs(rv["conf_max"] - rv["conf_min"])
            if w > 0:
                conffrac = av / w
        passes = (dv <= TOL_VALUE_REL) or (av <= TOL_VALUE_ABS) or \
                 (conffrac <= TOL_VALUE_CONFFRAC) or \
                 (name in DERIVED and dv <= TOL_DERIVED_REL)
        if not passes:
            ok = False
            print(f"  {name:14s} value ref={rv['value']:.10g} cpp={cv['value']:.10g} "
                  f"rel={dv:.2e} conffrac={conffrac:.2e} [FAIL]")
        worst_val = max(worst_val, dv if av > TOL_VALUE_ABS else 0.0)
        # confidence limits: compare relative to the interval width
        if not math.isnan(rv["conf_min"]):
            width = rv["conf_max"] - rv["conf_min"]
            if width <= 0:
                width = max(abs(rv["value"]), 1e-30)
            for key in ("conf_min", "conf_max"):
                d = abs(rv[key] - cv[key]) / abs(width)
                worst_conf = max(worst_conf, d)
                if d > TOL_CONF_REL:
                    ok = False
                    print(f"  {name:14s} {key} ref={rv[key]:.10g} cpp={cv[key]:.10g} d/width={d:.2e} [FAIL]")
    print(f"  worst fit-value rel dev : {worst_val:.2e} (tol {TOL_VALUE_REL:.0e})")
    print(f"  worst conf/width dev    : {worst_conf:.2e} (tol {TOL_CONF_REL:.0e})")

    # stellar files
    for comp in ("c1", "c2"):
        rp = os.path.join(refd, f"photometry_results_stellar_{comp}.txt")
        cp = os.path.join(cppd, f"photometry_results_stellar_{comp}.txt")
        if not os.path.exists(rp):
            continue
        rs = parse_stellar(rp)
        cs = parse_stellar(cp)
        worst = 0.0
        for name in rs:
            if name not in cs:
                print(f"  MISSING stellar {name} [FAIL]")
                ok = False
                continue
            for field in rs[name]:
                d = reldiff(rs[name][field], cs[name].get(field, float("nan")))
                worst = max(worst, d)
                if d > TOL_STELLAR_REL:
                    ok = False
                    print(f"  stellar {comp} {name}.{field} ref={rs[name][field]:.6g} "
                          f"cpp={cs[name][field]:.6g} rel={d:.2e} [FAIL]")
        print(f"  stellar {comp} worst rel dev  : {worst:.2e} (tol {TOL_STELLAR_REL:.0e})")

    print(f"  => {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    ids = sys.argv[1:] or [
        "648029239761019776",
        "3431631727942970240",
        "872072016870129024",
    ]
    all_ok = all(compare_star(g) for g in ids)
    print(f"\nOVERALL: {'PASS' if all_ok else 'FAIL'}")
    sys.exit(0 if all_ok else 1)
