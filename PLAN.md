# SED++ — C++ port of the ISIS photometry SED-fitting pipeline

**Handoff document.** Everything needed to finish the port is written down here,
including parity-critical details already extracted from the S-Lang/ISIS sources.
Read this fully before touching code.

## STATUS: Phase 1 COMPLETE — validated to parity (2026-07-03)

`sedfit` builds and passes `validation/compare.py` for all three reference
stars. Single-component fits match to ~1e-7; the composite fit matches to
≤1.3% of each parameter's 1σ interval (its reddening is degenerate at the E=0
boundary). See README.md. What remains is phase 2 (step 5 below): remote/SQL
querying, bulk + parallel multi-star fitting, plots/TeX/FITS outputs, and MC
speedup (the 2e6-sample R/M/L Monte-Carlo dominates the ~5–11s runtime; the fit
itself is <1s).

### Key discoveries made during validation (do not re-derive)

1. **`set_par` value rejection** (ISIS `do_set_par`): min/max/freeze are applied
   first; the value is assigned only if inside [min,max], else rejected (kept at
   default) with a warning — the min/max/freeze changes persist. This is why
   star 3's `c1_HE=-5.05` (with min=-5) stays frozen at the grid-mid default
   -2.546. Implemented in `FitFunction::set_par` (all three overloads).
2. **Confidence-loop restart** (ISIS `conf_loop.sl` `recv_slave_result`): when a
   single-parameter conf search re-fits the others and finds a lower statistic
   (strict `<`, `max_param_retries` defaults to 0 so collapsed-limit restarts
   are off), the improved fit is *adopted* and the whole conf loop restarts
   (param_version bump discards prior limits). This lets a fit that stopped
   short on a flat/degenerate surface settle into the true minimum before clean
   limits are computed. Implemented as the `conf_all` restart loop in
   `Fitter::run`. Without it the composite fit's conf limits collapse
   (conf_min==conf_max) and R/M/L are ~60% off.
3. **MC deltas use literal names only**: the template's
   `where(par_full.name==sprintf("c%d_%s",i,field))` does NOT match wildcard
   par_full entries (`c*_teff`), so wildcard-prescribed parameters contribute
   ZERO spread to the R/M/L Monte-Carlo. `main.cpp`'s `delta` lambda matches
   this (literal `c{i}_field` only).
4. **S-Lang `%S` float format** (`results.cpp` `slang_repr`): shortest
   round-trip digits, `%g`-style selection (fixed for decimal exponent in
   [-4,5], scientific otherwise), decimal point always present in fixed form,
   `-0.0`/`-nan` preserved.
5. **Residual composite-fit floor**: on the flat E_44m55 direction the C++ and
   ISIS mpfit fixed points differ by ~2.8e-4 in χ² (4 sig figs agree), giving
   ≤1.3%-of-1σ parameter offsets. This is irreducible without bit-identical
   model floats and is physically negligible (E consistent with 0).

## Goal

Rewrite `~/Projects/ISIS_install/src/stellar_isisscripts/templates/photometry.sl`
in C++ ("sedfit"), reaching output parity with the old script for three test
stars, then (step 2, later) add a bulk-safe query layer + local SQL source and
parallel multi-star fitting at 100s of fits/second.

**Phase 1 scope (current):** fit engine only. Inputs (`photometry.dat`,
astrometry from the `photometry_results.txt` header) are read from cached files
produced by the old pipeline. No remote queries, no plots, no TeX/PDF, no
SED_results.fits. Outputs: `photometry_results.txt` +
`photometry_results_stellar_c*.txt`.

## Current state (what is DONE)

- **Validation fixtures** in `validation/ref_<gaiaid>/` for the three stars
  (Gaia DR3 648029239761019776, 3431631727942970240, 872072016870129024 —
  the last one is a composite sdB+late-type fit with two grids).
  Each dir contains: `photometry.sl` (the exact script run), `photometry.dat`
  (cached photometry), `photometry_results.txt` + `photometry_results_stellar_c*.txt`
  (= reference outputs, produced by a **second** run of the old script so they
  are self-consistent with the cached inputs incl. the %.3g-rounded parallax
  in the header), `run.log`, `run2.log`, and `*_run1.txt` copies of the
  first-run outputs.
  - Run the old script via `cd <dir> && ~/Projects/ISIS_install/bin/isis-script photometry.sl`
    (needs network only if `photometry.dat` is missing).
  - **Noise floor between runs** (run1 vs run2, identical inputs except
    full-precision vs %.3g parallax): fit values + conf limits bit-identical;
    MC-derived R/M/L differ by ~0.1–0.8% (MC noise + parallax rounding).
- **All C++ sources written** (not yet compiled!) in `src/`:
  `util`, `extinction`, `fitsio_util`, `modelgrid`, `passbands`, `synthmag`,
  `fitfun`, `fitter`, `photometry_table`, `stellar`, `texval`, `results`,
  `config`, `main.cpp`; vendored `extern/mpfit/{mpfit.c,mpfit.h}` (copied from
  ISIS source, same LM engine); auto-generated `src/zp_table.inc` (289
  zero-point rows) and `src/f19_table.inc` (Fitzpatrick+2019 tables), both
  extracted programmatically from the .sl sources (do NOT hand-edit).
- `CMakeLists.txt` (C++20, cfitsio + gsl + nlohmann_json via pkg-config,
  static mpfit lib). All deps are installed on this machine.

## Remaining steps

### 1. Build (`task #2/#3/#4` wrap-up)
```
cd ~/Projects/SEDplusplus && cmake -B build && cmake --build build -j
```
Expect minor compile fixes (the code was written in one pass). Known loose ends:
- `src/fitter.hpp/cpp`: remove the dead `newton_norm` helper and the unused
  `friend class ConfSearch`.
- clangd flagged unused includes in several files — harmless.
- `mp_par`/`mp_config` come from `extern/mpfit/mpfit.h` (CMake adds the include
  dir via the static lib's PUBLIC include).

### 2. Write the three JSON configs (`configs/<gaiaid>.json`)
Mirror the user-supplied headers (they are embedded at the top of each
`validation/ref_*/photometry.sl`). Schema is `src/config.cpp`. Example:
```json
{
  "star": "Gaia DR3 648029239761019776",
  "coordinates": {"ra": 145.977, "dec": 27.7829},
  "par":      {"name": ["c*_xi","c*_z"], "value": [0,0], "freeze": [1,1]},
  "par_full": {"name": ["c*_HE","c*_logg","c*_teff","R_55"],
               "value": [-3.7499743,5.5231117,29232.555,3.02],
               "freeze": [1,1,1,1],
               "min": [-3.8441735,5.4471857,28811.428,2.5],
               "max": [-3.6557751,5.5990378,29653.683,6]},
  "griddirectories": ["processed/"],
  "bpaths": ["./","/home/fabian/ISIS_models/sdB/","/home/fabian/ISIS_models/"],
  "conf_level": 0, "apply_ZPO_corr": 1, "remove_outliers": 5, "nMC": 2000000,
  "refdata": "/home/fabian/Projects/ISIS_install/src/stellar_isisscripts/refdata",
  "workdir": "validation/ref_648029239761019776",
  "outdir":  "validation/cpp_648029239761019776"
}
```
Star 2: ra=90.1291 dec=29.1486, par_full values from its header.
Star 3 (composite): `griddirectories: ["processed/", "/home/fabian/ISIS_models/Phoenix_late_type_stars_photometry_v2.0/processed/"]`,
`par` freezes c1_xi,c1_z,c2_xi,c2_z,c2_HE(=-1.05) and leaves c2_logg(=4.0),
c2_teff(=6000) thawed; par_full prescribes c1_HE,c1_logg,c1_teff,R_55.
**outdir must differ from workdir** or the C++ run overwrites the fixture.

### 3. Validate (task #5)
Run `./build/sedfit configs/<id>.json` for each star and compare
`validation/cpp_<id>/photometry_results.txt` against
`validation/ref_<id>/photometry_results.txt`:
- **Fit parameter values** (logtheta, E_44m55, c2_teff, c2_logg, c2_sur_ratio):
  should agree to ~1e-6 relative (same minimizer, same numerics); accept up to
  ~1e-4 relative or ≪ conf interval width.
- **conf_min/conf_max**: agree to ~1e-3 of the interval width (conf tolerance
  is 1e-3·Δχ²); logtheta/c1_teff conf after error propagation similar.
- **chisqr_reduced (header, %.3f)** and **norm_chi_red (%.4f)**: match printed
  digits (star2 has norm_chi_red=0.0087 — exercises the excess-noise path).
- **nmag_good**: exact match (16 / 14 / 28).
- **Stellar R/M/L** (mode/HDI/median/quantiles): within ~1% (MC noise floor,
  see above). dummy_2..4 in the results table: should match ref run2 values
  (which used the %.3g-rounded parallax) to ~1e-6.
- Debug order if values disagree: (a) synthetic magnitudes for the best-fit
  reference parameters vs `photometry_fit_mag.txt` in the ref dirs (column
  `magnitude_model` is printed by the old script — compare per-filter; this
  isolates grid interpolation / extinction / filter integration from the fit),
  (b) chi² at reference best-fit params, (c) the fit itself, (d) conf limits,
  (e) propagation, (f) MC.

### 4. Cleanup + docs
- README.md: build, config schema, what is/isn't implemented.
- Commit-worthy state; user may want git init (ask).

### 5. Step 2 (later, separate effort — design already agreed with user)
- Query layer in C++ (libcurl): VizieR/Gaia TAP/IRSA with request throttling +
  bulk TAP uploads (avoid per-star hammering), optional local SQL (sqlite/
  postgres) photometry source; multi-star batch mode with OpenMP; keep
  fit engine unchanged. The S-Lang reference for catalog logic is
  `src/photometry/query_photometry.sl` (163 KB — port incrementally).

## Parity-critical implementation knowledge (already baked into the C++ code —
## verify against this if debugging)

### Pipeline flow (photometric_fitting + template)
1. Read astrometry from `photometry_results.txt` header
   (`# key = value; ...`, needs parallax/parallax_error/RA/DEC). Parallax cuts:
   NaN unless `parallax>0` and `parallax/parallax_error>=1`.
2. Read `photometry.dat` (fixed 8-column whitespace table; header lines
   `# RA = %.10f DEC = %.10f`, `# meanSFD = ... stdSFD = ...`, etc.).
   Entries with flag==0 and uncertainty<1e-5 or NaN get uncertainty=0.025.
3. E(44-55) initial guess: E3d=meanStilism (absent here), E2d=meanSandF else
   meanSFD; outside_stilism if parallax/err>2 && parallax<=0.5;
   E = (E3d if !outside else E2d)/0.91, capped at 2; **prepended** to par_full
   as free param, range [0,5] — only if par_full has no E_44m55 entry.
4. Parameter table (order matters, 1-based indices in output):
   logtheta(-9, free, [-18,-2], relstep 1e-4), E_44m55(0, free, [0,10] default,
   relstep 1e-3), R_55(3.02, free, [2.5,6], relstep 1e-3); per component:
   vsini/zeta/vrad (frozen), teff/logg/xi/z/HE (mid-of-grid default value
   rounded through %.0f/%.3f/%.2f/%.2f/%.3f; frozen iff grid dimension has a
   single node; min/max = grid extremes; hard limits extended by one node
   spacing, teff/logg/xi clamped ≥0), sur_ratio (comp1: frozen 1 [1,1];
   others: free 1 [0,1500]); bb/bb2 teff+sur_ratio (frozen 0); dummy_1..4 when
   good_astrometry (dummy_1 = parallax, [0,200]; dummy_2..4 are derived params
   → min/max ±DBL_MAX, freeze 1, excluded from fitting):
   dummy_2 = 10^logtheta/(2e-3·parallax)·4.4353565926280975e7 (R in Rsun),
   dummy_3 = R²(c1_teff/5772)⁴, dummy_4 = 10^c1_logg·R²·3.6469715273112305e-5.
5. `set_par` semantics: wildcard `*` matching; validation **errors** (caught →
   warning, param unchanged), no clamping; min==max==0 means "use hard
   limits"; freeze<0 = keep.
   Order: `par` (name/value/freeze) first, then `par_full` (+min/max).
6. Data selection: an observed entry is fittable if flag==0 and its
   system+passband has a non-NaN synthetic value. mag_ind = union of needed
   magnitude rows (colors pull in their dependencies).
7. Initial guesses (only if not frozen and not named in par/par_full):
   E from median of (obs−synth) B−V/b−y color offsets (max(0,·));
   logtheta = median over magnitude-type entries of
   log10(2) + 0.5·(−0.4·(m_obs − m_synth_noTheta)).
8. Fit: mpfit (vendored, identical to ISIS "mpfit" default engine) with
   ftol=xtol=gtol=1e-10, stepfactor=100, epsfcn=DBL_EPSILON, maxiter=200;
   bounded params; residual (data−model)/err; err = sqrt(unc² + ZP_err² + norm²).
9. Outlier removal (remove_outliers=5): loop {refit; worst |chi|; remove
   (flag=-2) if >5σ and not the last magnitude; len>nfree+1 required}; final refit.
10. norm_chi_red: if red-χ²>1, Newton iteration on norm so red-χ²=1
    (initial guess = avg of the two analytic roots; refit every iteration;
    threshold 2e-6, cap 100 iterations; if new<thresh → new=norm/2).
11. Confidence limits (conf_level 0/1/2 → Δχ² = 1.00/2.71/6.63, tol 1e-3):
    exact port of ISIS `fit-chisqrconf.c` (bracket from 0.5·(pbest+bound),
    halving toward the bound, then quadratic→linear→bisect interpolation with
    the sub-fit re-minimizing all other free params at each trial value;
    accept_tol = tol·Δχ²; improved-fit detection at 0.5·tol·Δχ²).
    If norm≠0: after conf of all params, if red-χ² < 0.99·old → re-normalize
    (Newton, threshold 1e-6, start=current norm) and repeat conf loop.
12. Prescribed-parameter error propagation into logtheta and c1_teff
    (runs only if logtheta has conf limits): for comp i and fields
    teff/logg/xi/z/HE/sur_ratio (+E_44m55 for i==1, using the **last** par_full
    E entry): if par_full has an entry with the **literal** name `c{i}_{field}`
    (wildcards like `c*_teff` do NOT match!) and freeze==1 → set param to
    par_full.min, refit, accumulate logtheta/c1_teff shifts in quadrature on
    the corresponding side; same for max; restore best params. Update
    conf_min/max, buf_*, tex (sci=6) of logtheta always and of c1_teff (tex
    only if free). This is why star1/2 (c*_ wildcards) show NO propagation and
    star3 (explicit c1_ names) DOES propagate.
13. MC stellar params per component (nMC=2e6): asymmetric-Gaussian draws of
    logtheta (deltas from final conf), teff/logg/sur_ratio (from conf if free;
    else from param min/max if prescribed frozen via par_full — the template
    checks par_full membership by explicit name; C++ additionally accepts the
    wildcard spelling, which is *needed* for star1/2 where c*_teff etc. are
    prescribed — this matches the old script's behaviour because there
    p.min/p.max carry the prescribed range and `temp=where(par_full.name==...)`
    ... **CHECK**: in the template, for star1/2 `where(par_full.name=="c1_teff")`
    is empty, so delta_teff falls back to 0 only if p.freeze!=0 branch...
    Actually the template's `else if(length(temp)==1 && par_full.freeze[temp[0]]==1)`
    fails for wildcard names → delta_teff_plus/minus stay 0 for star1/2!
    Compare ref stellar file spreads before trusting either variant; if the
    reference c1_R spread matches parallax+logtheta errors only, remove the
    wildcard acceptance in main.cpp (`prescribed()` lambda).
    Filter parallax>1e-4 && sur>0, positivity per quantity; mode_and_HDI with
    p=0.68268 (conf 0), 0.9/0.99 for 1/2 with parallax_error scaled by
    1.645/2.576. R = 10^logtheta/(2e-3·parallax)·pc/Rsun·sqrt(sur);
    M = 10^logg·R²·Rsun²/GMsun; L = R²(teff/5772)⁴.
    Constants: pc=3.0856775814913674e18 cm, Rsun=6.957e10 cm,
    GMsun=1.3271244e26 cgs, Teffsun=5772.
14. mode_and_HDI: 10 iterations; histogram with min(202,n) initial bin edges
    over [min−ε, max+ε]; iterations 1..8 re-bin around HDI±0.5·max(HDI width,
    quantile width) with min(200,n) points unioned with the data extremes;
    drop first/last (overflow) bins; final iteration: Savitzky-Golay smoothing
    of the **nonzero-compressed** histogram (nl=nr=nint(0.75·min((mode−lo)/step,
    (hi−mode)/step)), p=min(4,nl), only if count(nonzero)>2·nl), then linear
    interpolation onto 100× finer grid (÷100); HDI via bisection on the height
    threshold (sum of bins ≥ p). median = ISIS median (lower-middle element);
    quantile(p) = sorted[int(p·(n−1))] (no interpolation).
    If |p_a−p_s|>0.05 → mode/HDI = NaN (keep median/quantiles).

### Model evaluation (fit_theta_interext_atmparams)
- Grid interpolation (interpol_syn): directory tree
  `<grid>/ATLAS12/Z%.2f/HE%.3f/X%.2f/G%.3f/T%.0f.fits`, extension [1],
  column `c` (calibrated flux, erg/s/cm²/Å); wavelength grid
  `<grid>/ATLAS12/lambda.fits[1]` column `l` (5574 points, R≈1000, 800–180000 Å
  for sdB). Multilinear: per dimension the two enclosing nodes
  (wherefirst(x<arr); i==0→1; NULL→len−1), exact node match short-circuits;
  recursion order Z→HE→X→G→T; 2-point Lagrange weights
  `y0·(x−x1)/(x0−x1) + y1·(x−x0)/(x1−x0)`.
  Interpolated-spectrum cache keyed
  `ATLAS12Z%.4fHE%.5fX%.4fG%.5fT%.2fc[1]`, LRU capacity ncomp·6 (this
  quantizes nothing at mpfit's step sizes but replicate anyway).
- Multi-component: l_uni = union of grid lambdas restricted to
  [max(l_min_i), min(l_max_i)]; f_uni = Σ sur_ratio_i · linear-interp of f_i
  onto l_uni. Single component: f_uni = f on that grid's lambda. l_uni is
  FIXED during the fit (vrad/vsini/zeta all 0 in scope).
- Extinction (Fitzpatrick+2019): k(x=1e4/λ) = k0(x) + s(x)·(R55−3.02) with
  GSL **Akima** splines through the 102-point Table-3 data (tables in
  `src/f19_table.inc`); linear extrapolation for x≥8.7 from the last two
  points; factor = 10^(−0.4·E·(k+R55)). GSL was active in the reference runs
  (verified) — do not fall back to linear.
- Synthetic AB magnitude per filter: curves from
  `refdata/filter_passbands.fits.gz[<system>_<passband>]` (cfitsio reads .gz
  directly; label aliases DECaPS/DES_DR2/DELVE→DES, SMASH→SDSS,
  VegaHST/HST05/HST18→HST, VegaINT→INT). Curves are pre-normalized
  ∫SE d(−1/λ)=1. Box filters (system "box", passband "lo_hi") and TD1
  1565/1965/2365 (±165 Å) are analytic two-point boxes with
  f = 1/(1/lo − 1/hi). Integration: ind = points of l inside the filter,
  prepend one point left (the S-Lang append-right branch is dead code);
  u = union(l_filter, l[ind]); filter interpolated onto u; flux linearly
  interpolated (ISIS interpolate_d semantics: linear extrapolation from the
  two nearest subgrid points at the ends) onto u; trapezoid;
  mag = −2.5·log10(int) − 2.407948242680184, except Stroemgren
  Hbeta_narrow/wide which omit the constant. In C++ this is collapsed into
  fixed sparse weights on the flux array (src/synthmag.cpp).
- Then: mag += 1.505149978319906 − 5·logtheta (all rows); colors recomputed
  from the shifted mags (21 fixed rows: Johnson UmB/BmV/VmR/VmI/RmI,
  Stroemgren bmy/umb/c1/m1/Hbeta_B/Hbeta_AF (2.525+1.305·Hβ, 2.506+1.368·Hβ),
  Geneva UmB/B1mB/B2mB/V1mB/VmB/GmB/X/Y/Z with the Cramer&Maeder
  coefficients, AstroSat F148WmF169M); then mag += ZP.
- ZP table: `src/zp_table.inc`, 289 rows
  "system passband ZP ZP_err zpo_corr zpo_corr_err type"; with
  apply_ZPO_corr=1: ZP−=zpo_corr, ZP_err=sqrt(ZP_err²+zpo_corr_err²).
  Box rows get ZP=0, ZP_err=0.015. Color ZPs = linear combinations of
  magnitude ZPs (see src/passbands.cpp).

### Output formats
- `photometry_results.txt`: header lines exactly
  `# RA = %.8f; DEC = %.8f; norm_chi_red = %.4f; chisqr_reduced = %.3f;`,
  `# nmag_good = %d; grid = <last two path components of grid dir>;`,
  `# ruwe = %.3g; ipd_gof_harmonic_amplitude = %.3g; visibility_periods_used = %.3g; parallax = %.3g; parallax_error = %.3g; `,
  `# meanSFD = %.3g; meanSandF = %.3g; ` — then the print_struct table with
  columns index name value freeze min max conf_min conf_max buf_below
  buf_above tex; right-aligned, 3-space separator; doubles printed as S-Lang
  %S = shortest round-trip repr with ".0" appended to integral values,
  NaN prints as `-nan`; tex column via TeX_value_pm_error(...; sci=6)
  (ports in src/texval.cpp; `\ldots` for params without conf).
- `photometry_results_stellar_c{i}.txt`: name value conf_min conf_max
  median_value median_conf_min median_conf_max for c{i}_R/M/L
  (mode, HDI_lo, HDI_hi, median, quantile_lo, quantile_hi).

### Known open risks (things NOT yet verified — check during validation)
1. The MC-delta wildcard question flagged in step 13 above (star1/2 c1_teff
   deltas: 0 in the template vs prescribed-range in current main.cpp —
   compare against ref stellar spreads and make C++ match the reference).
2. mpfit convergence path: values should match to ~1e-6; if they differ more,
   check relstep/step wiring into mp_par and that frozen/derived params are
   excluded from the free set.
3. conf-limit port: verify star2's asymmetric E_44m55 interval and star3's
   c2_logg lower limit AT the grid bound (2.0) reproduce; bounds behave as
   prange endpoints (EVAL_INVALID → conf stays at the bound → tex shows red).
4. Star-3 c1_teff conf_min/max must equal the prescribed range exactly
   (29419.157 / 30300.247) via the propagation block.
5. `evaluate()` computes ALL prepared filters; prepared set = all 289 — fine
   numerically; if too slow later, prepare(mag_ind) after the guesses.
6. The `-nan` repr: S-Lang _NaN has the sign bit set; validation compares
   numerically so cosmetic differences are OK — don't chase string equality
   of the whole file, compare parsed numbers.

## File map

| file | contents |
|---|---|
| src/util.* | ISIS linear interp/extrapolation, trapz, union, S-Lang quantile, ISIS median (lower middle), string helpers |
| src/extinction.* | F19 Akima splines (GSL), extinction factor |
| src/fitsio_util.* | read double column via cfitsio (handles vector single-row cols) |
| src/modelgrid.* | grid.fits coverage (cols t,g,x,z,HHE), lambda, node cache, interpol_syn port, search_grid_dirs |
| src/passbands.* | ZP table + ZPO, filter curves (fits.gz + boxes + aliases), 21 colors + deps + color ZPs |
| src/synthmag.* | per-filter sparse integration weights on fixed l grid |
| src/fitfun.* | parameter table w/ ISIS defaults, set_par (glob, validate-no-clamp), model flux (multi-grid), evaluate → observed-system mags |
| src/fitter.* | photometric_fitting flow: guesses, mpfit, outliers, norm_chi_red Newton, conf port of fit-chisqrconf.c |
| src/photometry_table.* | photometry.dat reader, results-header (astrometry) reader |
| src/stellar.* | MC draws, mode_and_HDI + Savitzky-Golay, stellar constants |
| src/texval.* | round_err/round_conf (DIN 1333)/TeX_value_pm_error |
| src/results.* | S-Lang %S repr, print_struct-style tables, results/stellar writers |
| src/config.* | JSON schema (rejects unimplemented niche options: mass_can, derive_logg, hb_distance, derive_logg_c2, derive_sr, R1) |
| src/main.cpp | driver: astrometry → photometry → E guess → grids → fit → tex → propagation → dummies → writers → MC |
| src/zp_table.inc, src/f19_table.inc | auto-generated (scripts inline in session; regenerate by parsing photometric_magnitudes.sl / interstellar_extinction.sl) |
| extern/mpfit/ | vendored from ~/Projects/ISIS_install/src/isis/src/ |

## Reference sources (read these when in doubt)

- Template: `.../stellar_isisscripts/templates/photometry.sl`
- Fit: `.../src/photometry/photometric_fitting.sl` (photometric_fitting,
  initialize_grid_fit_photometry incl. the eval'd fit function + defaults)
- Magnitudes: `.../src/photometry/photometric_magnitudes.sl` (ZP table at top,
  integration at bottom), `filter_passbands.sl`, `photometric_colors.sl`
- Grid: `.../src/spectroscopy/interpol_syn.sl`
- Extinction: `.../src/photometry/interstellar_extinction.sl`
- Table I/O: `.../src/photometry/photometric_table.sl`,
  `.../src/isisscripts/print_struct.sl`
- Stats: `.../src/miscellaneous/mode_and_HDI.sl`, `savitzky_golay.sl`,
  `.../src/isisscripts/{round_err,round_conf,TeX_value_pm_error,quantile}.sl`
- ISIS engine: `~/Projects/ISIS_install/src/isis/src/mpfit-isis.c` (config),
  `fit-chisqrconf.c` (conf algorithm), `fit-params.c` (set_par validation),
  `share/fit-cmds.sl` (Δχ² levels 1.00/2.71/6.63, tol 1e-3),
  `src/util.c` interpolate_d (linear interp semantics),
  `src/math.c` find_median.
