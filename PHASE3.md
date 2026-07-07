# SED++ Phase 3 — Optional output products (TeX, SED_results.fits, plot)

**Handoff/spec document for Phase 3.** Phase 1 (fit engine) and Phase 2 (query
layer incl. astrometry/SIMBAD/spectra) are complete and validated — see
`PLAN.md`, `PHASE2.md`, `PHASE2_SPECTRA.md`. Phase 3 recreates the optional
output products the ISIS template writes after a fit:

| product | ISIS reference | fixture (validation target) |
|---|---|---|
| `photometry_fit_mag.txt`, `photometry_fit_col.txt`, `photometry_fit.txt` | `templates/photometry.sl:1306-1476` | `validation/ref_*/photometry_fit_*.txt` |
| `photometry_results.tex` (+ pdflatex → `.pdf`) | `photometry.sl:845-1304` | `validation/ref_*/photometry_results.tex` |
| `SED_results.fits` | `photometry.sl:1585-1664` | `validation/ref_*/SED_results.fits` |
| `photometry_SED.pdf` (Python matplotlib, replaces xfig) | `photometry.sl:1380-1525`, `src/plotting/xfig_photometry.sl` | `validation/ref_*/photometry_SED.pdf` (visual) |

All ISIS reference paths below are under
`~/Projects/ISIS_install/src/stellar_isisscripts/`.

## Config toggles (ALL OFF by default — bulk fitting must stay lean)

Add to `src/config.{hpp,cpp}` (Config already has `write_model`/`save_MC`
fields; change defaults and add the new ones):

- `write_model` (default **false**; ISIS default 1): build the model spectrum
  and write `photometry_fit.txt`, `photometry_fit_mag.txt`,
  `photometry_fit_col.txt`; also adds the `filters`/`colours`/`spectrum_fit`
  extensions to `SED_results.fits` when `write_fits` is on.
- `write_fits` (default **false**; ISIS: always on): write `SED_results.fits`
  (primary header + `params_fit` + `stellar_c*`; model extensions iff
  `write_model`).
- `save_MC` (default false, same as ISIS): add `MC_c1`/`MC_c2` extensions
  (R/L/M raw MC arrays) to `SED_results.fits`. Only meaningful with
  `write_fits`.
- `write_tex` (default **false**; ISIS: always on): write
  `photometry_results.tex` and run `pdflatex` (skip with a warning if pdflatex
  is not in PATH).
- `plot` (default **false**): after writing `SED_results.fits`, run the Python
  plot script to produce `photometry_SED.pdf`. Implies `write_fits` +
  `write_model` internally.

Validation configs `configs/<gaiaid>.json` get all toggles enabled.

## Stage 1 — model spectrum + write_model text products

**New files:** `src/spectrum.{hpp,cpp}` (a DRAFT already exists in the tree —
copied from a parked worktree; it ports `convolve_syn`'s C `convolve` kernel,
the S-Lang `[a:b:dx]` range semantics, and the multi-component union. VERIFY
it against the S-Lang before trusting it), plus a new module for
`magnitudes_to_flux` (suggest `src/magflux.{hpp,cpp}`).

Template flow (photometry.sl:1306-1476), scope = branches actually taken
(no user SED, no blackbody components — throw on bb_teff!=0&&bb_sur_ratio!=0):

1. `res_offset = 0, res_slope = 1./6.0; xmin = 1100, xmax = 59000`.
2. `(l,f) = syn_spec(c*teff, c*logg, c*xi, c*_z, c*HE, grid_info;
   sur_ratio=c*sur_ratio, metals=ATLAS12, res_offset, res_slope,
   xmin=xmin-10, xmax=xmax+10)` — reference `src/spectroscopy/syn_spec.sl`
   (262 lines) + `src/spectroscopy/convolve_syn.sl` (416 lines). Each
   component is interpolated on its own grid (reuse
   `ModelGrid::interpolate`), restricted to [xmin-10, xmax+10], convolved to
   constant FWHM = 1/res_slope = 6 Å (Gaussian), then components are combined
   on the union grid weighted by sur_ratio (mirrors the fit's l_uni logic).
3. Uncertainty inflation for the figure/tables (photometry.sl:1344): if
   `norm_chi_red != 0`, `uncertainty = sqrt(uncertainty² + norm_chi_red²)` —
   apply to the photometry copy used for mout/plot ONLY (results.txt and the
   stellar MC are already written by then; keep it that way).
4. `mout = magnitudes_to_flux(l, f, entries; theta=10^logtheta, E_44m55,
   R_55, errbar, flagged)` — reference `src/plotting/xfig_photometry.sl:3-600`.
   Output tables (fixture schemas):
   - `mout.mag`: lambda_min, lambda, lambda_max, flux_min, flux, flux_max,
     diff, diff_err, passband, system, flag, VizieR_catalog (29 rows for
     star 1 = ALL entries incl. flagged, `flagged` qualifier).
   - `mout.col`: diff, diff_err, passband, system, flag, VizieR_catalog.
   Written to `photometry_fit_mag.txt` / `photometry_fit_col.txt` (col file
   only if `length(mout.col.diff)>0`) via the existing `write_table` +
   `slang_repr` (print_struct format, exactly like the existing writers).
5. `sf = interstellar_extinction(l, E_44m55, R_55) * (0.5*10^logtheta)^2`;
   `sout = {l, f*sf}`; for multi-component fits per-component curves
   `f_c{i} = interpol(l, lt_i, ft_i) * sf` where `(lt_i, ft_i)` = syn_spec
   with sur_ratio vector zeroed except component i (photometry.sl:1429-1452).
   Written to `photometry_fit.txt` (columns l, f, then f_c1, f_c2 when
   present; struct_combine order).
6. Wire into `main.cpp` after the stellar-MC block, gated on
   `cfg.write_model || cfg.plot` (keep mout/sout in scope for stages 3/4).

**Validation:** run `./build/sedfit configs/<id>.json` for the three ref
stars (outdir != workdir!) and compare against
`validation/ref_<id>/photometry_fit*.txt`: the `l`/`lambda*` columns and all
string/flag columns should match byte-for-byte; `f`/`flux*`/`diff*` to ~1e-6
relative (fit parameters differ at the 1e-7 level from ISIS — PLAN.md).
Star 3 (872072016870129024) exercises f_c1/f_c2 and the colours file.

## Stage 2 — TeX output (photometry.sl:845-1304)

**Prereq refactor** in `src/stellar.{hpp,cpp}` + `main.cpp`: the tex table
needs, per component, the same MC that already produces R/M/L, PLUS:
- `vgrav = M_Msun*Const_GMsun/(R_Rsun*Const_Rsun*c) / 1e5` km/s and
  `vesc = sqrt(2*10^logg_draw * R_Rsun*Const_Rsun) / 1e5` km/s
  (c = 2.99792458e10 cm/s), computed AFTER the parallax/sur_ratio filter but
  BEFORE the per-quantity positivity filters, then each filtered >0
  independently and run through mode_and_HDI — exactly parallel to R/M/L
  (photometry.sl:1176-1222). Extend StellarMCResult with vgrav/vesc.
- optionally-retained filtered R_Rsun/L_Lsun/M_Msun arrays (for `save_MC`,
  stage 3). Gate retention on a flag in StellarMCInput to avoid 2e6-double
  copies in bulk mode.
- Gaia distance MC (photometry.sl:936-945): draws
  `parallax + grand()*parallax_error`, `mode_and_HDI(1e3/draws; p=confidence)`
  → mode/HDI and median/quantile rows. (Separate small helper; only needs the
  parallax pair + confidence + nMC.)

**Writer** (new `src/texout.{hpp,cpp}`), exact port of the fixture layout —
read `validation/ref_648029239761019776/photometry_results.tex` (single
component; "(fixed)" and free rows, no propagation) and
`validation/ref_872072016870129024/photometry_results.tex` (two components,
"(prescribed)" rows via par_full literal names + sci=6, red-at-bound values
like `$2.0^{+1.0}_{-{\color{red}0.0}}$`, per-component multicolumn headers)
against the template code before writing any C++. Key semantics:

- Table head: `Object: <star> & 68\% confidence interval\\` (90/99 for
  conf_level 1/2; conf_level 1/2 also scales astrometry.parallax_error by
  1.645/2.576 ONCE, affecting the parallax row, distance MC, and stellar MC —
  the C++ stellar path already applies perr_scale; keep single-scaling).
- Reddening rows: SFD and S&F via
  `TeX_value_pm_error(mean, mean-std, mean+std)` (default sci), with the
  exact `\href` strings from the fixture. Stilism/Bayestar rows only when
  those fields exist (they don't in our fixtures). No rstilism row (NULL).
- E_44m55 / R_55 / logtheta rows (photometry.sl:898-922): tex column if not
  `\ldots`; else free → `$%g$` (wrapped `\color{red}{...}` when value sits at
  min/max); else `(fixed) & $%g$`.
- Parallax row: `Parallax $\varpi$ ({\it Gaia}, $\text{RUWE}=%.2f$) & %s\,mas`
  with `TeX_value_pm_error(round_conf(parallax, parallax_error))` (round_conf
  is already in texval.cpp). The parallax_offset branch is NOT taken (our
  Astrometry has no such field).
- Distance mode/median rows (sci=3), only when good_astrometry && parallax>0.
- Per-component blocks (photometry.sl:973-1238): field order teff, logg, xi,
  z, HE, sur_ratio (skip sur_ratio for component 1); "(prescribed)" =
  par_full has the LITERAL name `c{i}_{field}` with freeze==1 →
  `TeX_value_pm_error(value, min, max; sci=6)`; then R/M/L mode+median rows
  (sci=3, exact strings incl. `\phantom{...}` from the fixtures — note the
  c2 radius label carries `(A_\mathrm{eff}/A_{\mathrm{eff,}1})^{1/2}`),
  vgrav/vesc mode rows (sci=3).
- Blackbody blocks: gate is `bb_teff!=0 && bb_sur_ratio!=0` — always false
  here; implement the gate only.
- Footer: `Generic excess noise ... $%.3f$\,mag`, `Reduced $\chi^2$ ... $%.2f$`,
  `\hline`, `\end{tabular}`, `\end{document}`.
- pdflatex invocation mirroring photometry.sl:1303 (system(), grep error,
  rm aux/log/out), gated on `write_tex`, warn+skip if pdflatex absent.

**Validation:** diff C++ tex vs the three fixtures. Non-MC lines must match
byte-for-byte. MC-derived lines (Distance/R/M/L/vgrav/vesc) may differ in the
last printed digit (~0.1-1% MC noise floor, PLAN.md); anything larger is a
bug. vgrav/vesc have no txt-file precedent — check them against the fixture
tex values (they come from the same validated MC arrays, so agreement should
be at the same ~1% level).

## Stage 3 — SED_results.fits (photometry.sl:1585-1664)

**Writer helpers** in `src/fitsio_util.{hpp,cpp}` (cfitsio): create-with-
clobber, empty primary image HDU (fixture: BITPIX=32, NAXIS=1, NAXIS1=0) with
keywords, and a generic binary-table writer (double / int / string columns;
string TFORM width = max strlen over the column, e.g. `13A`/`21A` in the
fixture's params_fit).

Primary keywords, in fixture order: RA, DEC, NORM_CHI_RED, CHISQR_RED,
NMAG_GOOD, GRID_SHORT, then Gaia fields (ruwe, ipd_gof_harmonic_amplitude,
visibility_periods_used, parallax, parallax_error — only those that are
non-NaN, photometry.sl:1595-1600), then MEANSFD, MEANSANDF, MEANSTILISM.
Long names go through HIERARCH automatically. Notes:
- Reddening fields are created for all of hfields_redd and left UNDEFINED
  when NULL (fixture shows `MEANSTILISM =` with no value) — write an
  undefined/null keyword in that case.
- `NORM_CHI_RED = 0` in the fixture is an integer (S-Lang norm_chi_red stays
  Integer_Type 0 when no renormalization ran) — write as int when ==0.
- Match double formatting against the fixture via astropy repr (e.g.
  `CHISQR_RED = 0.4715880324209`); check how the ISIS cfitsio module writes
  double keys (isis `fits_write_key`) if digits don't line up.

Extensions, in order:
1. `params_fit`: index (J, 1-based), name (A), value (D), freeze (J), min,
   max, conf_min, conf_max, buf_below, buf_above (D), tex (A) — straight from
   FitResults.
2. `stellar_c{i}` (one per grid, only when good_astrometry): the StellarRow
   columns (name/value/conf_min/conf_max/median_value/median_conf_min/
   median_conf_max) + `GRID` keyword = grid_short(griddirectories[i]).
3. iff write_model: `filters` (mout.mag), `colours` (mout.col, only if
   nonempty), `spectrum_fit` (sout: l, f [, f_c1, f_c2]).
4. `IUE` (wavelength/flux/error=sigma + FILE keyword): only when an IUE
   spectrum was loaded for the plot (photometry.sl:1487-1518 selects the
   file from box-filter entries: VizieR_catalog `*VI/110/inescat` directly,
   or `Average` → the flag==2 row whose magnitude is closest; file read from
   the local `IUE/` cache that src/spectra.cpp already populates, rows
   filtered to flux>0 && quality==0). NO fixture exercises this (the ref
   stars have no box rows) — implement faithfully, mark unvalidated.
5. iff save_MC: `MC_c1`/`MC_c2` with columns R_Rsun, L_Lsun, M_Msun (the
   filtered arrays retained by stage 2's refactor). The three arrays can
   have different lengths after per-quantity filtering; check how ISIS
   `fits_write_binary_table` handles ragged struct fields before choosing
   (likely pads); no fixture has these extensions (save_MC=0) — document
   whatever is chosen.

**Validation:** python3/astropy structural + numeric compare against the
three fixture files: HDU names/order, column names/TFORMs, header keywords;
values exact for strings/ints, ~1e-6 rel for fit-dependent doubles, ~1% for
MC-derived stellar values. Suggest a `validation/compare_fits.py`.

## Stage 4 — Python SED plot (replaces xfig_photometry rendering)

`scripts/plot_sed.py` (python3 + matplotlib + astropy, all installed):
reads `SED_results.fits` alone and writes `photometry_SED.pdf`. Replicate the
xfig_photometry figure (reference `src/plotting/xfig_photometry.sl:601-1149`
for layout semantics; `validation/ref_*/photometry_SED.pdf` for ground
truth):
- Main panel: x = λ [Å], log scale, range [xmin=1100, xmax=59000-1];
  y = f·λ^3 / 10^factor where factor is chosen power-of-ten (xfig `factor`
  autoscale); model spectrum (spectrum_fit l, f) as black line; per-component
  curves f_c1 (skyblue), f_c2 (pink) when present; observed points from
  `filters`: y from flux·λ³, x = lambda, horizontal bars lambda_min→lambda_max
  (filter width), vertical error bars flux_min→flux_max, colored per
  photometric system (xfig `colored` qualifier assigns per-system colors —
  eyeball the fixture PDFs), flagged (flag!=0) rows omitted the way the
  template does (it passes `flagged` to magnitudes_to_flux but the plot only
  shows good rows — check xfig_photometry's handling of the flag column).
- χ panel(s) below (`chi` qualifier): (obs−model)/err = diff/diff_err per
  magnitude row vs λ, same log x; a second χ panel for colours when the
  `colours` extension exists (fig length 2 vs 3, photometry.sl:1521-1522).
- IUE overplot in violet when the `IUE` extension exists.
- y-axis label `λ³·f (10^factor Å³ erg s⁻¹ cm⁻² Å⁻¹)`-style, matching the
  fixture's labeling.
- Wiring: `plot` toggle in the config makes sedfit invoke the script
  (`python3 <script> <SED_results.fits path> <output pdf path>`); locate the
  script relative to the sedfit binary with a `plot_script` config override.
**Validation:** render fixture-vs-C++ PDFs to PNG (pdftoppm) and compare
visually (an agent can Read images) for the three stars: same curves, point
positions, panel structure. Exact typography parity is NOT required.

## Working notes (agents: append status here)

- 2026-07-07: PHASE3.md written; draft `src/spectrum.{hpp,cpp}` copied from
  parked worktree `.claude/worktrees/agent-a281e0290b6344cf7` (unverified
  against convolve_syn.sl — verify in stage 1).
- Build: `cmake -B build && cmake --build build -j` (add new .cpp files to
  CMakeLists.txt). Grids live under `/home/fabian/ISIS_models/`; refdata path
  comes from the configs. Run fits from the repo root.
- The three validation configs are `configs/648029239761019776.json` etc. if
  present from Phase 1 — check `configs/`; create them per PLAN.md §2 if the
  outdir/cpp_* variants are missing, and ENABLE the new toggles there (never
  point outdir at the ref_* fixture dirs).

- 2026-07-07: **Stage 1 DONE + validated.** Files: verified draft
  `src/spectrum.{hpp,cpp}` against `convolve_syn.sl`/`syn_spec.sl` +
  `slirp/c_functions.c` `convolve_single` (faithful; the branch taken is
  Longslit res_offset=0 with macro-broadening OFF → clambda equidistant step
  0.5/(2*res_slope)=1.5 Å, glambda 6*max(nint(sig/0.015),15)+1 pts, then
  `intpol(...; akima)` back onto union(lambda,clambda) — draft matches). New
  `src/magflux.{hpp,cpp}` ports `magnitudes_to_flux` incl. the ISIS
  `complement`/`intersection` (setfuns.sl) row-ordering quirk (all-measurable →
  original order; else value-sorted) and the `m.VizieR_catalog[i]` loop-index
  indexing bug. Writers `write_mag_txt`/`write_col_txt`/`write_spectrum_txt`
  added to `src/results.{hpp,cpp}` (reuse `write_table`+`slang_repr`). Toggles
  added to `src/config.{hpp,cpp}` (write_model/write_fits/write_tex/plot all
  default false; plot⇒write_fits+write_model). Wired in `main.cpp` after the
  MC block, gated on `write_model||plot`. Added the 2 .cpp to CMakeLists.
  Comparator: `validation/compare_fit_txt.py`.
  * **GOTCHA (important for stages 3/4):** `src/synthmag.cpp` (the FIT engine)
    builds its integration grid `u` from the PREPEND-extended `ind`, whereas
    `photometric_magnitudes.sl:767-768` builds `u` from the UN-extended `ind0`
    (the left-prepend point is only for the flux interpolation, not added to
    `u`). This makes SynthMag extrapolate the filter one node past its blue
    edge → ~3e-5 mag error. It is invisible to the Phase-1 fit validation
    (data errors ≫ 3e-5) but WRONG for the plot magnitudes. magflux therefore
    computes its synthetic AB mags with its OWN faithful local port
    (`synth_ab_mag`) rather than reusing SynthMag — the fit engine was left
    untouched per the task constraint. If a later stage needs model AB mags,
    use magflux's port, NOT SynthMag (or fix SynthMag's `u` and re-validate
    photometry_results.txt).
  * **Validation (`./build/sedfit configs/<id>.json`, cpp_* vs ref_*):**
    all three PASS. Per column, worst deviations:
    - 648 (single): flux/flux_min/flux_max rel ≤2.4e-15; lambda_eff rel
      1.2e-8; diff abs 3.3e-7 mag; diff_err/lambda_min/lambda_max/strings/flag
      EXACT; fit.txt l EXACT, f rel med 1.8e-7 max 5.4e-7. col file (2 rows)
      diff abs 6.1e-8.
    - 3431 (single, norm_chi_red=0.0087): flux rel ≤7.1e-7; diff abs 3.1e-7;
      diff_err abs 7.7e-7 (from norm_chi_red matching only to its 4-digit
      print); lambda_eff rel 2.3e-8; fit.txt f rel med 2.8e-7 max 1.1e-6.
    - 872 (composite, exercises f_c1/f_c2): flux data points rel ≤2.9e-13
      (obs-driven, so insensitive to the composite floor!); lambda_min/max,
      diff_err, strings, flag EXACT; but diff abs up to 1.3e-3 mag, lambda_eff
      rel 4e-5, and fit.txt f/f_c1 max rel ~2-3% (median 4e-5/1.3e-3). These
      larger model-side deltas track the ACCEPTED Phase-1 composite floor:
      E_44m55 is degenerate at the E=0 boundary (ref 8.1e-5 vs cpp 4.7e-5),
      logtheta/c2_teff/c2_sur_ratio differ at 1e-5..1e-4. Phase-1
      `compare.py` still PASSes all three (fit engine unchanged).
  * The obs-flux data points are model-scaling-invariant: flux_obs =
    10^(-0.4*(m_obs-n_mag))·flux_weighted, and the model f cancels between
    flux_weighted (∝∫f·SE) and 10^(0.4·n_mag) (∝1/∫f·SE) — hence flux columns
    hit machine precision even when the model spectrum differs (872). The
    lambda_min/lambda_max come straight from the filter curve (FWTM) → exact.
  * Deviations from spec: none functional. The single quirk is the SynthMag
    avoidance above.

- 2026-07-07: **Stage 2 DONE + validated.** Files: `src/stellar.{hpp,cpp}`
  (vgrav/vesc ModeHDI in StellarMCResult computed from the EXISTING draws —
  no new RNG calls, sequence unchanged; `retain_arrays` flag keeps the
  filtered R/L/M draws for stage 3's MC_c* extensions; new
  `gaia_distance_mc()` with its own RNG stream, seed = cfg.mc_seed, comps use
  mc_seed+1..n); new `src/texout.{hpp,cpp}` (exact port of
  photometry.sl:845-1304 for our scope + `run_pdflatex()` that cd's into
  outdir, pipes through `grep 'error'`, rm -f aux/log/out, warns+skips if
  pdflatex absent); `src/main.cpp` (MC loop now collects
  `std::vector<StellarMCResult> mc_all`, stellar-txt writing unchanged; tex
  block gated on cfg.write_tex after the MC loop, before write_model —
  matches ISIS ordering); `write_tex: 1` added to the three ref configs;
  texout.cpp added to CMakeLists.
  * **BUG FIX in `src/texval.cpp` round_conf:** `std::min` → `std::fmin`.
    A negative error (value outside [lo,hi], tolerated by round_conf.sl:29
    "live_dangerously") makes round_err return digit=NaN; S-Lang
    `_min(NaN,x)==x` but `std::min(NaN,x)==NaN`, which poisoned the whole
    string. Only reachable for negative errors (never hit by previously
    validated outputs; fmin==min otherwise). Exercised by 872's prescribed
    c1_HE row `$-5.05^{+\phantom{0}0.06}_{--0.05}$` — now byte-identical.
  * **Validation (diff cpp vs ref photometry_results.tex):**
    - 648: all non-MC lines byte-identical. MC diffs (last digit only):
      Distance mode `$640^{+23}_{-22}$`→`$639^{+23}_{-22}$`; Distance median
      `+23`→`+24`; Radius mode `$0.196^{+0.008}_{-0.007}$`→`$0.197 \pm
      0.007$` (rounding-boundary collapse to symmetric form). PDF builds
      cleanly (42 kB, PDF 1.7); aux/log/out removed.
    - 3431: all non-MC lines byte-identical; ONE MC diff: Distance median
      `+28`→`+29`. Note the vesc line's `\left(1.045 \pm
      0.017\right)\times10^{3}` scientific form reproduced exactly.
    - 872: prescribed rows (c1_teff/logg/HE incl. the negative-error HE),
      red-at-bound c2_logg, multicolumn headers, phantom paddings all
      byte-identical. MC last-digit diffs: c1 R mode 0.142→0.141 (-0.010→
      -0.009), c1 M mode -0.029→-0.030, c1 L mode 14.3→14.2 (+2.6→+2.7),
      c2 M mode +0.0048→+0.0049 (-0.0011→-0.0010), c2 L mode +0.14→+0.13,
      c2 vgrav 0.0020→0.0019 (+0.0036→+0.0037). Three fit-level diffs, NOT
      tex bugs: Surface ratio 35.6→35.7, E(44-55)
      `$0.00008...$`→`$0.00005...$`, logtheta `+0.0330/-0.0026`→
      `+0.0331/-0.0025` — all byte-identical to the tex column of cpp's own
      photometry_results.txt, i.e. the ACCEPTED Phase-1 composite floor
      (E degenerate at the E=0 bound).
  * Regression: `validation/compare.py` all three PASS;
    `validation/compare_fit_txt.py validation/ref_<id> validation/cpp_<id>`
    (note: takes DIRS, not ids) all three PASS. Stellar txt values unchanged
    (RNG sequence untouched).
  * Gotchas for stages 3/4: retained MC arrays are per-quantity filtered so
    R_arr/L_arr/M_arr lengths can differ (see stage-3 ragged-column note);
    they are only populated when `save_MC` — configs must set it for the
    MC_c* extensions. TexInput holds a POINTER to the mc vector (avoids
    copying 3x2e6 doubles when save_MC is on; keep mc_all alive through the
    stage-3 FITS write). `gaia_distance_mc` does NOT filter parallax>0 draws
    (faithful to photometry.sl:939-940).

- 2026-07-07 (Stage 4 plot, done EARLY by the manager session, out of order):
  `scripts/plot_sed.py` (python3 + astropy + matplotlib, no other deps)
  renders photometry_SED.pdf from SED_results.fits alone:
  `python3 scripts/plot_sed.py <SED_results.fits> <out.pdf>`.
  * Ports xfig_photometry semantics: exponent=3, xmin/xmax=1100/58999,
    factor=floor(log10(ymax)) power-of-ten factoring, +-0.1-range y padding,
    chi panes clamped to +-(2.25+4%), dashed 0 / dotted +-1 guides, flag==0
    rows only, per-system colors + label offsets + passband prettification
    from filter_color_xypos (xfig_photometry.sl:187-600), f_c1 skyblue /
    f_c2 pink3 component curves, IUE ext overplot in violet (unexercised),
    side chi_color pane only when a colours ext exists.
  * Validated VISUALLY against all three ref fixtures by rendering the
    fixture SED_results.fits with the script and comparing PNGs side by
    side with the reference photometry_SED.pdf: layout, colors, point/label
    positions, filter-width bars, component curves, chi patterns all match;
    typography (txfonts vs stix mathtext) and tick density differ, accepted
    per spec.
  * STILL OPEN (belongs to Stage 3 integration): sedfit does not yet invoke
    the script — wire the `plot` toggle after the SED_results.fits writer
    exists (plot implies write_fits+write_model; invoke
    `python3 <repo>/scripts/plot_sed.py <outdir>/<basename>SED_results.fits
    <outdir>/<basename>photometry_SED.pdf`, locate the script relative to
    the sedfit binary with a `plot_script` config override, warn+skip if
    python3/astropy is missing).

- 2026-07-07 (Stage 3, done by the manager session): SED_results.fits writer
  + `plot` toggle wiring. PHASE 3 COMPLETE.
  * Files: `src/fitsio_util.{hpp,cpp}` (+FitsWriter: clobbering create,
    BITPIX=32/NAXIS1=0 primary, %.15G double keys via cfitsio defaults,
    undefined-key support, generic D/J/A binary tables with max-strlen
    string widths); `src/fitsout.{hpp,cpp}` (SED_results.fits per
    photometry.sl:1585-1664 + IUE selection port of :1487-1518);
    `src/main.cpp` (stellar_all collection, hoisted mout/spec products,
    write_fits block after tex, plot invocation via
    /proc/self/exe/../../scripts/plot_sed.py with `plot_script` config
    override); `src/config.{hpp,cpp}` (plot_script); CMakeLists.txt;
    `validation/compare_fits.py` (new comparator); configs got
    `"write_fits": 1, "plot": 1`.
  * Faithfulness details: NORM_CHI_RED written as int when 0 (S-Lang
    Integer_Type); VISIBILITY_PERIODS_USED written as int when integral
    (ISIS astrometry struct carries it as Integer — fixture card `= 19`);
    MEANSFD/MEANSANDF/MEANSTILISM cards always present when any reddening
    exists, unset ones as undefined keywords; GRID keyword per stellar_c*;
    FILE keyword on the IUE ext; MC_c* columns R_Rsun/M_Msun/L_Lsun
    (struct-field order), ragged arrays padded with NaN (ISIS behaviour
    unverified — save_MC=0 in every fixture).
  * Validation: `compare_fits.py ref cpp` PASS for 648 + 3431 (structure,
    TFORMs incl. string widths, keyword presence/format, values); 872 PASS
    with FIT_RTOL=0.04 (accepted Phase-1 composite E=0 floor: CHISQR_RED
    1.3e-5 rel, faint-component spectrum ~2-3%, tex rows 0/1/20 — all
    byte-identical to cpp's own results.txt). Card-image notes are the same
    fit-floor value digits, not formatting differences. Regression:
    compare.py PASS, compare_fit_txt.py PASS x3. save_MC smoke test (100k
    draws, scratch outdir): MC_c1 ext with correct columns and values.
    plot toggle end-to-end: photometry_SED.pdf produced for all three
    stars from the C++ FITS; composite plot visually matches the fixture.
  * UNVALIDATED: IUE extension (no fixture star has box rows); the template
    pairs the FILE keyword with whatever IUE_file was last SET (even if
    unreadable) — C++ pairs it with the last successfully READ spectrum.
