# SED++ speedup plan (for tens of millions of stars)

## Round 2 (2026-07-10): Phoenix-grid bulk fits, 2.6-9 s/star baseline

Profiled on 10 real bulk stars (Phoenix grid, 47k lambda points, conf_level=0,
nMC=20k). Time is dominated by the fit itself: thousands of model evaluations
(mostly conf-limit re-fits), each interpolating/extincting/integrating the
full 47k-point spectrum, with the extinction pow() loop alone ~50% of wall
time (invisible to gprof: libm is not instrumented).

Bit-exact changes (default path, all 10 bench stars byte-identical):

9.  **Subset fit path.** Only ~50-77% of the grid wavelengths lie inside any
    fitted filter (38% of Phoenix points are UV below every filter). The fit
    now interpolates the spectrum, computes the extinction factor and gathers
    the band integrals ONLY at those wavelengths (FitFunction::prepare builds
    the subset; ModelGrid::interpolate_sub mirrors the rounded-key LRU cache
    of the full path exactly, and the caches persist across the re-prepare
    after outlier removal -- both are required for bit-exactness). The full-
    spectrum path is untouched and still serves model-spectrum output;
    SEDFIT_NO_SUBSET=1 falls back to it.
10. **Band-integral memo.** Integrals are a pure function of (spectrum cache
    generation, E, R); a per-star memo makes logtheta-only evaluations
    (1 of 4 Jacobian columns) free.
11. **Single-pass grid interpolation.** interpol_fast_sub replaces the
    recursive Lagrange combine (string keys, temporaries, one array pass per
    dimension) with one fused pass over the <=8 corner nodes, replicating the
    per-element rounding order. Falls back to the recursion for >3
    interpolating dimensions.

Result: 10-star bench 92.8 s -> 49.4 s (1.9x), byte-identical outputs.

Opt-in, NOT byte-exact (documented accuracy trade; all default OFF):

- `fast_ext: 1` -- extinction factor via a vendored, deterministic,
  AVX2-vectorized 10^x (src/fastext.cpp) instead of libm pow: ~3e-16 relative
  difference per factor (vs mag errors of 1e-3..1e-1). On chaotic/degenerate
  conf loops even 1-ulp noise can flip a restart decision, so conf limits of
  degenerate stars may differ visibly; well-constrained fits are unaffected.
- `conf_tol: 1e-2` -- conf-limit search tolerance (default 1e-3 = ISIS
  parity); limits become ~1% coarse, noticeably fewer re-fits per limit.
- `max_conf_restarts: 1` -- CAUTION: on one bench star this truncated the
  conf improve-adopt loop early enough that c1_teff pegged at the grid edge
  (3152 -> 15000). Validate on a pixel before lowering below the campaign's
  current 10.

All three combined: 10-star bench 92.8 s -> 15.0 s (6.2x, 1.5 s/star; the
bench set is fit-heavy -- typical accepted MS-first stars are faster).
`SEDFIT_STATS=1` prints per-star evaluation/memo/ext counters.

## Round 1 (original)

Baseline (Gaia DR3 3431631727942970240, conf_level=0, nMC=2e6, no tex/plot):
**8.4 s total** = fit 1.0 s + Monte-Carlo statistics ~6.7 s + outputs ~0.7 s.

All items below are bit-exact (identical outputs) unless marked otherwise.
Validated after each step against reference outputs of the unmodified binary.

## Approaches, in implementation order

1. **`mode_and_HDI`: sort once.** The 2M-sample array is re-sorted in every
   one of the 10 refinement iterations, and `median`/`quantile` sort further
   copies. Sort once on entry, reuse everywhere. (~75 % of total wall time.)
   Status: DONE (benchmark star 8.4 s -> 2.0 s; outputs byte-identical)

2. **Prepare only the needed passbands.** `Fitter::select_data` calls
   `prepare(all 289 magnitude rows)`, so every `evaluate()` integrates every
   filter in the database and every filter curve is loaded, although only the
   observed bands (plus color dependencies) are ever read. Build the needed
   row set from the photometry table first and prepare only that.
   Cuts per-evaluate cost and filter-curve loading ~10x.
   Status: DONE (fit stage 1.0 s -> 0.1 s; outputs byte-identical)

3. **`PassbandDB::compute_colors`: precompute row indices** in the
   constructor instead of 20 linear string `find()`s per model evaluation.
   Status: DONE (outputs byte-identical)

4. **Cache the extinction factor keyed on (E_44m55, R_55).** mpfit varies one
   parameter at a time, so most evaluations reuse the reddening vector
   (5574 `pow()`s each). Recompute only when E or R changed; same arithmetic
   when recomputed.
   Status: DONE (outputs byte-identical)

5. **Avoid the full-spectrum copy in `model_flux`** (single-grid case copies
   44 KB per evaluation). Return a reference to the cached spectrum.
   Status: DONE (outputs byte-identical)

6. **Load filter curves in one archive pass.** `PassbandDB::filter()` opens
   `filter_passbands.fits.gz[EXT]` per filter; cfitsio decompresses the whole
   archive each time. Open once, cache all extensions (single inflate).
   Status: DONE (14-star set 28.6 s -> 15.6 s; outputs byte-identical)

7. **`ModelGrid` interpolation cache: `list::splice` on hit** (no node
   allocation / string copy per evaluation), larger default capacity.
   Status: DONE splice only. Raising the capacity is NOT bit-exact: the
   cache key is fixed-precision, so capacity is observable behaviour (a
   larger cache reuses spectra for parameters differing below the key
   precision). Kept the ISIS capacity of 6 per component.

8. **Batch mode (`--multi <configlist.txt>`): amortize per-process setup.**
   One config path per line ('#' comments allowed); grids, node FITS files,
   and the filter archive are loaded once per process; errors in one star are
   reported and the batch continues (exit code 1 if any star failed).
   The rounded-key spectrum interpolation cache is cleared between stars, so
   each star's outputs are byte-identical to a single-config run.
   For parallelism, split the list into chunks and run one process per core
   (bit-exact, no shared mutable state between processes).
   Status: DONE (outputs byte-identical to single-config runs)

## Measured results (14-star validation set, conf_level=0)

| scenario | before | after |
|---|---|---|
| original validation configs (nMC=2e6, write_model=1) | 120.8 s | 15.0 s |
| production-style (nMC=2e5, write_fits only), single | 23.0 s | 1.6 s |
| production-style, --multi | - | 1.5 s (~105 ms/star) |

Benchmark star (Gaia DR3 3431631727942970240, nMC=2e6): 8.4 s -> 1.1 s,
fit stage 1.0 s -> <0.05 s. All outputs (results/stellar/fit-mag/fit-col txt
and every SED_results.fits table) byte-identical to the unmodified binary
across all 14 stars, in both single and --multi mode.

At ~105 ms/star, 10^7 stars are ~12 CPU-days -> ~9 h on 32 cores (split the
config list into per-core chunks and run one `sedfit --multi` each). The
remaining per-star cost is dominated by the stellar-parameter Monte Carlo
(scales linearly with nMC) and FITS output writing.

## Not done on purpose (would change results)

- Relaxing mpfit ftol/xtol/gtol (1e-10) or maxiter.
- Lowering nMC / parallelizing the MC draws (would change the RNG sequence).
  If statistical (not bit-exact) reproducibility is acceptable, drawing with
  per-thread RNGs would cut the remaining ~1 s MC cost by the core count;
  so would simply setting nMC=2e5 in the config (MC error of the
  mode/HDI scales ~1/sqrt(n); 2e5 is usually plenty).
- `-ffast-math` / reordering FP arithmetic.


## Round 3: bulk fast paths (NOT ISIS-parity; all opt-in, bulk_fit defaults)

- `error_mode: "covar"`: parameter errors from the mpfit covariance
  (xerror * sqrt(delta_chi2)), skipping the conf-limit + renormalize loops.
- `nMC: 0`: analytic asymmetric R/M/L error propagation instead of the MC.
- `use_mag_grid: 1`: precomputed per-band mag grids (`sedfit --premag
  <griddir> <refdata> [box...]` -> <griddir>/mag/ with band integrals +
  extinction moments; manifest hash of the filter archive gates staleness,
  flux-path fallback when a needed band is missing). Model mags agree with
  the flux path to ~1 mmag at E<=0.3, few mmag at E~1.
- bulk_fit.py: jobs "auto" = all cores, premag preprocessing phase,
  error-aware MS mass-radius acceptance (mr_max_rel_err gate).

Measured (Phoenix MS star, single fit): 3.0 s (conf) -> 0.4 s (covar) ->
0.03 s (covar + mag grid). 40-star pixel-0 campaign, 16 jobs: fits
57.1 s -> 2.1 s (27x). Parity path untouched: validation/compare.py PASS.
Caveats: without the conf loop's restart-on-improvement a few degenerate
stars settle in different (occasionally worse, occasionally better) minima;
covar intervals on degenerate stars are wider/more honest than capped-conf.
