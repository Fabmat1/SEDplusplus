# SED++ speedup plan (for tens of millions of stars)

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

