# SED++ (`sedfit`)

A C++ reimplementation of the `stellar_isisscripts` photometric SED-fitting
pipeline (`templates/photometry.sl`), used to derive angular diameters, colour
excesses, atmospheric parameters, and stellar radii/masses/luminosities of hot
subdwarfs and similar objects from broad-band photometry.

This is **phase 1**: the fit engine, reaching output parity with the original
S-Lang/ISIS pipeline. Photometry and astrometry are read from cached files
produced by the original pipeline (`photometry.dat`, and the header of
`photometry_results.txt`). Remote catalogue querying, a local SQL photometry
source, bulk/parallel multi-star fitting, plots, and the TeX/PDF and
`SED_results.fits` outputs are **not** part of this phase (see
[PLAN.md](PLAN.md), step 5).

## Status

Validated to parity against the original pipeline for three Gaia DR3 stars
(two single-component sdB fits and one composite sdB + late-type fit) — see
`validation/` and run `python3 validation/compare.py`. Agreement:

- Fit parameter values: single-component fits to ~1e-7 relative; the composite
  fit to ≤1.3% of each parameter's 1σ confidence interval (its reddening is
  degenerate, pinned against the E=0 boundary).
- Confidence limits: to ~1e-5 of the interval width (single); ≤6% (composite).
- Reduced χ² and generic excess noise: to the printed precision.
- MC-derived radius/mass/luminosity: to <1% (Monte-Carlo noise floor).

## Build

Requires a C++20 compiler, CMake ≥ 3.20, and (via pkg-config) `cfitsio`, `gsl`,
and `nlohmann_json`. `mpfit` (the same MINPACK-derived Levenberg-Marquardt
engine ISIS uses) is vendored in `extern/mpfit/`.

```sh
cmake -B build
cmake --build build -j
```

## Run

```sh
./build/sedfit config.json
```

The JSON config mirrors the header block of `photometry.sl`. See
`configs/*.json` for the three validated examples. Key fields:

| field | meaning |
|---|---|
| `star` | object identifier (informational in phase 1) |
| `coordinates` | `{ra, dec}` in degrees (informational in phase 1) |
| `par` | `{name, value, freeze}` — parameters to set (names may use `*` wildcards) |
| `par_full` | `{name, value, freeze, min, max}` — parameters to set with ranges |
| `griddirectories`, `bpaths` | model-grid subdirectories and base paths to search |
| `conf_level` | −1 (none), 0 (68%), 1 (90%), 2 (99%) |
| `apply_ZPO_corr` | apply empirical zero-point-offset corrections |
| `remove_outliers` | σ threshold for outlier rejection (0 disables) |
| `nMC` | Monte-Carlo trials for R/M/L |
| `refdata` | path to the `stellar_isisscripts/refdata` directory (filter curves) |
| `workdir` | directory holding the input `photometry.dat` / `photometry_results.txt` |
| `outdir` | directory for outputs (default: `workdir`) |
| `mc_seed` | RNG seed for the Monte-Carlo error propagation |

Outputs written to `outdir`: `photometry_results.txt` (fit parameter table with
confidence limits, in the original `print_struct` format) and
`photometry_results_stellar_c{i}.txt` (R/M/L per component).

Niche options of the original template (`mass_can`, `derive_logg`,
`hb_distance`, `derive_logg_c2`, `derive_sr`, `R1`) and black-body components
are not yet implemented; the config loader rejects them explicitly.

## Debugging aids

Two environment variables help diagnose parity issues (both read a file of
`name value` lines):

- `SEDFIT_EVAL=<file>`: after the fit, set those parameters, print the reduced
  χ² at that point (using the final dataset) alongside the best-fit reduced χ²,
  and exit. Use it to check whether another code's best fit is a lower point in
  this model (convergence deficiency) or equal (genuine degeneracy).
- `SEDFIT_INIT=<file>`: override starting parameter values after the initial
  guesses, to probe convergence-path sensitivity.

## Validation

```sh
python3 validation/compare.py            # all three reference stars
python3 validation/compare.py <gaiaid>   # one star
```

Compares `validation/cpp_<id>/` against the frozen reference fixtures in
`validation/ref_<id>/` (produced by the original ISIS script) with
physically-motivated tolerances documented at the top of the script.

## Layout

See [PLAN.md](PLAN.md) for the full file map and the parity-critical algorithm
details extracted from the S-Lang/ISIS sources.
