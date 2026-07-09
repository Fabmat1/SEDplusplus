# SED++ (`sedfit`)

A C++ reimplementation of the `stellar_isisscripts` photometric SED-fitting
pipeline (`templates/photometry.sl`), used to derive angular diameters, colour
excesses, atmospheric parameters, and stellar radii/masses/luminosities of hot
subdwarfs and similar objects from broad-band photometry.

Implemented in three validated phases:

- **Phase 1 — fit engine** (`PLAN.md`): mpfit-based SED fitting, confidence
  limits, prescribed-parameter error propagation, Monte-Carlo R/M/L.
- **Phase 2 — data acquisition** (`PHASE2.md`, `PHASE2_SPECTRA.md`): live
  photometry queries (VizieR bulk + TAP services), IRSA reddening, Gaia DR3
  astrometry (with Lindegren/El-Badry corrections), SIMBAD name resolution,
  IUE/MAST spectrophotometric box magnitudes.
- **Phase 3 — output products** (`PHASE3.md`): `photometry_fit*.txt`,
  `photometry_results.tex` (+ PDF via pdflatex), `SED_results.fits`, and the
  SED plot (`scripts/plot_sed.py`), all behind config toggles that default to
  off so bulk fitting stays lean.
- **Bulk local fitting** (`scripts/bulk_fit.py` + `scripts/sed_local.py`,
  campaign example `configs/bulk_campaign.json`): fully offline fitting of
  the local per-HEALPix SED parquet tables. Per pixel it joins the local
  Gaia store (corrected parallaxes, GSP-Phot priors) and the Wang et al.
  2025 3D dust map, classifies band columns via the photometry pipeline's
  `catalogs_use.txt`, and fits with an optional MS-first strategy: the
  main-sequence candidate runs first (Teff from GSP-Phot, ZAMS logg init);
  stars with a good, mass-radius-plausible MS fit skip the remaining
  candidates. `sedfit` reads `refdata/zeropoint_offsets.txt` when present
  (else its builtin table), accepts `filter_passbands.fits[.gz]` with
  either DES or DECam extension names, and caps pathological
  confidence-loop restarts via the `max_conf_restarts` config key (default
  1000 = ISIS-equivalent; bulk uses 10).

## Status

Validated to parity against the original pipeline for three Gaia DR3 reference
stars (two single-component sdB fits, one composite sdB + late-type fit) plus
13 additional hot subdwarfs for the query layer. Agreement:

- Fit parameter values: single-component fits to ~1e-7 relative; the composite
  fit to ≤1.3% of each parameter's 1σ confidence interval (its reddening is
  degenerate, pinned against the E=0 boundary).
- Confidence limits: to ~1e-5 of the interval width (single); ≤6% (composite).
- Query layer: every fit-relevant photometry row byte-identical to ISIS.
- Reduced χ² and generic excess noise: to the printed precision.
- MC-derived radius/mass/luminosity: to <1% (Monte-Carlo noise floor).
- Phase-3 outputs: text/FITS products at the same tolerances (structure
  byte-exact); TeX byte-identical outside Monte-Carlo digits; plots visually
  equivalent to the xfig originals.

## Installation

Requirements:

- a C++20 compiler and CMake ≥ 3.20;
- via pkg-config: `cfitsio`, `gsl`, `nlohmann_json`, `libcurl`, `zlib`
  (Arch: `pacman -S cfitsio gsl nlohmann-json curl zlib`;
  Debian/Ubuntu: `apt install libcfitsio-dev libgsl-dev nlohmann-json3-dev
  libcurl4-openssl-dev zlib1g-dev`);
- `mpfit` (the same MINPACK-derived Levenberg-Marquardt engine ISIS uses) is
  vendored in `extern/mpfit/` — nothing to install;
- model grids and the `stellar_isisscripts` `refdata` directory (filter
  curves), referenced from the config;
- optional, only for the corresponding output toggles: `pdflatex`
  (`write_tex`) and `python3` with `astropy` + `matplotlib` (`plot`).

```sh
cmake -B build
cmake --build build -j
```

The binary is `build/sedfit`; run it from the repository root (or set
`plot_script` in the config so the plot helper is found from elsewhere).

## Run

```sh
./build/sedfit config.json
```

Given only a star name or coordinates, the pipeline is self-sufficient: if
`workdir` has no cached `photometry.dat`, photometry (incl. IUE/MAST box
magnitudes) and IRSA reddening are queried live and cached there; astrometry
comes from the cached `photometry_results.txt` header, `fix_distance`, or a
live Gaia DR3 query, in that order; a missing `coordinates` field is resolved
via SIMBAD from `star`.

Standalone helper modes:

```sh
./build/sedfit --query <ra> <dec> <out.dat>   # photometry query only
./build/sedfit --query "HZ 44" <out.dat>
./build/sedfit --astrometry <ra> <dec>        # Gaia astrometry only
./build/sedfit --astrometry "Feige 34"
```

The JSON config mirrors the header block of `photometry.sl`. See
`configs/*.json` for the three validated examples. Key fields:

| field | meaning |
|---|---|
| `star` | object identifier; resolved via SIMBAD when no coordinates are available |
| `coordinates` | `{ra, dec}` in degrees |
| `fix_distance`, `fix_distance_err` | assume a distance (kpc) instead of Gaia astrometry |
| `par` | `{name, value, freeze}` — parameters to set (names may use `*` wildcards) |
| `par_full` | `{name, value, freeze, min, max}` — parameters to set with ranges |
| `griddirectories`, `bpaths` | model-grid subdirectories and base paths to search |
| `conf_level` | −1 (none), 0 (68%), 1 (90%), 2 (99%) |
| `apply_ZPO_corr` | apply empirical zero-point-offset corrections |
| `remove_outliers` | σ threshold for outlier rejection (0 disables) |
| `nMC` | Monte-Carlo trials for R/M/L |
| `refdata` | path to the `stellar_isisscripts/refdata` directory (filter curves) |
| `workdir` | directory holding (or receiving) the input `photometry.dat` / `photometry_results.txt` |
| `outdir` | directory for outputs (default: `workdir`) |
| `mc_seed` | RNG seed for the Monte-Carlo error propagation |

Optional output products (all default **off**; enable per config):

| toggle | products |
|---|---|
| `write_model` | `photometry_fit.txt` (best-fit spectrum), `photometry_fit_mag.txt` / `photometry_fit_col.txt` (per-filter fluxes and residuals) |
| `write_tex` | `photometry_results.tex` result table + PDF via pdflatex (skipped with a warning if absent) |
| `write_fits` | `SED_results.fits` (fit parameters, stellar parameters, and — with `write_model` — filter fluxes and model spectrum) |
| `save_MC` | raw Monte-Carlo R/M/L draws as `MC_c*` extensions in `SED_results.fits` (large!) |
| `plot` | SED figure `photometry_SED.pdf` via `scripts/plot_sed.py` (implies `write_fits` + `write_model`); `plot_script` overrides the script path |

Always written to `outdir`: `photometry_results.txt` (fit parameter table with
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
python3 validation/compare.py                       # fit results, all three reference stars
python3 validation/compare.py <gaiaid>              # one star
python3 validation/compare_fit_txt.py <ref> <cpp>   # photometry_fit*.txt products
python3 validation/compare_fits.py <ref> <cpp>      # SED_results.fits
```

Compares `validation/cpp_<id>/` against the frozen reference fixtures in
`validation/ref_<id>/` (produced by the original ISIS script) with
physically-motivated tolerances documented at the top of each script.

## Layout

See [PLAN.md](PLAN.md) for the full file map and the parity-critical algorithm
details extracted from the S-Lang/ISIS sources; [PHASE2.md](PHASE2.md),
[PHASE2_SPECTRA.md](PHASE2_SPECTRA.md), and [PHASE3.md](PHASE3.md) document the
query layer, the spectra subsystem, and the output products.
