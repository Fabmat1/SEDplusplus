# SED++ Phase 2 — Remote query layer (port of `query_photometry.sl`)

**Handoff document for Phase 2.** Phase 1 (the fit engine) is complete and
validated — see `PLAN.md`. Phase 2 adds the data-acquisition layer so `sedfit`
can produce its own `photometry.dat` from (ra, dec) instead of depending on a
cached file from the old pipeline.

## Goal of this phase

Port `query_photometry()` (and its helpers) from
`~/Projects/ISIS_install/src/stellar_isisscripts/` to C++ so that, given
coordinates, we reproduce the cached `validation/ref_*/photometry.dat` fixtures
(same rows, magnitudes, flags, angular distances, and the `meanSFD`/`meanSandF`
reddening header). Then feed that table straight into the existing fitter.

Network access is available in this environment; VizieR TAP, MAST, and IRSA all
respond, so validation is against **live** queries compared to the fixtures
(catalogs are largely stable for these three bright stars).

## Reference sources (S-Lang)

| concern | file |
|---|---|
| main orchestrator + flag logic | `src/photometry/query_photometry.sl` (2445 lines) |
| catalog registry (308 rows) | same file, `photometric_catalogs()` @ line 2062 |
| VizieR bulk query + VOTable parse | `src/miscellaneous/query_vizier_catalogue.sl` (`query_vizier_catalogues_votable`, `parse_vizier_votable_multi`) |
| IRSA reddening (SFD/SandF) | `src/miscellaneous/query_reddening.sl` (hits IRSA DUST cgi, XML parse) |
| table container + `.dat` write | `src/photometry/photometric_table.sl` |
| SIMBAD name resolve (not needed — we pass ra/dec) | `src/miscellaneous/resolve_simbad.sl` |

## Data flow (as the template drives it)

1. `query_photometry(ra, dec)` → collects every registry catalog whose
   `[dec_min, dec_max]` contains dec, issues **one** VizieR VOTable request for
   all VizieR catalogs at once (max radius; filtered per-catalog afterwards),
   plus separate TAP requests for the non-VizieR ones (MAST PS1_DR2, NOAO
   datalab, Skymapper, JPLUS…). Per-catalog flag postprocessing then runs.
2. Separately, the template calls `query_reddening(ra, dec)` → IRSA DUST →
   fills `interstellar_reddening` (meanSFD/stdSFD/meanSandF/stdSandF).
   `query_reddening()` inside `query_photometry` is **commented out** — it is a
   distinct step.
3. `photo.write("photometry.dat")` emits the fixture format (header
   `%.4f` reddening + fixed 8-col table).

## VizieR bulk query mechanics (port target)

URL: `http://vizier.u-strasbg.fr/viz-bin/votable?-source=<c1,c2,...>&-c=<ra>%20<+dec>&-c.rs=<maxRadiusArcsec>&-out.add=_r&-out.all`
(TAPVizieR sync endpoint `https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync`
also works and returns cleaner VOTable; decide during impl — the CGI votable
endpoint is what the S-Lang uses and groups tables by `ID="II_335_galex_ais"`).
Response is a multi-`TABLE` VOTable; per catalog:
- locate `TABLE` by `ID="<cat with / → _>"` (fallback `name="<cat>"`),
- `FIELD name="..."` → column names; `TABLEDATA`/`TR`/`TD` → string cells,
- compute `dist = 3600·sqrt((Δra·cos(dec))² + Δdec²)`, keep rows ≤ this
  catalog's radius (`max(search_radius, angular_accuracy/3600)`, arcsec),
- pick requested columns (filter_colname + error_colname), case-insensitive
  fallback; empty cell → NaN.

## Registry row schema (308 rows)

`catalogue  ra_colname  dec_colname  system  filter_colname  error_colname  passband  type  dec_min  dec_max  angular_accuracy`
Auto-extract into `src/catalog_table.inc` by parsing `photometric_catalogs()`
(do NOT hand-type; same approach as `zp_table.inc`/`f19_table.inc`). Commented
(`%`) rows are excluded. `error_colname == "NaN"` means no uncertainty column.

## Flag semantics (from the docstring, applied in the 721–2062 body)

```
-5 not useful / poorly calibrated      -1 serious catalog flags
-4 suspected blending (WISE in disk)    0 good
-3 outdated (superseded DR present)      1 redundant (same datum, 2 catalogs)
                                         2 replaced by an average
```
The per-catalog quirks (GALEX Wall+2019 bright correction, Gaia EDR3 G
correction, WISE disk `disk_b`/`disk_l` blending, DR redundancy dedup) live in
the main body and must be ported to match the fixtures. **This is the bulk of
task #4 and still needs a close read of lines 721–2062.**

## Proposed C++ module layout

| file | contents |
|---|---|
| `src/http.{hpp,cpp}` | libcurl RAII wrapper: GET/POST (form + raw), timeouts, global init, connectivity probe, min-interval throttle |
| `src/votable.{hpp,cpp}` | multi-table VOTable parser (by table ID), FIELD/TABLEDATA → rows of string cells; string-based like the S-Lang, tolerant of `<TD/>` |
| `src/catalog_registry.{hpp,cpp}` + `src/catalog_table.inc` | the 308-row registry + accessors, dec-range selection, per-catalog column list |
| `src/reddening.{hpp,cpp}` | IRSA DUST query → meanSFD/stdSFD/meanSandF/stdSandF |
| `src/query.{hpp,cpp}` | orchestrator: bulk VizieR, MAST PS1_DR2, flag postprocessing, angular distance, redundancy dedup → `PhotometryTable` |
| driver | `sedfit --query <config>`: if `photometry.dat` absent, query + write it, then fit (mirrors the template's `stat_file(...)==NULL` gate) |

Writer: extend `PhotometryTable` with a `write()` matching
`photometric_table.sl` (header `%.4f` reddening + the 8-col body; the reader
already exists in `photometry_table.cpp`).

## Scope of the FIRST validated slice

Target = reproduce the three `ref_*/photometry.dat` fixtures.
- **In:** VizieR bulk query (covers all reference-star catalogs except PS1_DR2),
  MAST PS1_DR2 TAP (star 1 uses it), IRSA reddening, per-catalog flag logic.
- **Deferred (no reference star needs them):** NOAO datalab (DES/DECaPS/SMASH/
  DELVE), Gaia XP box filters (`box_from_gaia`), IUE/MAST spectra downloads,
  SIMBAD name resolution (we always pass ra/dec).

## Status (updated 2026-07-03)

- [x] Study `query_photometry.sl` + helpers, confirm network + libcurl (8.21.0).
- [x] **http foundation** — `src/http.{hpp,cpp}`: libcurl RAII wrapper (GET/POST,
      timeouts, throttle, probe, url_encode). Live-verified.
- [x] **votable parser** — `src/votable.{hpp,cpp}`: multi-TABLE parser (by ID),
      FIELD/TABLEDATA → string rows, `<TD/>` + entity handling. Live-verified:
      a bulk VizieR query for star 1 returns all 6 test catalogs parsed cleanly.
- [x] **catalog registry** — `src/catalog_table.inc` (308 rows, auto-generated
      by `scratchpad/extract_catalogs.py`) + `src/catalog_registry.{hpp,cpp}`
      (parse, `is_vizier_catalogue`, dec-range selection). 76 catalogues,
      58 VizieR; 173 VizieR rows in-range at star 1's dec.
- [x] Confirmed: VizieR rounds `_r` to 3 dp, so angular distance must be
      **recomputed** from the ra/dec columns (matches fixture full precision),
      exactly as `parse_vizier_votable_multi` does.
- [ ] **reddening** (IRSA DUST → meanSFD/meanSandF) — `src/reddening.{hpp,cpp}`.
- [ ] **orchestrator + flag logic** — `src/query.{hpp,cpp}`. THE big remaining
      piece: still needs a close read of `query_photometry.sl` lines 721–2062
      (GALEX Wall+2019 correction, Gaia G correction, WISE disk blending, DR
      redundancy dedup, per-catalog flags). Registry column-collection per
      catalogue and the single bulk request are understood; the postprocessing
      is not yet ported.
- [ ] **MAST PS1_DR2 TAP** (star 1 needs PS1_DR2 rows; VizieR has II/349/ps1).
- [ ] **`.dat` writer** (extend `PhotometryTable`) + **`--query` driver** in
      `main.cpp` (gate on `photometry.dat` absent, like the template).
- [ ] **validate** produced `.dat` vs the three fixtures.

Temporary `qtest` target (`src/qtest.cpp`, `EXCLUDE_FROM_ALL`) exercises the
foundation on live data; remove once the driver is wired in.
