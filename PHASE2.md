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

## Status (COMPLETE — 2026-07-04)

The query layer is implemented, integrated into `sedfit`, and validated against
the three fixtures: **every fetchable fixture row reproduces byte-for-byte**
(flag, magnitude at full precision, uncertainty, angular distance, catalogue).

- [x] Study `query_photometry.sl` + helpers, confirm network + libcurl (8.21.0).
- [x] **http foundation** — `src/http.{hpp,cpp}`: libcurl wrapper (GET/POST,
      timeouts, throttle, probe, url_encode).
- [x] **votable parser** — `src/votable.{hpp,cpp}`: multi-TABLE parser (by ID),
      FIELD name+datatype, TABLEDATA → string rows, `<TD/>` + entity handling.
- [x] **catalog registry** — `src/catalog_table.inc` (308 rows) +
      `src/catalog_registry.{hpp,cpp}`.
- [x] **reddening** — `src/reddening.{hpp,cpp}` (IRSA DUST). meanSFD/stdSFD/
      meanSandF/stdSandF match all three fixture headers exactly.
- [x] **orchestrator + flag logic** — `src/query.{hpp,cpp}`. One bulk VizieR
      VOTable request for all in-dec VizieR catalogues, byte-sorted per-catalogue
      processing (matches `array_sort` order, which the `-3`/redundancy rules
      depend on), closest-in-radius selection (flat cut like the parser, then
      true `angular_separation`), and the full per-catalogue quality logic:
      2MASS Rflg/Cflg/Xflg/snr; Gaia `E(BP/RP)Corr` gate + saturation correction
      + `round2(·, round_err(e).digit-1)`; GALEX Wall+2019 (Eq. 5) + artefact
      bits; APASS→`-5` with 0.04 quadrature; SDSS `q_mode`/`Q` + `mode==1`
      primary filter; PS1-DR1 `flag=1` redundant + Qual bits; unWISE flux→mag +
      background-subtraction error; CatWISE `-3` when unWISE present; AllWISE
      always `-3`; IGAPS `errBits`; plus the flag==0 duplicate resolution.
      `round_err`/`round2` reused from `texval` (round2 kept in `util`, S-Lang
      multiply form, for bit-exact rounding). `radec2lb`/`angular_separation`
      ported for disk-blending + distances.
- [x] **non-VizieR TAP** — `tap_query` (POST REQUEST/LANG/FORMAT/QUERY to
      `<url>/sync`, parse VOTable): PS1_DR2 (MAST), DELVE_DR2 (NOAO datalab),
      GPS/LAS/GCS/DXS/UDS_DR11 (ROE WSA). MAST float32 columns re-quantised to
      their shortest float decimal (`Row::dstore`) so `%S` matches the fixture.
- [x] **`.dat` writer** — `PhotometryTable::write()` (reuses `slang_repr` +
      `write_table` from `results`; header `%.10f`/`%.4f`, `-nan` for masked).
- [x] **`--query` driver** — `sedfit --query <ra> <dec> <out.dat>` (standalone)
      and an auto-query gate when `photometry.dat` is absent and ra/dec known
      (mirrors the template's `stat_file(...)==NULL`).
- [x] **validate** — see below. `qtest` target removed; `src/qtest.cpp` deleted.

### Validation result (live queries vs `validation/ref_*/photometry.dat`)

All rows that can be fetched match byte-for-byte. The only diffs are external:

- **UKIDSS** (`GPS_DR11` star1 J/H/K, `LAS_DR11` star3 Y/J/H/K): the ROE WSA TAP
  (`http://tap.roe.ac.uk/wsa`) currently returns **zero rows for anonymous
  access even on an unfiltered `TOP 3`** — a service/auth state, not a query
  bug. The code path is implemented identically to the S-Lang and will return
  data when the service does.
- **star2 `DELVE_DR2`** (4 extra rows): the bright target is not in DELVE
  (saturated); a faint neighbour at ~12.4″ is now the closest source in the
  ±36″ box, so it is matched (i,z at flag 0, g null→`-1`). This is exactly what
  `query_photometry.sl` produces against today's DELVE catalogue; the fixture
  predates those faint detections (or datalab was unreachable at capture). NB:
  this would perturb a *fresh* star2 fit — it is data drift, not a port defect.

### Not ported (deferred; no fixture needs them)

Gaia XP / IUE / MAST box-filter spectra (`box_from_gaia`, `average_boxes`),
SIMBAD name resolution (we pass ra/dec), and the surveys absent from the three
fixtures (Skymapper, JPLUS, SPLUS, DES/DECaPS/SMASH, VISTA/VHS/VIKING/…, HSC,
Spitzer/SWIFT/XMM, …). The per-catalogue dispatch in `process_filter` is
structured so these slot in as additional branches.
