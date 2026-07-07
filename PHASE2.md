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

~~Gaia XP / IUE / MAST box-filter spectra~~ — **PORTED & byte-validated
2026-07-06** (`src/spectra.{hpp,cpp}`; enabled by default in `--query` and the
auto-query path, mirroring the template's `; IUE, MAST` call): all 33 box rows
of the HZ 44 fixture (`validation/spectra/HZ44_1473687671071803520/`)
reproduce byte-for-byte, FOS/Gaia-XP zero-row quirks replicated — see
[`PHASE2_SPECTRA.md`](PHASE2_SPECTRA.md) §10 for the implementation notes and
validation. SIMBAD name resolution was ported in Phase 4 (below). Still
absent: only the inactive-registry catalogues listed in the gap list.

## Status — Part 2: extended catalogue branches (validated 2026-07-06)

The branches added in the (previously uncommitted) `query.cpp` diff were validated
against the 13-star fixture set in `validation/newstars/` (ISIS `photometry.dat`)
using `validation/newstars/compare_query.py`. C++ `--query` output lives in
`validation/newstars_cpp/<id>/query.dat`.

**Result: every fit-relevant (flag-0) row that both sides can fetch matches
byte-for-byte on all 13 stars.** `good✓ == C++0` for every star (each C++ flag-0
row reproduces the ISIS magnitude exactly). Remaining diffs are all external
(service access / fixture-capture), never a port defect.

### Branches now ported + validated (new dispatch + TAP layouts)

- [x] **SkyMapper DR4** (`Skymapper_DR4`, TAP cone search, `api.skymapper.nci.org.au`).
- [x] **JPLUS DR3** (`JPLUS_DR3`, vector-cell TAP: 12 filters in one row via
      `process_jplus`, `mag_iso_worstpsf`/`mag_err_psfcor`/`flags`).
- [x] **DES DR1** (`II/357/des_dr1`, VizieR; always `-3`/`-1`, quality bits — the
      "DES redundancy" rule: DR1 is deprecated in favour of DR2, `flag=-3`).
- [x] **DES DR2 / DECaPS DR1 / SMASH DR2** (NOAO datalab TAP, `datalab.noirlab.edu`).
- [x] **HSC** (STScI HSCv3 TAP; always `-1`).
- [x] **XMM-OM** (`II/356/xmmom41s`), **Swift/UVOT** (`II/339/uvotssc1`).
- [x] **Spitzer** SAGE (`II/305/catalog`), GLIMPSE (`II/293/glimpse`), SpiKeS
      (`J/ApJS/254/11/spikes`), M31 (`J/MNRAS/459/1403/catalog`), SEIP
      (`II/368/sstsl2`, µJy→Vega).
- [x] **VISTA** VHS/VIKING/VMC/VVV/VIDEO (VSA TAP `wfaudata.roe.ac.uk`), **UHS DR3**
      (WSA TAP).
- [x] **S-PLUS DR4** (`II/380/splusdr4`), **ZTF** (`J/ApJS/249/18/table3`), **KiDS**
      (`II/347/kids_dr3`), **VST ATLAS** (`II/350/vstatlas`), **KIS**
      (`J/AJ/144/24/kisdr2`), **Edinburgh-Cape** (`J/MNRAS/*`), **BATC**
      (`II/262/batc`), **Geneva** (`II/169/main`), **Ducati** (`II/237/colors`),
      **TD1** (`II/59B/catalog`), **FAUST** (`J/ApJS/96/461/table2`), **DENIS**,
      **Kilkenny sdB** (`III/137/catalog`), **2MASS-in-GC** (`J/AJ/150/176/table3`),
      **LMC-PS** (`J/AJ/128/1606/lmcps`), **SPM4**/**Tycho-2** (`-5`).
- [x] **VizieR mirror fall-over** (`vizier_bulk`): cfa.harvard → u-strasbg →
      cds.unistra; byte-identical multi-TABLE documents, transport-only fall-over.

### Fixes made in this pass

1. **unWISE saturated-quality rule** (`II/363/unwise`, star `85707760314515712`).
   The S-Lang uses a *greedy* decomposition subtracting `2^(k-1)` for k=7…0 and
   flags at the `k==6` step — this is **not** a plain `Flags & 32` test. e.g.
   `FlagsW1 = 128` → −64 → 64 → (−32 at k==6) → `flag=-1`. The old `& 32L` mask
   gave `128 & 32 = 0` and left W1 at flag 0. Replaced with the faithful loop.
   (`src/query.cpp`, unWISE branch.) → star `85707` W1 now `-1`, matches ISIS.
2. **DES DR1 redundancy rule** (`II/357/des_dr1`, star `2533806305483837440`).
   Verified: DR1 defaults to `flag=-3` ("exists in DR2"), `-1` on bad quality
   bits, never flag 0 (the flag==0 DR2-supersession block is unreachable — same
   in S-Lang). C++ reproduces all 5 DR1 rows (`-3,-3,-1,-3,-1` for g,r,i,z,Y)
   byte-for-byte against the fixture. Rule already covered; confirmed correct.
3. **J/AcA/50/307/phot** (OGLE) — closed the one active-registry gap: flag V/B/I
   with <4 measurements (`Vdat`/`Bdat`/`Idat` ≤ 3). No test star reaches its
   dec range (−73…−64) so it is unexercised by this fixture set, but the branch
   now matches the S-Lang.

### Residual diffs (all external; not port defects)

| star | catalogue | kind | explanation |
|------|-----------|------|-------------|
| 2533806305483837440 | `LAS_DR11` (UKIDSS) | onlyISIS ×4, flag −1 | ROE WSA TAP returns 0 rows for anonymous access (known). Not fit-relevant (−1). |
| 595128265015393152 | `LAS_DR11` (UKIDSS) | onlyISIS ×4, **flag 0** | Same ROE WSA anonymous-access limit; here the ISIS rows are flag 0, so C++ is missing 4 fit-relevant bands until the service serves anon requests. |
| 682072009543970688 | `LAS_DR11` (UKIDSS) | onlyISIS ×4, **flag 0** | idem. |
| 3355897092143865088 | `GPS_DR11` (UKIDSS) | onlyISIS ×3, **flag 0** | idem (GPS WSA table). |
| 2533806305483837440 | `DES_DR2` | onlyC++ ×5, flag 0 | C++ correctly retrieves DES_DR2 (datalab `des_dr2.mag`) at 0.046″. **No ISIS fixture captured DES_DR2 at all** (that datalab table was unavailable when fixtures were generated, though DELVE_DR2 served fine and matches identically). C++ is more complete than the fixture; behaviour is faithful to the S-Lang (DES_DR2 has no supersession, stays flag 0). |
| 1196375273586921472 | `II/328/allwise` | onlyC++ ×4, flag −3/−1 | AllWISE matches byte-for-byte on the other 11/13 stars. Here the nearest AllWISE source is a 12.12″ neighbour (a different, brighter star) which C++ correctly matches and flags −3/−1; the ISIS fixture has no AllWISE row for this one position (transient empty AllWISE sub-table in ISIS's bulk VizieR fetch at capture). Not fit-relevant (−3/−1). |

All other stars: `onlyISIS = onlyC++ = flagDiff = 0`, `good✓ = ISIS0 = C++0`.

### Remaining gap list (S-Lang dispatch branches with no C++ equivalent)

Method: compared `phot_cat.catalogue[index]=="…"` branches in `query_photometry.sl`
against `cat == "…"` in `src/query.cpp` and against the active registry
(`catalog_table.inc`, 308 rows).

- **Closed this pass:** `J/AcA/50/307/phot` (was the only *active-registry*
  catalogue lacking its flag branch).
- **Inactive — no gap in output** (present as S-Lang branches but **not** in the
  active registry, i.e. `%`-commented out of `photometric_catalogs()`, so never
  queried): `SPLUS_DR2`, `SPLUS_DR3`, `SPLUS_DR4`, `II/319/las9`, `II/319/gcs9`,
  `II/319/dxs9`, `II/316/gps6`. If any are re-enabled in the registry they will
  need branches + (for SPLUS) their TAP layout.
- **Coverage confirmed:** every active-registry non-VizieR (TAP) catalogue has a
  `tap_url_for` entry (SkyMapper, JPLUS, PS1, DELVE/DES/DECaPS/SMASH datalab,
  UHS + GPS/LAS/GCS/DXS/UDS WSA, VHS/VIKING/VMC/VVV/VIDEO VSA). `HSC` has a C++
  branch + TAP url but is not in the active registry (harmless dead code).
- **Out of scope (unchanged):** Gaia astrometry, SIMBAD, Gaia XP / IUE / box
  spectra, output products.

---

## Phase 3: `query_astrometry` port (live Gaia DR3 astrometry)

Closes the last from-scratch gap: with only ra/dec, ISIS gets parallax etc.
LIVE from `query_astrometry`; the C++ previously required a cached
`photometry_results.txt` header. Ported into `src/astrometry.{hpp,cpp}`,
wired into `main.cpp` as a third precedence tier (cached header → fix_distance →
**live query**).

### Provenance (verified against `src/kinematics/query_astrometry.sl` + `templates/photometry.sl:210-249`)

**Call-site (photometry.sl:213-215)** — qualifiers `quality_warnings`,
`corrected_values`, `search_radius=4./60/60` (= 0.00111111 deg = 4″),
`columns="designation,parallax,parallax_error,ruwe,ipd_gof_harmonic_amplitude,ipd_frac_multi_peak,visibility_periods_used"`.
catalog defaults to `gaiadr3.gaia_source`.

**Service / endpoint.** `corrected_values`+DR3 makes the S-Lang probe
`https://gaia.ari.uni-heidelberg.de/tap/` (ARI Heidelberg, DaCHS). If it
answers it becomes the *sole* URL and the query gains
`INNER JOIN gaiadr3.gaia_source_corrections USING (source_id)` plus a
`parallax_corr` column (Lindegren-2020 zero-point applied server-side). Only if
Heidelberg is offline does it fall back to `https://gaia.aip.de/tap` then
`https://gea.esac.esa.int/tap-server/tap/` — those have **no** `parallax_corr`,
so the raw parallax/parallax_error are used (see cut logic below).

**Columns actually sent** (after the S-Lang's column-injection logic):
`designation,parallax,parallax_error,ruwe,ipd_gof_harmonic_amplitude,ipd_frac_multi_peak,visibility_periods_used,phot_g_mean_mag,parallax_corr`
(`phot_g_mean_mag` injected because `parallax_error` is present + corrected_values;
`parallax_corr` injected because Heidelberg online). quality_warnings adds
nothing new (its DR3 indicators ruwe/ipd_gof/ipd_frac_multi_peak/vpu already listed).

**ADQL** (radius r deg):
`SELECT <cols> FROM gaiadr3.gaia_source INNER JOIN gaiadr3.gaia_source_corrections USING (source_id) WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', <ra>, <dec>, <r> ))`
POST to `<url>/sync` with `REQUEST=doQuery&LANG=ADQL&FORMAT=votable/td&QUERY=…`.
NOTE: Heidelberg/DaCHS defaults to base64 `BINARY` serialization; must request
`FORMAT=votable/td` to get `<TABLEDATA>` (the shared votable parser is TABLEDATA-only).
Closest-row selection = min true angular separation among in-radius rows (closest_object).

**Correction formulas.**
- `parallax_corr` — taken verbatim from the Heidelberg corrections table
  (Lindegren et al. 2020 zero-point; not recomputed in C++).
- `parallax_error_corr` (El-Badry, Rix & Heintz 2021 Eq. 16, `query_astrometry.sl:214`):
  `parallax_error * (0.21*exp(-((G-12.65)/0.9)^2) + 1.141 + 0.0040*G - 0.00062*G^2)`,
  G = `phot_g_mean_mag`. Computed locally in C++.
- pmra/pmdec corrections (Cantat-Gaudin & Brandt 2021) — **not applicable** here
  (no pmra/pmdec columns requested), so not ported.

**NaN-masking / final value selection (photometry.sl:229-248).** After the query:
1. if raw `parallax<=0` → parallax=NaN; if raw `parallax/parallax_error<1` → parallax=NaN.
2. if both `parallax_corr` and `parallax_error_corr` exist:
   - if `parallax_corr<=0` → parallax_corr=NaN; if `parallax_corr/parallax_error_corr<1` → parallax_corr=NaN;
   - if parallax_corr valid&>0 AND parallax_error_corr valid&>0 →
     **use parallax=parallax_corr, parallax_error=parallax_error_corr** (this is the
     normal Heidelberg path; the corrected pair is used *together* or not at all).
3. `good_astrometry = parallax not-NaN and >0`.
Header line-3 fields ruwe/ipd_gof/vpu/parallax/parallax_error are printed with `%.3g`.

**Verified reproduction (pre-implementation, curl):**
- ref_3431631727942970240 (90.1291,+29.1486): raw plx 1.11287, plx_corr 1.14606→`1.15`;
  ruwe 0.99195→`0.992`; ipd_gof 0.12586→`0.126`; vpu `16`; G 13.167337;
  plx_err 0.0290380 → corr 0.035922→`0.0359`. All match fixture.
- 800817856594433152 (144.0471,+38.3324): plx_corr 0.27112→`0.271`; ruwe `0.967`;
  ipd_gof `0.0525`; vpu `21`; plx_err 0.0254725 → corr 0.027131→`0.0271`. All match fixture.

### Implementation (done, validated 2026-07-06)

- **`src/astrometry.{hpp,cpp}`** — `query_astrometry(ra, dec, radius=4/3600)`:
  probes Heidelberg (`http::probe`, 4s/6s); if up, single-URL query with the
  corrections JOIN, else AIP→ESA fallback without it. POSTs to `<url>/sync`
  with `FORMAT=votable/td` (DaCHS defaults to base64 BINARY, which the shared
  TABLEDATA-only parser can't read — this is the one deliberate deviation from
  query.cpp's `tap_query`, hence a local `tap_query_td`). Closest row by true
  angular separation; quality warnings printed (ruwe>1.4, ipd_gof>0.1,
  ipd_frac_multi_peak>10); EDR3 error inflation computed locally; NaN-masking +
  corrected-pair substitution per photometry.sl:216,229-248 done inside, so the
  returned `parallax`/`parallax_error` are final. Returns `good` flag.
- **`src/main.cpp`** — third astrometry tier after cached header and
  fix_distance: if ra/dec known, live query; `a.ra/a.dec` stay the *input*
  coordinates (header `# RA/DEC` are the requested position, matching ISIS).
  The pre-existing generic parallax masking downstream is idempotent on the
  already-masked live values. Also added standalone debug mode
  `sedfit --astrometry <ra> <dec>` printing the %.3g header-line-3 fields.
- **`CMakeLists.txt`** — added `src/astrometry.cpp`.

### 16-star parity (live query vs fixture header line 3, byte-exact compare)

Compared `ruwe; ipd_gof_harmonic_amplitude; visibility_periods_used; parallax;
parallax_error` as printed (%.3g) — full line string equality:

| star | result |
|------|--------|
| ref_3431631727942970240, ref_648029239761019776, ref_872072016870129024 | PASS |
| 2533806305483837440, 3355897092143865088, 865191719781344384, 595128265015393152 | PASS |
| 4028265945931548544, 363795760175504768, 800817856594433152, 1842216206036695296 | PASS |
| 4450708930484104960, 85707760314515712, 919963474904219776, 682072009543970688, 1196375273586921472 | PASS |

16/16 byte-exact. (Ref-star coordinates taken from fixture line 1, newstars
from `validation/newstars/stars.tsv`.)

### From-scratch smoke test

Fresh scratch workdir/outdir, config copy of `configs/new_800817856594433152.json`,
no cached files: run auto-queried astrometry then photometry, fit completed,
`photometry_results.txt` header lines 1-3 byte-identical to the fixture,
`dummy_1` (parallax) identical at full precision (0.2711215662559987 —
confirming ISIS stored the same Heidelberg `parallax_corr`), R/M/L
(dummy_2..4) populated and matching the fixture to ~4e-6 relative (MC noise).

### Caveats

- The Heidelberg corrections service is a hard dependency for *corrected*
  parallaxes (same as ISIS): if it is down, the AIP/ESA fallback returns raw
  parallax/parallax_error (no Lindegren-2020 zero point, no EdR3 error
  inflation) — identical behaviour to the S-Lang, but values then differ from
  fixtures generated while Heidelberg was up.
- `parallax_corr` is taken verbatim from `gaiadr3.gaia_source_corrections`
  (server-side Lindegren zero point), not recomputed — exactly like ISIS.
- pmra/pmdec (Cantat-Gaudin & Brandt 2021) corrections not ported (not
  requested by the photometry template).

## Phase 4: SIMBAD name resolution (port of `resolve_simbad.sl`)

### S-Lang semantics (from `src/miscellaneous/resolve_simbad.sl` + template use)

`resolve_simbad(String name)` -> Struct with fields `oid, ra, dec, main_id,
coordreference, nbreferences, parallax, radialvelocity`, or `NULL` on no match.

- Service: TAP sync endpoint, primary `http://simbad.u-strasbg.fr/simbad/sim-tap/sync`.
  A reachability probe (`curl -s --connect-timeout 4 -m 1 <url>`; offline iff
  exit code != 0) selects the Harvard mirror
  `http://simbad.cfa.harvard.edu/simbad/sim-tap/sync` when the primary is down.
  Only one URL is used per call (no cross-mirror retry on empty result).
- POST params: `REQUEST=doQuery PHASE=RUN FORMAT=csv LANG=ADQL` and
  url-encoded `query=`. ADQL:
  `SELECT basic.OID, RA, DEC, main_id, coo_bibcode AS "CoordReference",
   nbref AS "NbReferences", plx_value AS "Parallax",
   rvz_radvel AS "RadialVelocity" FROM basic JOIN ident ON oidref = oid
   WHERE id = '<name>';`
- Parsing: read 2 lines (header + first data row) via CSV decoder. If the query
  returns nothing (0 lines) -> print message + return NULL. If there is a header
  but no data row (unresolvable name) -> the `rcv[1]` access / CSV decode throws
  and the `catch` returns NULL. So **unresolvable name => NULL (not an error)**;
  the caller decides. Typecasts: ra/dec/parallax/radialvelocity via atof,
  oid/nbreferences via atoi.
- Template (`templates/photometry.sl`) name handling before resolution:
  - line 15: `star = strreplace(strtrim(star), "_", " ")` (underscores -> spaces).
  - line 13 (CLI arg only): if `_slang_guess_type(star)==Integer_Type`, i.e. a
    bare integer such as a Gaia source_id, prefix `"Gaia DR3 "`.
  - lines 197-208: coordinates are resolved via SIMBAD **only if still NULL**
    after the cached-header / explicit-coordinate paths; a NULL result throws
    `DataError "Object '<star>' could not be resolved and no coordinates were
    provided either."`.
  (`coord_from_identifier` / SALT-J parsing is a *separate* code path, not part
  of `resolve_simbad`; not ported here.)

### Implementation (done)

`src/simbad.{hpp,cpp}` — `normalize_star_name()` (trim, bare-integer ->
"Gaia DR3 <id>", underscore -> space) and `resolve_simbad()` returning
`std::optional<SimbadResult>`. Mirrors astrometry.cpp: reuses `sed::http` +
the VOTable parser (requests `FORMAT=votable` instead of the S-Lang CSV path,
same as the rest of the port), probe-selects primary vs Harvard mirror, single
ADQL query, `nullopt` on empty/unreachable. ADQL single-quotes in the name are
doubled for safety.

Wiring (numeric-coordinate paths left byte-identical):
- `sedfit --astrometry <ra> <dec>` unchanged; `sedfit --astrometry <name>`
  (argv not two parseable numbers) resolves the name first.
- `sedfit --query <ra> <dec> <out>` unchanged; `sedfit --query <name> <out>`
  resolves the name first.
- Config: after the cached-`photometry_results.txt` header block in main.cpp,
  if `cfg.ra`/`cfg.dec` are still unset (no `coordinates` in JSON, no cached
  header) the `star` field is resolved via SIMBAD; failure throws the same
  DataError-style message as the template. Configs that carry `coordinates`
  never hit the network here.

Gotcha found + fixed during validation: SIMBAD's sim-tap `sync` endpoint returns
base64 `<BINARY>` VOTable for `FORMAT=votable` (our parser is TABLEDATA-only), so
the request must ask for `FORMAT=votable/td` (url-encoded `votable%2Ftd`) — same
serialization issue as the Heidelberg astrometry service. The reachability probe
correctly treats the endpoint's bare-GET HTTP 400 as "reachable" (any status).

### Validation (live, 2026-07-06)

- Resolution of `HZ 44`, `Feige 34`, `BD+28 4211` all succeed; e.g. HZ 44 ->
  RA=200.89692893566996 DEC=36.13320670564 (matches SIMBAD published position).
- `sedfit --astrometry "HZ 44"` produces a byte-identical astrometry line to
  `sedfit --astrometry 200.89692893566996 36.13320670564` (resolved coords feed
  the downstream Gaia query correctly). Same object via the bare-integer Gaia
  source_id form `sedfit --astrometry 1473687671071803520` (normalized to
  "Gaia DR3 1473687671071803520") and the underscore form `Feige_34` -> same
  results.
- `sedfit --query "HZ 44" out.dat` writes `# RA = 200.8969289357 DEC = 36.1332067056`
  — the resolved coordinates flow into the photometry + reddening queries.
- Config path: a copy of `configs/new_1196375273586921472.json` with
  `coordinates` removed and an empty `workdir` resolves
  `Gaia DR3 1196375273586921472` -> RA=237.701242 DEC=16.440366 (PG 1548+166),
  matching the config's stored coordinates to ~0.01 arcsec, then proceeds to the
  astrometry+photometry queries.
- Unresolvable name (`ZZ NoSuchObject 99999`): CLI prints the DataError-style
  message and exits 1; config path throws the same message — matching the
  S-Lang's NULL-then-`throw DataError` behaviour (a miss is not an HTTP error).
