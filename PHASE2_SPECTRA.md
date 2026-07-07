# SED++ Phase 2 (addendum) — Spectra → box-filter pseudo-photometry

**STATUS: IMPLEMENTED & VALIDATED (2026-07-06) — see §10 at the bottom for
the implementation notes, validation results, and confirmed quirks.**
Code: `src/spectra.{hpp,cpp}` (+ hooks in `src/query.{hpp,cpp}`,
`src/main.cpp`, `src/fitsio_util.{hpp,cpp}`, `CMakeLists.txt`).

This document specifies the *spectra subsystem*: turning Gaia DR3 XP-sampled,
IUE, and MAST (FUSE/FOS) spectra into **box-filter** magnitude rows in
`photometry.dat`. Everything else in the query layer (all VizieR/TAP
catalogs, reddening, Gaia astrometry, SIMBAD) was already ported &
byte-validated.

All S-Lang lives in **one** file:
`~/Projects/ISIS_install/src/stellar_isisscripts/src/photometry/query_photometry.sl`
(helpers in `src/photometry/photometric_magnitudes.sl`,
`src/photometry/filter_passbands.sl`, `src/photometry/photometric_table.sl`,
`src/isisscripts/{round_err,round2,integrate_trapez}.sl`).

Reference fixture: `validation/spectra/HZ44_1473687671071803520/` (see §7).

---

## 0. Where it plugs into the pipeline / what enables it

Inside `query_photometry()` (`query_photometry.sl:401`), *after* the catalog
loop ends (line 1849), the "Construct box filters" section runs (lines
1850–2052):

| feature | enable condition | file:line |
|---|---|---|
| Gaia DR3 XP-sampled boxes | **always run** (unconditional) | `query_photometry.sl:1854-1859` |
| IUE boxes (VI/110/inescat) | `qualifier_exists("IUE")` | `:1862-1941` |
| MAST FUSE/FOS boxes | `qualifier_exists("MAST")` | `:1942-2052` |

The template drives it with **both** qualifiers on:
`validation/newstars/<id>/photometry.sl:191`:
```slang
photo = query_photometry(coordinates.ra, coordinates.dec; IUE, MAST, search_radius=4./60/60);
```
So the "standard settings" = `IUE, MAST, search_radius=4/3600 deg`. Gaia XP is
attempted for every object regardless of qualifiers.

Final step (line 2054): entries are re-sorted by `system` string
(`struct_filter(...array_sort(...system))`). Reddening query at line 2057 is
commented out (done separately by the template, see PHASE2.md).

### CRITICAL gotcha — Gaia XP is effectively OFF in practice
`query_gaia_spectrum()` (line 1) does the XP retrieval. Two facts make it
almost always yield `NULL` (⇒ **no** `Gaia_DR3_XP` rows), which is why none of
the 13 `validation/newstars/*/photometry.dat` fixtures contain XP boxes:

1. **Mirror ordering + guard** (lines 85-122): it probes TAP mirrors in order
   `["https://gaia.aip.de/tap", "https://gea.esac.esa.int/tap-server/tap"]`,
   picks the *first reachable* as `url_tap`. Then line 121:
   `if(string_match(url_tap,"gea.esac.esa.int")==0) return NULL;` — i.e. if the
   AIP mirror answered (it usually does), the function bails and returns NULL.
2. **Broken download URL** (lines 125-127): even when `url_tap` == esac, the
   datalink URL is built as `url_tap + "data?..."` =
   `https://gea.esac.esa.int/tap-server/tapdata?...` which currently **404s**
   (verified 2026-07-06). The working endpoint is
   `https://gea.esac.esa.int/data-server/data?ID=Gaia+DR3+<id>&RETRIEVAL_TYPE=XP_SAMPLED&format=FITS`
   (commented-out `url` at line 124). A 404 → HTML → `fits_read_table` throws
   `FitsError` → `s = NULL` (line 180).

**Implication for the C++ port:** to match the *current ISIS behavior* byte-for-
byte, `box_from_gaia` produces **nothing** (matches fixtures). The XP algorithm
is documented in §2 in case a future ISIS fix (or the C++ author deliberately)
enables it. The reference fixture in §7 therefore is authoritative: reproduce
exactly the rows ISIS actually emits. **Confirm against the fixture which
`Gaia_DR3_XP` rows (if any) are present before implementing XP.**

---

## 1. Row / file format (target output)

`photo.write("photometry.dat")` →
`write_photometric_table` (`photometric_table.sl:253`). Header lines:
```
# RA = %.10f DEC = %.10f
# meanSFD = %.4f stdSFD = %.4f
# meanSandF = %.4f stdSandF = %.4f
```
then `print_struct(fp, photometric_entries)` dumps the struct as whitespace-
aligned columns in this field order (set in `photometric_entry`,
`photometric_table.sl:52-61`):

```
flag  system  passband  magnitude  uncertainty  type  angu_dist_arcsec  VizieR_catalog
```

Box rows are ordinary rows with **`system == "box"`** and **`type ==
"magnitude"`**. The `passband` is the box string `"<wmin>_<wmax>"` (e.g.
`"1300_1800"`). `VizieR_catalog` is the provenance string (see each source
below). Numbers print in S-Lang default `print_struct` precision (magnitude &
uncertainty as `%S` of the Double — i.e. full repr like `13.481` or
`0.002857`; angu_dist_arcsec full double repr). Match the C++ writer already
used for the catalog rows — box rows go through the same `add_entry`/`write`.

Every box entry is created via
`photometric_entry(flag, "box", box, magnitude, uncertainty, "magnitude", angu_dist_arcsec, catalog_name)`
then `table.add_entry(...)`.

---

## 2. Gaia DR3 XP-sampled boxes — `box_from_gaia` (`:184-275`)

Driver (`:1854-1859`):
```slang
if(qualifier_exists("gaia_id") && typeof(qualifier("gaia_id"))==String_Type)
  sgaia = query_gaia_spectrum(qualifier("gaia_id"); id);
else
  sgaia = query_gaia_spectrum(table.ra, table.dec);
table = box_from_gaia(table, sgaia);
```

### 2a. `query_gaia_spectrum` (`:1`) — acquisition
- Resolve DR3 source_id from (ra,dec): TAP `gaiadr3.gaia_source`, cone radius
  `0.01` deg, columns `designation,ra,dec`, keep closest match, `id =
  strtok(designation)[-1]` (numeric id).
- Download (the *intended* URL; note §0 caveat):
  `GET https://gea.esac.esa.int/data-server/data?ID=Gaia+DR3+<id>&RETRIEVAL_TYPE=XP_SAMPLED&format=FITS`
  → FITS binary table. Columns used: `wavelength`, `flux`, `flux_error`.
- **Unit / error post-processing (must replicate if XP enabled)** (`:132-177`):
  - `wavelength *= 10.` (nm→Å).
  - `flux *= 100.; flux_error *= 100.` (W m⁻² nm⁻¹ → erg cm⁻² s⁻¹ Å⁻¹).
  - `flux_error *= 1.30` (blanket +30%).
  - u-band bumps (add fractions of the *median positive error* in each range):
    - `wavelength < 3800`: `flux_error += 0.5 * median(flux_error[good])`
    - `wavelength < 3660`: `flux_error += 0.6 * median(...)`
    - `wavelength < 3500`: `flux_error += 0.6 * median(...)`
    (each computed over its own sub-range, `good` = flux_error>0, applied in
    sequence so the <3500 region receives all three additions.)

### 2b. `box_from_gaia` (`:184`) — box magnitudes
1. `struct_filter(sgaia, where(flux>0 and flux_error>0))` (`:229`).
2. Standard boxes (`:230-234`), 14 of them:
   ```
   3400_3660 3660_4000 4000_4300 4300_4600 4600_4740 4740_5200
   5200_5600 5600_6000 6000_6400 6400_6800 6800_8400
   8400_8800 8800_9200 9200_10200
   ```
   (overridable via `boxes` qualifier; template does not override).
3. `photo = photometric_magnitudes(wavelength, flux; boxes=boxes)` → magnitudes
   (see §5 for the exact box-magnitude formula).
4. For each box where `magnitude` is not NaN:
   - `sscanf(box,"%F_%F",&l_start,&l_end)`.
   - `ind = where(l_start < wavelength < l_end)` (strict, endpoints excluded).
   - `u = union(l_start, wavelength[ind], l_end)` (sorted unique grid incl. box
     edges).
   - `int = gsl->interp_linear_integ(u, interpol(u, wavelength, flux), u[0], u[-1])`
     (linear-interp trapezoid integral; fallback `integrate_trapez(u, interpol(...))`).
   - **Error** (Gaussian propagation of trapezoid, `:262`):
     `int_error = sqrt( sum( (0.5*(u[[2:]]-u[[:-3]]) * flux_error[ind])^2 ) )`.
     Note: `u[[2:]]-u[[:-3]]` = for interior grid point k, `u[k+1]-u[k-1]`;
     length matches `flux_error[ind]` (the interior points).
   - `mag_error = round_err( abs( -2.5/ln(10) / int * int_error ) )` (§6).
   - `flag = 0`; **if `l_end < 3700` then `flag = -1`** (`:268`).
   - `catalog_name = "Gaia_DR3_XP"`, `angu_dist_arcsec = 0.` (`:265-266`).
   - `add_entry(flag,"box",box, round2(magnitude, mag_error.digit-1),
     mag_error.value, "magnitude", 0., "Gaia_DR3_XP")`.
   No `average_boxes` for Gaia XP (single spectrum).

`boxes_from_spectrum` (`:277`) is an identical helper for arbitrary CALSPEC-
style spectra; **not used by the pipeline** — skip unless needed.

---

## 3. IUE boxes (`:1862-1941`) — enabled by `IUE` qualifier

### 3a. Catalog query
`radius_iue = max(radius, 70/3600 deg)` (70″; `force_search_radius` ⇒ use
`radius`). Query VizieR **VI/110/inescat** (INES merged IUE log) via
`query_vizier_catalogue(ra, dec, radius_iue, "VI/110/inescat";
ra_colname="RAJ2000", dec_colname="DEJ2000")`. Columns used from result:
`Object, Camera, Image, Disp, Aper, AbnCode, RAJ2000, DEJ2000`.

Live TAP equivalent (verified working 2026-07-06):
```
https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync
  request=doQuery lang=ADQL format=csv
  query=SELECT * FROM "VI/110/inescat"
        WHERE 1=CONTAINS(POINT('ICRS',RAJ2000,DEJ2000),CIRCLE('ICRS',<ra>,<dec>,<radius_deg>))
```

### 3b. Per-spectrum flag from AbnCode/Aper (`:1874-1877`)
`flag = 0`; set `flag = -1` if `AbnCode` contains any of `A B C R T Z`
(`string_match` != 0), **or** `Aper == "S"` (small aperture).
(Ref: http://cdsarc.u-strasbg.fr/viz-bin/ReadMe/VI/110)

### 3c. Filename + download (`:1878-1887`)
```
filename = sprintf("%s%05.0f%s%s.FITS", Camera, Image, strreplace(Disp,"H","R"), Aper)
```
i.e. `Camera` + zero-padded 5-digit `Image` + `Disp` with **H→R** + `Aper`
(e.g. Camera=SWP,Image=3432,Disp=L,Aper=L → `SWP03432LL.FITS`;
high-disp Disp=H → `R`). Download into `./IUE/` subdir:
```
mkdir -p -m 0750 <cwd>IUE
# only if <cwd>IUE/<filename>.gz absent:
wget -q -O <cwd>IUE/<filename> "http://sdc.cab.inta-csic.es/ines/jsp/SingleDownload.jsp?filename=<filename>"
gzip -f <cwd>IUE/<filename>          # ⇒ <filename>.gz
```
Then `filename += ".gz"`. (Endpoint verified live 2026-07-06; returns a FITS
BinTable with cols `WAVELENGTH, FLUX, SIGMA, QUALITY`.) Corrupt files are
`try/catch`-skipped (`:1890`).

### 3d. Box magnitudes (`:1893-1932`)
- `angu_dist_arcsec = 3600 * angular_separation(ra, dec, RAJ2000[i], DEJ2000[i])`.
- **IUE boxes** (`:1895`): `["1300_1800", "2000_2500", "2500_3000"]`.
- Bad-pixel filter (`:1896`): `struct_filter(IUE, where(flux>0 and quality==0))`
  (quality flag must be exactly 0; see
  https://archive.stsci.edu/iue/manual/dacguide/node60.html).
- `photo = photometric_magnitudes(IUE.wavelength, IUE.flux; boxes=boxes)`.
- Per box (same integral recipe as §2b but using `IUE.sigma` for errors):
  - `int` via `gsl->interp_linear_integ(u, interpol(u,wavelength,flux),u[0],u[-1])`.
  - `int_error = sqrt(sum((0.5*(u[[2:]]-u[[:-3]])*IUE.sigma[ind])^2))`.
  - `mag_error = abs(-2.5/ln(10)/int*int_error)`;
    **then add generic uncertainty**: `mag_error = sqrt(mag_error^2 + 0.04^2)`
    (`:1928`); `mag_error = round_err(mag_error)`.
  - `add_entry(flag, "box", box, round2(magnitude, mag_error.digit-1),
    mag_error.value, "magnitude", angu_dist_arcsec,
    strreplace(filename,".FITS.gz","") + "_VI/110/inescat")`.
    ⇒ provenance like `SWP03432LL_VI/110/inescat`.

### 3e. Averaging (`:1935-1938`)
After all IUE spectra processed:
```
thres_iue = 0.32
table = average_boxes(table, boxes, thres_iue; mag_err_sys=0.03)
```
`boxes` here is the IUE box list. See §4 for `average_boxes`.

---

## 4. `average_boxes` (`:354-386`) — collapse multiple spectra per box

For each `box` in the given box list:
- `ind = where(passband==box and flag==0)` over **all current table entries**
  (only good, unaveraged rows).
- `ndata = length(ind)`. If `ndata > 1`:
  - **Mark all those rows `flag = 2`** ("replaced by an average"; they stay in
    the table).
  - `mag_box = magnitude[ind]`, sort ascending (`idx_sort = array_sort`).
  - Take the brightest fraction: `idx_cut = [0 : max(nint(ndata*thres), 1)]`
    (**inclusive** range ⇒ `max(nint(ndata*thres),1)+1` elements; always ≥2).
  - `mom = moment(mag_box[idx_cut])`; `mag_error = sqrt(mom.sdev^2 +
    mag_err_sys^2)` then `round_err`.
  - Add averaged row:
    `add_entry(0, "box", box, round2(mom.ave, mag_error.digit-1),
    mag_error.value, "magnitude",
    mean(angu_dist_arcsec[ind][idx_sort][idx_cut]), "Average")`.
  - `mag_err_sys` default 0.005; IUE passes 0.03, MAST uses default.
  Notes: `moment().sdev` is the sample std-dev (N-1). `nint` = round-half-to-
  even? (S-Lang `nint` rounds to nearest, .5 away from zero — verify against
  fixture). If `ndata <= 1` the single row keeps flag 0 (no Average row).

---

## 5. Box magnitude core — `photometric_magnitudes` + `filter_passbands`

The *magnitude value* for each box comes from `photometric_magnitudes(l, f;
boxes=...)` (`photometric_magnitudes.sl:82`), **not** from the local `int` in
§2b/§3d (that local `int` is only used for the error scaling). Recipe:

1. Box handling (`photometric_magnitudes.sl:686-716`): for each box string,
   appends an entry with `system="box"`, `type="magnitude"`, `passband=box`,
   `ZP=0`, `ZP_err=0.015`.
2. Filter shape (`filter_passbands.sl:120-133`): for `system=="box"`,
   `sscanf(passband,"%F_%F",&l_start,&l_end)`, response is a 2-point top-hat
   ```
   l = [l_start, l_end];  f = [1,1] / (1/l_start - 1/l_end)
   ```
   normalized so ∫ f d(−1/λ) = 1.
3. Integration (`photometric_magnitudes.sl:786-813`):
   - `ind = where(l_filter[0] <= l <= l_filter[-1])`, padded by one point each
     side to avoid extrapolation.
   - `u = union(l_filter, l[ind])`; `fu_filter = interpol(u, l_filter, f_filter)`
     (top-hat resampled); `fu = interpol(u, l[ind], f[ind])` (flux resampled).
   - `int = gsl->interp_linear_integ(u, fu*fu_filter, u[0], u[-1])`
     (fallback `integrate_trapez`).
   - `magnitude = -2.5*log10(int) + int_ref`, with
     **`int_ref = -2.407948242680184`** (= `2.5*log10(c·1e8·10^(48.6/-2.5))`,
     the AB zero-point constant). (Stroemgren-Hbeta exception does not apply to
     boxes.)
   - Equivalent closed form (see `magnitude_from_box`,
     `photometric_magnitudes.sl:40-80`):
     `norm = 1/(1/wmin - 1/wmax)`;
     `mag = -2.5*log10( norm * ∫flux dλ over [wmin,wmax] ) - 2.407948242680184`.
   - NaN if the box is outside the spectrum's wavelength coverage (`ind` empty
     or filter edges outside `l` range) ⇒ that box produces no row.

The ZP/ZP_err (0/0.015) on box entries are **not** written to photometry.dat
(only used inside the fitter later); the `.dat` uncertainty column is the
`mag_error.value` computed in §2b/§3d/§4.

---

## 6. Rounding helpers (exact — must match byte-for-byte)

`round_err(x)` (`isisscripts/round_err.sl`): rounds error per DIN 1333.
- `digits = floor(log10(x))`; if `round2(x/10^digits, -14) < 3` then `digits--`.
- `x = round(x*10^(-digits+2))/10^(-digits+2)`; `x = ceil(x*10^(-digits))/10^(-digits)`;
  then reformat via `sprintf("%.<−digits>f")` (if digits<0) or `%.0f`.
- Returns `{value=<rounded err>, digit=digits}`. `digit` is the decimal place;
  the magnitude is then rounded to `digit-1` (one extra digit than the error).

`round2(val, dig)` (`isisscripts/round2.sl`): `round(val*10^(-dig))*10^dig`.
Used as `round2(magnitude, mag_error.digit-1)`.

`integrate_trapez(x,y)` (`isisscripts/integrate_trapez.sl`):
`sum(0.5*(y[:-2]+y[1:])*(x[1:]-x[:-2]))`. (GSL `interp_linear_integ` gives the
same for linear interpolation; port either — they agree to ~1e-15.)

---

## 7. MAST FUSE/FOS boxes (`:1942-2052`) — enabled by `MAST` qualifier

(Lower priority — HZ 44 fixture may or may not have MAST rows; document for
completeness.)

- `radius_mast = max(radius, 60/3600 deg)` (60″).
- SSAP query: `ssapquery("https://archive.stsci.edu/ssap/search2.php", ra, dec,
  radius_mast)` (`src/miscellaneous/ssapquery.sl`). Result fields used:
  `url, telescop, instrume, aperture, cr_ident, fluxcal, der_snr, object,
  ang_sep, spectralaxisname, fluxaxisname`.
- Filters (`:1953-1959`): keep `telescop=="FUSE" or instrume=="FOS"`; then keep
  `der_snr>=2 and object!="WAVE" and aperture!="RFPT"`.
- Flag −1 (`:1967-1972`) if `fluxcal != "ABSOLUTE"`, or aperture
  `"1.25x20"`/`"4x20"` (FUSE HRS/MRS), or aperture contains "unknown".
- Filename: `sprintf("%s-%s-%s_%s.FITS", telescop, instrume, aperture,
  cr_ident)` with spaces stripped; download via `wget` from `query.url[i]` into
  `./MAST/`, then `gzip`.
- Read FITS BinTable; `filetype = header key FILETYPE` (skip `"WAV"`). Spectral
  axis = `get_struct_field(mspec, strlow(spectralaxisname))`, flux axis =
  `strlow(fluxaxisname)`, `sigma = mspec.sigma`. Filter `flux>0, sigma>0,
  wavelength>920` (drop EUV).
- Boxes:
  - FUSE (`:2005`): `["950_1000","1000_1080","1091_1120"]`, **flag forced −1**
    (unreliable pointing).
  - FOS (`:2011-2017`): 4-point grid `bbox=[wmin+15 : wmax-15 : #4]`, boxes =
    consecutive `"%.0f_%.0f"` pairs (3 boxes).
- Box magnitudes: same integral recipe as §2b (uses `mspec.sigma`), **no**
  generic-uncertainty addition, `mag_error = round_err(abs(-2.5/ln10/int*int_error))`.
  Provenance `<filename without .FITS.gz>_MAST`.
- Averaging (`:2047-2049`): `thres_mast = 0.16`, `average_boxes(table, boxes,
  thres_mast)` (default `mag_err_sys=0.005`). NOTE: `boxes` here is whatever the
  last spectrum set — a latent S-Lang quirk; verify against fixture which boxes
  actually get averaged.

---

## 8. Suggested C++ port order & test plan

1. Box magnitude core (`photometric_magnitudes` box path + `filter_passbands`
   box top-hat + `int_ref`) — pure numeric, unit-testable against §5.
2. `round_err`/`round2` — port exactly; test against §6 examples.
3. IUE path (highest value; the endpoints all work live): VizieR VI/110 query →
   filename build → INES download+gunzip → FITS read (WAVELENGTH/FLUX/SIGMA/
   QUALITY) → box mags + generic 0.04 err + `average_boxes(thres=0.32,
   mag_err_sys=0.03)`.
4. `average_boxes`.
5. Gaia XP (§2) — implement but expect the fixture to have **no** XP rows under
   current ISIS (see §0). If the C++ author *fixes* the URL to
   `data-server/data?`, XP boxes will appear — that would DIVERGE from the ISIS
   fixture, so gate it behind an off-by-default flag or regenerate the fixture
   with a patched ISIS.
6. MAST (§7) — lowest priority.

**Validate** against `validation/spectra/HZ44_1473687671071803520/photometry.dat`
(§9). The box rows a C++ impl must reproduce byte-for-byte are exactly the
`system==box` rows there (IUE per-spectrum rows, plus the `Average` rows, plus
MAST FUSE rows).

---

## 9. Reference fixture — `validation/spectra/HZ44_1473687671071803520/`

**Star:** HZ 44 = Gaia DR3 1473687671071803520 (iHe-sdOB, G=11.6). Chosen
because it has **both** `has_xp_sampled=true` (verified via Gaia TAP) and
abundant IUE data (SWP/LWR/LWP, low+high dispersion), plus FUSE & FOS in MAST.
Coordinates used: ra=200.8969289357, dec=36.1332067056.

**Generated by** (saved alongside as `gen_fixture.sl`):
```
cd validation/spectra/HZ44_1473687671071803520
~/Projects/ISIS_install/bin/isis-script gen_fixture.sl   # rc=0, ~run.log captured
```
`gen_fixture.sl` calls exactly the template's standard spectra query:
`query_photometry(200.8969289357, 36.1332067056; IUE, MAST, search_radius=4./60/60)`
then `query_reddening` + `photo.write("photometry.dat")`.

**Artifacts in the dir:** `photometry.dat` (the fixture), `gen_fixture.sl`
(driver), `run.log` (full ISIS stdout), `IUE/*.FITS.gz` (15 downloaded IUE
spectra), `MAST/*.FITS.gz` (8 downloaded FUSE+FOS spectra).

**Box rows actually produced (2026-07-06 live run):**
| source | rows | flag | notes |
|---|---|---|---|
| Gaia XP (`Gaia_DR3_XP`) | **0** | — | confirms §0 gotcha: AIP mirror answers first / esac URL 404s ⇒ `sgaia==NULL` |
| IUE (`*_VI/110/inescat`) | 21 | all **2** | 15 spectra; boxes {1300_1800,2000_2500,2500_3000}; all originally flag 0 (none tripped AbnCode/small-aper) ⇒ all set to 2 by `average_boxes` |
| IUE `Average` | 3 | **0** | one per IUE box; uncertainty 0.04/0.05; the rows the fitter actually uses |
| MAST FUSE (`FUSE-*_MAST`) | 9 | all **−1** | 3 FUSE spectra × boxes {950_1000,1000_1080,1091_1120}; forced −1 |
| MAST FOS | **0** | — | 5 FOS spectra downloaded+read but yielded **no rows** (quirk — see below) |

**Rows a C++ port must reproduce byte-for-byte:** every `system==box` row in
`photometry.dat` — the 21 IUE flag-2 rows, the 3 IUE flag-0 `Average` rows, and
the 9 MAST FUSE flag-−1 rows. Magnitudes/uncertainties/angu_dist strings must
match (e.g. `Average 1300_1800 10.057 0.04`, `Average 2500_3000 10.558 0.05`).

**FOS quirk (documented, unresolved):** the 5 HST/FOS spectra
(`HST-FOS-0.3_*`) were downloaded and read without error, but produced zero box
rows (no NaN/error message in `run.log`). Likely their dynamic
`[wmin+15:wmax-15:#4]` boxes evaluate to NaN in `photometric_magnitudes`, or the
FILETYPE/`spectralaxisname` handling skips them. Since MAST/FOS is the
lowest-priority path, a C++ port only needs to match "FOS ⇒ no rows" here;
investigate against these exact FITS files if FOS is later ported.

**Reproducibility caveat:** IUE/MAST catalog contents are stable, but if a
future run finds the esac datalink reachable AND the URL fixed, Gaia XP rows
would newly appear and diverge from this fixture — regenerate then.

---

## 10. IMPLEMENTATION STATUS (2026-07-06) — done & validated

Implemented in **`src/spectra.{hpp,cpp}`** (`add_spectra_boxes`), called from
`query_photometry` (`src/query.cpp`) after the catalogue loop and before the
final by-system stable sort — exactly where the S-Lang runs it. Hooks:

- `QueryOptions` gained `iue`/`mast` (default **on**, mirroring the template's
  standard `; IUE, MAST` call) and `spectra_dir` (parent of the `IUE/` and
  `MAST/` download caches; `main.cpp` passes the output file's directory —
  the ISIS equivalent of `$cwd`). Both the standalone `--query` path and the
  config-driven auto-query enable IUE+MAST.
- `vizier_fetch()` exported from `query.cpp` (single-catalogue VOTable query
  used for VI/110/inescat, same mirror fall-over as the bulk query).
- `fits_read_key_string()` added to `src/fitsio_util.{hpp,cpp}` (FILETYPE;
  missing key ⇒ "" like ISIS's NULL, open failure ⇒ throw ⇒ spectrum skipped).
- MAST SSAP via the exact STILTS command `ssapquery.sl` uses
  (`stilts cone servicetype=ssa serviceurl=... ofmt=csv`, popen + CSV parse).
- Downloads cached as gzipped FITS via zlib (`PkgConfig::ZLIB` in
  CMakeLists.txt); a file is fetched only if `<dir>/<name>.gz` is absent —
  cfitsio reads .gz natively.
- `average_boxes` replicates ISIS `moment` exactly (**Welford's online
  mean/variance**, `isis/src/math.c mean_stddev_doubles`; sample sdev N−1) —
  the naive two-pass mean can differ by 1 ulp and flip `round2` at a boundary.
  `nint` = round-half-away-from-zero (verified against the ISIS binary);
  `array_sort` is stable (verified) ⇒ `std::stable_sort`.

### Quirks confirmed & replicated (beyond §0-§9)

1. **`photometric_magnitudes` right-pad bug** (the one §5 didn't know about):
   at `photometric_magnitudes.sl:775` the right-padding condition is
   `ind[-1] < length(ind)-1` — the last *index* compared against the length of
   `ind` itself, which is never true. So the flux resample
   `interpol(u, l[ind], f[ind])` linearly **extrapolates at u == wmax** from
   the last two in-box points instead of interpolating across the edge
   (left edge is padded and interpolates fine). Observable: LWP07276LL
   2000_2500 = 10.319 (fixture) vs 10.326 with edge interpolation; also
   LWP07278LL (both boxes). The C++ `box_magnitude` replicates: inclusive
   `ind`, left-pad only, per-point `fu*fu_filter` (fu_filter is exactly
   `norm`), trapezoid, `-2.5*log10(int) + int_ref`.
   The **error** integral is a *different* recipe (driver-local): strict
   `ind`, `u = union(l_start, wave[ind], l_end)`, `interpol` over the **full**
   filtered arrays (interpolates at both edges). Both are implemented as
   specified; they intentionally disagree at the edges.
2. **FOS ⇒ 0 rows, root cause found**: the 5 HZ 44 HST/FOS spectra have
   `FILETYPE = 'WAV'` in the primary header (wavelength-calibration
   exposures) ⇒ skipped by the `filetype!="WAV"` guard *after* download.
   They are downloaded at all only because ssapquery's intended filters
   `string(query.object) != "WAVE"` / `... != "RFPT"` are inert S-Lang bugs:
   `string()` of an *array* yields `"String_Type[n]"`, so both comparisons are
   scalar-true no-ops. C++ replicates: filter FUSE/FOS + `der_snr>=2` only,
   download, then skip on FILETYPE=="WAV". (Also `string_match(fluxcal,
   "ABSOLUTE")!=1` means "does not match at position 1" ⇒ C++ tests
   starts-with.)
3. **Shared `boxes` variable**: one function-local in the S-Lang spans the IUE
   and MAST sections; `average_boxes` at the end of MAST sees whatever the
   last *processed* spectrum set (FUSE list here; would be the IUE list if no
   MAST spectrum passed the FILETYPE guard). Replicated with a shared vector.
4. **`angular_separation` 1-ulp**: the S-Lang converts rad→deg by
   `angusep *= 180./PI` (multiply by reciprocal); the C++ used
   `acos(temp)/d2r` (divide). 1-ulp difference, observable in HZ 44's 2MASS
   `angu_dist_arcsec` (…55 vs …53). Fixed in `query.cpp` (multiply form) —
   this also benefits ordinary catalogue rows.
5. **Gaia XP**: implemented faithfully to current ISIS (mirror probe order
   `aip.de` → `esac`; bail unless the *first reachable* mirror is esac; then
   the broken `<tap>data?...` datalink URL). Confirmed: zero `Gaia_DR3_XP`
   rows for HZ 44, matching the fixture. **Possible future opt-in
   divergence:** switching the download to the working
   `https://gea.esac.esa.int/data-server/data?ID=Gaia+DR3+<id>&RETRIEVAL_TYPE=XP_SAMPLED&format=FITS`
   endpoint would produce XP boxes ISIS doesn't currently emit — parity rules,
   so this is deliberately NOT done. The full XP post-processing (units,
   +30% error, u-band bumps) and `box_from_gaia` box math ARE implemented and
   would activate if ISIS's URL ever starts working.

### Validation — gate (a): HZ 44 fixture (byte-for-byte)

`./build/sedfit --query 200.8969289357 36.1332067056 query.dat` vs
`validation/spectra/HZ44_1473687671071803520/photometry.dat`:

| section | result |
|---|---|
| header (RA/DEC, meanSFD, meanSandF) | identical |
| 21 IUE flag-2 rows | **byte-identical, same order** |
| 3 IUE `Average` flag-0 rows | **byte-identical** |
| 9 MAST FUSE flag-−1 rows | **byte-identical** |
| Gaia XP rows | 0 = 0 ✓ |
| FOS rows | 0 = 0 ✓ (5 spectra downloaded+skipped, like ISIS) |
| 46 ordinary catalogue rows | **byte-identical** (after the angular_separation 1-ulp fix) |
| stderr messages | line-for-line match with the fixture run.log download/Found sequence |

Both cache-hit (pre-seeded `IUE/`+`MAST/`) and fresh-download runs validated
(INES + MAST endpoints live).

## Final status (2026-07-06, manager session)

Fresh-download path closed: `archive.stsci.edu/pub` (FOS/FUSE vospectra) fails
the default TLS 1.3 handshake (curl error 35; ISIS never saw it because wget
uses GnuTLS). Fix in `src/http.cpp::perform`: on `CURLE_SSL_CONNECT_ERROR`,
retry once capped at TLS 1.2. Re-validated live from an empty cache:
0 download failures, all 33 box rows + header byte-identical to the fixture.
Residual full-file diffs on HZ 44 are fixture-capture artifacts (ISIS's bulk
VizieR fetch missed II/349 PS1-DR1 + an APASS/SDSS sub-table at capture; C++
rows are correct and more complete — same class as the AllWISE/DELVE transients
documented in PHASE2.md). Single-star regression (1196375273586921472): all
returned flag-0 rows byte-identical; 4 missing DELVE rows were a transient
datalab outage during the check, not a code change (path untouched).

Spectra subsystem: COMPLETE. The query stream (Phase 2) is now fully ported:
catalogs, flags, reddening, Gaia astrometry, SIMBAD, IUE/MAST box spectra
(Gaia XP faithful-to-ISIS: zero rows while ISIS's datalink URL is dead).
