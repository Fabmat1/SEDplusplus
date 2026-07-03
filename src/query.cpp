#include "query.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "catalog_registry.hpp"
#include "http.hpp"
#include "texval.hpp"
#include "util.hpp"
#include "votable.hpp"

namespace sed {

namespace {

constexpr double kPi = 3.14159265358979323846;

bool is_nan(double x) { return std::isnan(x); }

double flat_dist_arcsec(double ra, double dec, double r, double d) {
  double dr = (r - ra) * std::cos(dec * kPi / 180.0);
  double dd = d - dec;
  return 3600.0 * std::sqrt(dr * dr + dd * dd);
}

// A single catalogue row, exposed by (case-insensitive-fallback) field name.
// Empty cells read back as "NaN", mirroring parse_vizier_votable_multi.
struct Row {
  const votable::Table* t = nullptr;
  int row = -1;
  double angu_dist_arcsec = 0.0;

  bool has(const std::string& name) const {
    return t->field_index(name) >= 0 || t->field_index_ci(name) >= 0;
  }
  std::string s(const std::string& name) const {
    int i = t->field_index(name);
    if (i < 0) i = t->field_index_ci(name);
    if (i < 0) return "NaN";
    const std::string& v = t->rows[row][i];
    return v.empty() ? "NaN" : v;
  }
  double d(const std::string& name) const { return std::atof(s(name).c_str()); }

  // Read a stored quantity (magnitude/uncertainty) honouring the column's
  // VOTable datatype: a single-precision "float" column is re-quantised to its
  // shortest float decimal, so "%S" prints e.g. 15.0711 rather than the widened
  // 15.071100234985352 (matches the S-Lang, which holds these as Float_Type).
  double dstore(const std::string& name) const {
    int i = t->field_index(name);
    if (i < 0) i = t->field_index_ci(name);
    if (i < 0) return std::nan("");
    const std::string& v = t->rows[row][i];
    if (v.empty()) return std::nan("");
    if (t->datatype(i) == "float") {
      float f = std::strtof(v.c_str(), nullptr);
      char buf[32];
      auto r = std::to_chars(buf, buf + sizeof(buf), f);
      return std::strtod(std::string(buf, r.ptr).c_str(), nullptr);
    }
    return std::atof(v.c_str());
  }
};

struct Ctx {
  double ra, dec, l, b;
  double disk_b, disk_l;
};

// Locate the closest in-radius row of a table, optionally pre-filtered.
// radius_arcsec is this catalogue's own radius; candidates pass the flat
// distance cut (as VizieR's parser does), the winner minimises the true
// angular separation (as closest_object does).
std::optional<Row> closest_row(const votable::Table* t,
                               const std::string& ra_col,
                               const std::string& dec_col, const Ctx& c,
                               double radius_arcsec,
                               const std::function<bool(const Row&)>& keep) {
  if (!t || t->rows.empty()) return std::nullopt;
  int ira = t->field_index(ra_col);
  if (ira < 0)
    for (const char* alt : {"_RA", "_RA.icrs", "RAJ2000", "RA_ICRS", "ra"}) {
      ira = t->field_index(alt);
      if (ira >= 0) break;
    }
  int ide = t->field_index(dec_col);
  if (ide < 0)
    for (const char* alt : {"_DE", "_DE.icrs", "DEJ2000", "DE_ICRS", "dec"}) {
      ide = t->field_index(alt);
      if (ide >= 0) break;
    }
  if (ira < 0 || ide < 0) return std::nullopt;

  int best = -1;
  double best_ang = 1e18;
  for (size_t r = 0; r < t->rows.size(); ++r) {
    double rr = std::atof(t->rows[r][ira].c_str());
    double dd = std::atof(t->rows[r][ide].c_str());
    if (flat_dist_arcsec(c.ra, c.dec, rr, dd) > radius_arcsec) continue;
    Row cand{t, (int)r, 0.0};
    if (keep && !keep(cand)) continue;
    double ang = angular_separation_deg(c.ra, c.dec, rr, dd);
    if (ang < best_ang) { best_ang = ang; best = (int)r; }
  }
  if (best < 0) return std::nullopt;
  Row out{t, best, 3600.0 * best_ang};
  return out;
}

// ---- Redundancy helpers over the entries assembled so far -------------------

// Is there already an entry from `cat` for this passband whose flag is > -100?
bool has_prior(const std::vector<PhotEntry>& e, const std::string& cat,
               const std::string& band) {
  for (const auto& x : e)
    if (x.vizier_catalog == cat && x.passband == band && x.flag > -100)
      return true;
  return false;
}

// Append an entry, applying the flag==0 duplicate resolution (lines 1820-1838):
// if an equal (system,passband,magnitude,flag==0) entry exists, keep whichever
// has the smaller uncertainty at flag 0 and flag the other as redundant (1).
void add_entry(std::vector<PhotEntry>& e, int flag, const std::string& system,
               const std::string& band, double mag, double err,
               const std::string& type, double dist, const std::string& cat) {
  if (flag == 0) {
    for (auto& x : e) {
      if (x.system == system && x.passband == band && x.magnitude == mag &&
          x.flag == 0) {
        if (!is_nan(err) && (is_nan(x.uncertainty) || x.uncertainty > err))
          x.flag = 1;  // old entry worse -> demote it, keep new at 0
        else
          flag = 1;  // new entry worse -> demote it
        break;
      }
    }
  }
  e.push_back({flag, system, band, mag, err, type, dist, cat});
}

// GALEX GR6+7 bright-magnitude correction (Wall et al. 2019, Eq. 5).
void galex_correct(const std::string& band, double& mag, double& err) {
  if (band == "NUV" && mag < 16.95) {
    double corr = 3.778 + std::sqrt(24.337 * mag - 241.018);
    auto re = round_err(std::sqrt(err * err + (corr - mag) * (corr - mag)));
    mag = round2(corr, (int)re.digit - 1);
    err = re.value;
  } else if (band == "FUV" && mag < 15.95) {
    double corr = 6.412 + std::sqrt(17.63 * mag - 192.135);
    auto re = round_err(std::sqrt(err * err + (corr - mag) * (corr - mag)));
    mag = round2(corr, (int)re.digit - 1);
    err = re.value;
  }
}

// Apply the per-catalogue quality/flag logic for one filter row and append it.
// Returns nothing; mutates `entries`.
void process_filter(const std::string& cat, const CatalogRow& reg,
                    const Row& row, const Ctx& c,
                    std::vector<PhotEntry>& entries) {
  double mag_error =
      (reg.error_colname != "NaN" && row.has(reg.error_colname))
          ? row.dstore(reg.error_colname)
          : std::nan("");
  double magnitude =
      row.has(reg.filter_colname) ? row.dstore(reg.filter_colname)
                                  : std::nan("");
  if (is_nan(magnitude)) return;  // S-Lang skips NaN-magnitude entries
  int flag = 0;
  const std::string& band = reg.passband;
  bool disk = std::abs(c.b) < c.disk_b && std::abs(c.l - 180.0) > c.disk_l;

  if (cat == "II/246/out") {  // 2MASS
    auto pad3 = [](std::string s) {
      while (s.size() < 3) s = "0" + s;
      return s;
    };
    std::string Rflg = pad3(row.s("Rflg"));
    std::string Cflg = pad3(row.s("Cflg"));
    int pos = band == "J" ? 0 : (band == "H" ? 1 : 2);
    int rf = Rflg[pos] - '0';
    std::string cf(1, Cflg[pos]);
    double Xflg = row.d("Xflg");
    if (cf != "0" || rf < 1 || (rf > 3 && rf != 6) || Xflg > 0 ||
        row.d(band + "snr") < 2)
      flag = -1;
    if (is_nan(mag_error) || is_nan(magnitude)) flag = -1;
  } else if (cat == "I/350/gaiaedr3") {  // Gaia EDR3
    double G = row.d("Gmag");
    double bp_rp = row.d("E(BP/RP)Corr");
    if (is_nan(bp_rp) ||
        std::abs(bp_rp) >
            2 * (0.0059898 + 8.817481e-12 * std::pow(G, 7.618399)))
      flag = -1;
    if (band == "G" && G >= 2.0 && G < 8)
      magnitude += -0.09892 + 0.059 * G - 0.009775 * G * G +
                   0.0004934 * G * G * G;
    if (band == "GBP" && G >= 2.0 && G < 3.94)
      magnitude += -0.9921 - 0.02598 * G + 0.1833 * G * G - 0.02862 * G * G * G;
    if (band == "GRP" && magnitude >= 2.0 && magnitude < 3.45)
      magnitude += -14.94 + 14.41 * magnitude - 4.657 * magnitude * magnitude +
                   0.503 * magnitude * magnitude * magnitude;
    magnitude = round2(magnitude, (int)round_err(mag_error).digit - 1);
  } else if (cat == "II/335/galex_ais") {  // GALEX
    galex_correct(band, magnitude, mag_error);
    if (is_nan(magnitude) || (magnitude < 10 && mag_error > 1.8)) flag = -1;
    // artefact flag: bits at k=2,3,7 -> values 2^1|2^2|2^6 = 70 (window/
    // dichroic reflections), matching the high-to-low decomposition.
    long afl = (long)row.d(band == "FUV" ? "Fafl" : "Nafl");
    if ((afl & 70L) != 0) flag = -1;
  } else if (cat == "II/363/unwise") {  // unWISE (flux -> mag)
    double snr = magnitude / mag_error;
    mag_error = std::sqrt(mag_error * mag_error +
                          std::pow(row.d("e_F" + band + "lbs"), 2));
    auto re = round_err(std::abs(-2.5 / std::log(10.0) * mag_error / magnitude));
    magnitude = round2(22.5 - 2.5 * std::log10(magnitude), (int)re.digit - 1);
    mag_error = re.value;
    if (magnitude > 16.5 && mag_error > 0.)
      mag_error = std::sqrt(mag_error * mag_error +
                            std::pow(0.0075 * magnitude, 2));
    if (disk) flag = -4;
    if (is_nan(magnitude) || row.d("nm" + band) == 0 || snr < 2 ||
        row.d("q_" + band) < 0.5)
      flag = -1;
    else if (((long)row.d("Flags" + band) & 32L) != 0)  // bit 6 (Saturated)
      flag = -1;
  } else if (cat == "II/365/catwise") {  // CatWISE2020
    if (magnitude > 16.5 && mag_error > 0.)
      mag_error = std::sqrt(mag_error * mag_error +
                            std::pow(0.0075 * magnitude, 2));
    if (disk) flag = -4;
    std::string kcol = band;
    kcol[0] = 'k';  // W1 -> k1
    if (row.d("snr" + band + "pm") < 3 || row.d("n" + band) == 0 ||
        row.s(kcol) == "NaN" || row.d(kcol) != 3)
      flag = -1;
    if (flag == 0 && has_prior(entries, "II/363/unwise", band)) flag = -3;
  } else if (cat == "II/328/allwise") {  // AllWISE
    if (magnitude > 16.5 && mag_error > 0.)
      mag_error =
          std::sqrt(mag_error * mag_error + std::pow(0.01 * magnitude, 2));
    if (disk) flag = -4;
    std::string snrcol = "snr" + band.substr(1);  // W1 -> snr1
    if (row.d("ex") > 2 || row.d(snrcol) < 2) flag = -1;
    if (flag == 0 && !entries.empty()) flag = -3;  // always demote AllWISE
  } else if (cat == "II/336/apass9") {  // APASS DR9
    flag = -5;
    mag_error = std::sqrt(mag_error * mag_error + 0.04 * 0.04);
    if (row.d("u_" + reg.error_colname) == 1) flag = -1;
    std::string efield = "g_" + reg.error_colname;
    if (row.has(efield) && row.d(efield) == 1) flag = -1;
  } else if (cat == "V/147/sdss12") {  // SDSS DR12
    if (row.s("q_mode") != "+" || row.d("Q") == 1) flag = -1;
  } else if (cat == "II/349/ps1") {  // Pan-STARRS DR1 (redundant with DR2)
    flag = 1;
    if (is_nan(mag_error))
      flag = -1;
    else {
      // Qual bits 6 (64) or 7 (128): suspect / poor-quality stack object.
      if (((long)row.d("Qual") & 192L) != 0) flag = -1;
      if (row.d("e_" + band + "mag") > 0.12 || (magnitude > 18.2 && band == "y"))
        flag = -1;
    }
  } else if (cat == "V/165/igapsdr1") {  // IGAPS
    if (row.d("errBits") > 0.99) flag = -1;
  } else if (cat == "PS1_DR2") {  // Pan-STARRS DR2 (MAST TAP)
    double err_c = row.d(band + "meanpsfmagerr");
    double mag_c = row.d(band + "meanpsfmag");
    if (mag_c < -998 || err_c < -998 || row.d(band + "flags") <= 0 ||
        err_c > 0.12)
      flag = -1;
    if (mag_c > 18.2 && band == "y") flag = -1;
    if (mag_c < 10.0) flag = -1;
  } else if (cat == "DELVE_DR2") {  // DECam Local Volume Exploration (datalab)
    if (row.d("wavg_mag_psf_" + band) > 98.9 ||
        row.d("wavg_magerr_psf_" + band) > 98.9 ||
        row.d("wavg_magrms_psf_" + band) > 0.1 ||
        row.d("class_star_" + band) < 0.5 || row.d("flags_" + band) > 3)
      flag = -1;
    // (DES_DR2 not queried -> no -3 supersession)
  } else if (cat == "GPS_DR11" || cat == "LAS_DR11" || cat == "GCS_DR11" ||
             cat == "DXS_DR11" || cat == "UDS_DR11") {  // UKIDSS (WSA)
    std::string errbits = reg.filter_colname;
    auto p = errbits.find("AperMag3");
    if (p != std::string::npos) errbits.replace(p, 8, "ErrBits");
    if (!(magnitude > -9) || !(mag_error > 0) || row.d(errbits) > 12 ||
        row.d("pStar") < 0.05)
      flag = -1;
  }
  // (other catalogues fall through with flag 0 and direct magnitude/error)

  add_entry(entries, flag, reg.system, band, magnitude, mag_error, reg.type,
            row.angu_dist_arcsec, cat);
}

// Issue a synchronous TAP query and return its first non-empty table.
// Mirrors tapquery.sl: POST REQUEST/LANG/FORMAT/QUERY to <url>/sync. We request
// VOTable (MAST ignores FORMAT and returns VOTable anyway) and reuse the
// VOTable parser instead of the S-Lang's CSV path.
std::optional<votable::Table> tap_query(const std::string& base_url,
                                        const std::string& adql) {
  std::string url = base_url;
  if (url.size() < 5 || url.compare(url.size() - 5, 5, "/sync") != 0) {
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/sync";
  }
  std::string body = "REQUEST=doQuery&LANG=ADQL&FORMAT=votable&QUERY=" +
                     sed::http::url_encode(adql);
  http::Options opt;
  opt.timeout_s = 180;
  sed::http::Response resp;
  try {
    resp = sed::http::post(url, body, "application/x-www-form-urlencoded", opt);
  } catch (const std::exception&) {
    return std::nullopt;
  }
  if (!resp.ok()) return std::nullopt;
  auto doc = votable::parse(resp.body);
  for (auto& t : doc.tables)
    if (!t.rows.empty()) return t;
  return std::nullopt;
}

// Build the ADQL box query for a non-VizieR catalogue, or "" if unsupported.
std::string tap_url_for(const std::string& cat) {
  if (cat == "PS1_DR2") return "https://mast.stsci.edu/vo-tap/api/v0.1/ps1dr2/";
  if (cat == "DELVE_DR2") return "https://datalab.noirlab.edu/tap";
  if (cat == "GPS_DR11" || cat == "LAS_DR11" || cat == "GCS_DR11" ||
      cat == "DXS_DR11" || cat == "UDS_DR11")
    return "http://tap.roe.ac.uk/wsa";
  return "";
}

std::string tap_adql_for(const std::string& cat, const CatalogRow& r,
                         double ra, double dec, double rad_deg) {
  double lo_ra = ra - rad_deg, hi_ra = ra + rad_deg;
  double lo_de = dec - rad_deg, hi_de = dec + rad_deg;
  char q[512];
  const char *table = nullptr, *rac = "RA", *dec_c = "DEC";
  if (cat == "PS1_DR2") {
    table = "dbo.MeanObjectView";
    rac = r.ra_colname.c_str();  // ramean
    dec_c = r.dec_colname.c_str();
  } else if (cat == "DELVE_DR2") {
    table = "delve_dr2.objects";
  } else if (cat == "GPS_DR11") {
    table = "UKIDSSDR11PLUS.gpsSource"; rac = "ra"; dec_c = "dec";
  } else if (cat == "LAS_DR11") {
    table = "UKIDSSDR11PLUS.lasSource"; rac = "ra"; dec_c = "dec";
  } else if (cat == "GCS_DR11") {
    table = "UKIDSSDR11PLUS.gcsSource"; rac = "ra"; dec_c = "dec";
  } else if (cat == "DXS_DR11") {
    table = "UKIDSSDR11PLUS.dxsSource"; rac = "ra"; dec_c = "dec";
  } else if (cat == "UDS_DR11") {
    table = "UKIDSSDR11PLUS.udsSource"; rac = "ra"; dec_c = "dec";
  } else {
    return "";
  }
  std::snprintf(q, sizeof(q),
                "SELECT * FROM %s WHERE (%s BETWEEN %f AND %f) AND (%s BETWEEN "
                "%f AND %f)",
                table, rac, lo_ra, hi_ra, dec_c, lo_de, hi_de);
  return q;
}

}  // namespace

void radec2lb(double ra, double dec, double& l, double& b) {
  const double d2r = kPi / 180.0;
  double at = ra * d2r, dt = dec * d2r;
  double a0 = 192.8594812065348 * d2r;
  double d0 = 27.12825118085622 * d2r;
  double l0 = 122.9319185680026 * d2r;
  double gal_b = std::asin(std::sin(dt) * std::sin(d0) +
                           std::cos(dt) * std::cos(d0) * std::cos(at - a0));
  double gal_l =
      l0 - std::atan2(std::cos(dt) * std::sin(at - a0),
                      std::sin(dt) * std::cos(d0) -
                          std::cos(dt) * std::sin(d0) * std::cos(at - a0));
  if (gal_l < 0.0) gal_l += 360.0 * d2r;
  l = gal_l / d2r;
  b = gal_b / d2r;
}

double angular_separation_deg(double ra1, double dec1, double ra2,
                              double dec2) {
  if (ra1 == ra2 && dec1 == dec2) return 0.0;
  const double d2r = kPi / 180.0;
  ra1 *= d2r; dec1 *= d2r; ra2 *= d2r; dec2 *= d2r;
  double temp = std::sin(dec1) * std::sin(dec2) +
                std::cos(dec1) * std::cos(dec2) * std::cos(ra1 - ra2);
  if (std::abs(temp) > 1)
    temp -= (temp > 0 ? 1 : -1) * 2.220446049250313e-16;
  return std::acos(temp) / d2r;
}

PhotometryTable query_photometry(double ra, double dec,
                                 const QueryOptions& opt) {
  PhotometryTable table;
  table.ra = ra;
  table.dec = dec;
  Ctx ctx{ra, dec, 0, 0, opt.disk_b, opt.disk_l};
  radec2lb(ra, dec, ctx.l, ctx.b);

  const auto& reg = catalog_registry();

  // Unique catalogues in byte-sorted processing order (array_sort of unique).
  auto uniq = unique_catalogues();
  std::sort(uniq.begin(), uniq.end());

  // Registry rows grouped by catalogue (registry appearance order preserved).
  std::map<std::string, std::vector<CatalogRow>> by_cat;
  for (const auto& r : reg) by_cat[r.catalogue].push_back(r);

  // radius (arcsec) for a catalogue: max(search_radius, angular_accuracy/3600).
  auto radius_arcsec = [&](const CatalogRow& r) {
    double deg = opt.force_search_radius
                     ? opt.search_radius_deg
                     : std::max(opt.search_radius_deg, r.angular_accuracy / 3600.0);
    return deg * 3600.0;
  };

  // ---- One bulk VizieR request for all in-dec VizieR catalogues -------------
  std::vector<std::string> vz_cats;
  double vz_maxr = 0;
  for (const auto& cat : uniq) {
    const auto& rows = by_cat[cat];
    const auto& r0 = rows.front();
    if (!is_vizier_catalogue(cat)) continue;
    if (dec < r0.dec_min || dec > r0.dec_max) continue;
    vz_cats.push_back(cat);
    vz_maxr = std::max(vz_maxr, radius_arcsec(r0));
  }
  votable::Document doc;
  if (!vz_cats.empty()) {
    std::string source;
    for (size_t i = 0; i < vz_cats.size(); ++i)
      source += (i ? "," : "") + vz_cats[i];
    char params[8192];
    std::snprintf(params, sizeof(params),
                  "-source=%s&-c=%.6f%%20%+.6f&-c.rs=%.4f&-out.add=_r&-out.all",
                  source.c_str(), ra, dec, vz_maxr);
    std::string url =
        std::string("http://vizier.u-strasbg.fr/viz-bin/votable?") + params;
    auto resp = sed::http::get(url);
    if (resp.ok()) doc = votable::parse(resp.body);
  }

  // ---- Process every catalogue in byte-sorted order -------------------------
  for (const auto& cat : uniq) {
    const auto& rows = by_cat[cat];
    const auto& r0 = rows.front();
    if (dec < r0.dec_min || dec > r0.dec_max) continue;

    std::optional<Row> row;
    std::optional<votable::Table> tap_table;  // must outlive process_filter
    if (is_vizier_catalogue(cat)) {
      const votable::Table* t = doc.find(cat);
      if (!t) continue;
      std::function<bool(const Row&)> keep;
      if (cat == "V/147/sdss12")  // primary detections only (mode==1)
        keep = [](const Row& rw) { return rw.d("mode") == 1; };
      row = closest_row(t, r0.ra_colname, r0.dec_colname, ctx,
                        radius_arcsec(r0), keep);
    } else if (!opt.skip_tap) {
      std::string url = tap_url_for(cat);
      if (url.empty()) continue;  // TAP catalogue not yet supported
      std::string adql =
          tap_adql_for(cat, r0, ra, dec, radius_arcsec(r0) / 3600.0);
      tap_table = tap_query(url, adql);
      if (!tap_table) continue;
      row = closest_row(&*tap_table, r0.ra_colname, r0.dec_colname, ctx,
                        radius_arcsec(r0), nullptr);
    } else {
      continue;
    }
    if (!row) continue;

    for (const auto& reg_row : rows)
      process_filter(cat, reg_row, *row, ctx, table.entries);
  }

  // Cosmetic: stable sort by photometric system (fit is order-independent).
  std::stable_sort(table.entries.begin(), table.entries.end(),
                   [](const PhotEntry& a, const PhotEntry& b) {
                     return a.system < b.system;
                   });
  return table;
}

}  // namespace sed
