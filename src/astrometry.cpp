#include "astrometry.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

#include "http.hpp"
#include "query.hpp"  // angular_separation_deg
#include "votable.hpp"

namespace sed {

namespace {

const char* kUrlCorr = "https://gaia.ari.uni-heidelberg.de/tap/";
const char* kUrlsFallback[] = {"https://gaia.aip.de/tap",
                               "https://gea.esac.esa.int/tap-server/tap/"};

// Columns exactly as the S-Lang assembles them for the photometry.sl call:
// the caller's list, +ra,dec (closest_object), +phot_g_mean_mag
// (parallax_error correction), +parallax_corr (Heidelberg only).
const char* kColsBase =
    "designation,parallax,parallax_error,ruwe,ipd_gof_harmonic_amplitude,"
    "ipd_frac_multi_peak,visibility_periods_used,ra,dec,phot_g_mean_mag";

// Synchronous TAP query returning the first non-empty table. Same protocol as
// query.cpp's tap_query but requests FORMAT=votable/td: the Heidelberg DaCHS
// service defaults to base64 BINARY serialization, which our TABLEDATA-only
// VOTable parser cannot read.
std::optional<votable::Table> tap_query_td(const std::string& base_url,
                                           const std::string& adql) {
  std::string url = base_url;
  if (url.size() < 5 || url.compare(url.size() - 5, 5, "/sync") != 0) {
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/sync";
  }
  std::string body = "REQUEST=doQuery&LANG=ADQL&FORMAT=votable%2Ftd&QUERY=" +
                     http::url_encode(adql);
  http::Options opt;
  opt.timeout_s = 180;
  http::Response resp;
  try {
    resp = http::post(url, body, "application/x-www-form-urlencoded", opt);
  } catch (const std::exception&) {
    return std::nullopt;
  }
  if (!resp.ok()) return std::nullopt;
  auto doc = votable::parse(resp.body);
  for (auto& t : doc.tables)
    if (!t.rows.empty()) return t;
  return std::nullopt;
}

// Cell access with empty -> NaN, like the query layer.
double cell_num(const votable::Table& t, int row, const char* name) {
  int i = t.field_index(name);
  if (i < 0) i = t.field_index_ci(name);
  if (i < 0) return std::numeric_limits<double>::quiet_NaN();
  const std::string& v = t.rows[row][i];
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  return std::atof(v.c_str());
}

bool has_field(const votable::Table& t, const char* name) {
  return t.field_index(name) >= 0 || t.field_index_ci(name) >= 0;
}

}  // namespace

std::optional<AstrometryResult> query_astrometry(double ra, double dec,
                                                 double search_radius_deg) {
  // URL / JOIN selection (query_astrometry.sl:140-153): with corrected_values
  // and a DR3 catalog, probe the Heidelberg corrections service; if it answers
  // it becomes the only URL and the query gains the corrections JOIN +
  // parallax_corr column. Otherwise fall back to the default URL list with a
  // plain gaia_source query.
  std::vector<std::string> urls;
  std::string cols = kColsBase;
  std::string join;
  if (http::probe(kUrlCorr, 4, 6)) {
    urls = {kUrlCorr};
    join = " INNER JOIN gaiadr3.gaia_source_corrections USING (source_id)";
    cols += ",parallax_corr";
  } else {
    urls.assign(std::begin(kUrlsFallback), std::end(kUrlsFallback));
  }

  char adql[1024];
  std::snprintf(adql, sizeof adql,
                "SELECT %s FROM gaiadr3.gaia_source%s WHERE 1=CONTAINS("
                "POINT('ICRS', ra, dec), CIRCLE('ICRS', %f, %f, %f ))",
                cols.c_str(), join.c_str(), ra, dec, search_radius_deg);

  std::optional<votable::Table> table;
  for (const auto& u : urls) {
    if (&u != &urls.front() || urls.size() > 1) {
      if (!http::probe(u, 4, 6)) {
        std::fprintf(stderr, "query_astrometry: URL offline: %s\n", u.c_str());
        continue;
      }
    }
    table = tap_query_td(u, adql);
    if (table) break;
    std::fprintf(stderr, "query_astrometry: query not successful for %s\n",
                 u.c_str());
  }
  if (!table) return std::nullopt;
  const votable::Table& t = *table;

  // Closest row by true angular separation (closest_object).
  int best = -1;
  double best_ang = 1e18;
  for (int r = 0; r < (int)t.rows.size(); ++r) {
    double rr = cell_num(t, r, "ra"), dd = cell_num(t, r, "dec");
    double ang = angular_separation_deg(ra, dec, rr, dd);
    if (ang < best_ang) {
      best_ang = ang;
      best = r;
    }
  }
  if (best < 0) return std::nullopt;

  AstrometryResult out;
  {
    int i = t.field_index("designation");
    if (i < 0) i = t.field_index_ci("designation");
    if (i >= 0) out.designation = t.rows[best][i];
  }
  out.ra = cell_num(t, best, "ra");
  out.dec = cell_num(t, best, "dec");
  double parallax = cell_num(t, best, "parallax");
  double parallax_error = cell_num(t, best, "parallax_error");
  double ruwe = cell_num(t, best, "ruwe");
  double ipd_gof = cell_num(t, best, "ipd_gof_harmonic_amplitude");
  double ipd_frac = cell_num(t, best, "ipd_frac_multi_peak");
  double vpu = cell_num(t, best, "visibility_periods_used");
  double G = cell_num(t, best, "phot_g_mean_mag");
  if (!std::isnan(ruwe)) out.ruwe = ruwe;
  if (!std::isnan(ipd_gof)) out.ipd_gof_harmonic_amplitude = ipd_gof;
  if (!std::isnan(vpu)) out.visibility_periods_used = vpu;

  // quality_warnings (query_astrometry.sl:195-205): DR3 thresholds.
  struct QI {
    const char* name;
    double value, threshold;
  } qis[] = {{"ruwe", ruwe, 1.4},
             {"ipd_gof_harmonic_amplitude", ipd_gof, 0.1},
             {"ipd_frac_multi_peak", ipd_frac, 10},
             {"visibility_periods_used", vpu, 1e10}};
  for (const auto& q : qis)
    if (q.value > q.threshold)
      std::fprintf(stderr,
                   "Warning in 'query_astrometry': The quality indicator '%s' "
                   "(%g) exceeds its recommended threshold (%g).\n",
                   q.name, q.value, q.threshold);

  // Parallax error inflation, El-Badry, Rix & Heintz 2021 Eq. (16)
  // (query_astrometry.sl:210-215).
  double parallax_error_corr = std::numeric_limits<double>::quiet_NaN();
  if (has_field(t, "parallax_error") && !std::isnan(G))
    parallax_error_corr =
        parallax_error * (0.21 * std::exp(-std::pow((G - 12.65) / 0.9, 2)) +
                          1.141 + 0.0040 * G - 0.00062 * G * G);
  double parallax_corr = has_field(t, "parallax_corr")
                             ? cell_num(t, best, "parallax_corr")
                             : std::numeric_limits<double>::quiet_NaN();

  // NaN-masking + corrected-value selection (photometry.sl:216, 229-248).
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  out.good = !std::isnan(parallax) && parallax > 0;
  if (out.good) {
    if (!(parallax > 0)) parallax = NaN;
    if (!(parallax / parallax_error >= 1.0)) parallax = NaN;
    if (has_field(t, "parallax_corr") && !std::isnan(parallax_error_corr)) {
      if (!(parallax_corr > 0)) parallax_corr = NaN;
      if (!(parallax_corr / parallax_error_corr >= 1.0)) parallax_corr = NaN;
      if (!std::isnan(parallax_corr) && parallax_corr > 0 &&
          !std::isnan(parallax_error_corr) && parallax_error_corr > 0) {
        // Use the Lindegren-2020 zero-point-corrected parallax together with
        // the inflated error (they are only ever substituted as a pair).
        parallax = parallax_corr;
        parallax_error = parallax_error_corr;
      }
    }
    out.good = !std::isnan(parallax) && parallax > 0;
  }
  out.parallax = parallax;
  out.parallax_error = parallax_error;
  return out;
}

}  // namespace sed
