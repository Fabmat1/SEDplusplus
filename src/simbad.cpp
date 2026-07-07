#include "simbad.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "http.hpp"
#include "votable.hpp"

namespace sed {

namespace {

// resolve_simbad.sl:80,83 — primary CDS service, Harvard mirror on offline.
const char* kUrlPrimary = "http://simbad.u-strasbg.fr/simbad/sim-tap/sync";
const char* kUrlFallback = "http://simbad.cfa.harvard.edu/simbad/sim-tap/sync";

double cell_num(const votable::Table& t, int row, const char* name) {
  int i = t.field_index(name);
  if (i < 0) i = t.field_index_ci(name);
  if (i < 0 || row < 0 || row >= (int)t.rows.size()) return std::nan("");
  const std::string& v = t.rows[row][i];
  if (v.empty()) return std::nan("");
  return std::atof(v.c_str());
}

std::string cell_str(const votable::Table& t, int row, const char* name) {
  int i = t.field_index(name);
  if (i < 0) i = t.field_index_ci(name);
  if (i < 0 || row < 0 || row >= (int)t.rows.size()) return "";
  return t.rows[row][i];
}

}  // namespace

std::string normalize_star_name(std::string name) {
  // strtrim
  size_t a = name.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = name.find_last_not_of(" \t\r\n");
  name = name.substr(a, b - a + 1);

  // Bare integer (e.g. a Gaia source_id) -> "Gaia DR3 <id>" (photometry.sl:13).
  bool all_digits = !name.empty();
  for (char c : name)
    if (!std::isdigit((unsigned char)c)) {
      all_digits = false;
      break;
    }
  if (all_digits) return "Gaia DR3 " + name;

  // Underscores -> spaces (photometry.sl:15).
  for (char& c : name)
    if (c == '_') c = ' ';
  return name;
}

std::optional<SimbadResult> resolve_simbad(const std::string& name) {
  // Escape single quotes for the ADQL string literal (double them).
  std::string id;
  for (char c : name) {
    id += c;
    if (c == '\'') id += '\'';
  }

  // ADQL exactly as resolve_simbad.sl:93-102 assembles it.
  std::string adql =
      "SELECT basic.OID, RA, DEC, main_id, "
      "coo_bibcode AS \"CoordReference\", nbref AS \"NbReferences\", "
      "plx_value AS \"Parallax\", rvz_radvel AS \"RadialVelocity\" "
      "FROM basic JOIN ident ON oidref = oid WHERE id = '" +
      id + "';";

  // Probe-select the mirror (resolve_simbad.sl:81-83): only one URL per call.
  const char* url =
      http::probe(kUrlPrimary, 4, 6) ? kUrlPrimary : kUrlFallback;

  // FORMAT=votable/td: SIMBAD's sync endpoint defaults to base64 BINARY
  // serialization for FORMAT=votable, which our TABLEDATA-only VOTable parser
  // cannot read (same reason astrometry.cpp requests votable/td).
  std::string body =
      "REQUEST=doQuery&PHASE=RUN&LANG=ADQL&FORMAT=votable%2Ftd&query=" +
      http::url_encode(adql);
  http::Options opt;
  opt.timeout_s = 120;
  http::Response resp;
  try {
    resp = http::post(url, body, "application/x-www-form-urlencoded", opt);
  } catch (const std::exception&) {
    return std::nullopt;
  }
  if (!resp.ok()) return std::nullopt;

  auto doc = votable::parse(resp.body);
  const votable::Table* table = nullptr;
  for (auto& t : doc.tables)
    if (!t.rows.empty()) {
      table = &t;
      break;
    }
  // No data row (0 lines / unresolvable name) -> NULL, like the S-Lang.
  if (!table) return std::nullopt;
  const votable::Table& t = *table;

  SimbadResult r;
  r.ra = cell_num(t, 0, "ra");
  r.dec = cell_num(t, 0, "dec");
  r.main_id = cell_str(t, 0, "main_id");
  r.coordreference = cell_str(t, 0, "CoordReference");
  {
    std::string s = cell_str(t, 0, "OID");
    if (s.empty()) s = cell_str(t, 0, "oid");
    r.oid = s.empty() ? 0 : std::atoll(s.c_str());
  }
  {
    std::string s = cell_str(t, 0, "NbReferences");
    r.nbreferences = s.empty() ? 0 : std::atoi(s.c_str());
  }
  r.parallax = cell_num(t, 0, "Parallax");
  r.radialvelocity = cell_num(t, 0, "RadialVelocity");

  // A resolved object without usable coordinates is treated as a miss.
  if (std::isnan(r.ra) || std::isnan(r.dec)) return std::nullopt;
  return r;
}

}  // namespace sed
