// SIMBAD/Sesame name resolution (port of
// src/miscellaneous/resolve_simbad.sl as used by templates/photometry.sl:
// coordinates are resolved from an object identifier only when neither a
// cached results header nor explicit coordinates are available).
#pragma once

#include <optional>
#include <string>

namespace sed {

struct SimbadResult {
  long long oid = 0;
  double ra = 0, dec = 0;  // decimal degrees
  std::string main_id;
  std::string coordreference;
  int nbreferences = 0;
  double parallax = 0;        // mas, NaN if absent (not necessarily Gaia)
  double radialvelocity = 0;  // km/s, NaN if absent
};

// Normalize a user-supplied object name the way templates/photometry.sl does
// for its `star` before resolution: trim, a bare integer becomes a
// "Gaia DR3 <id>" designation (photometry.sl:13), and underscores become
// spaces (photometry.sl:15).
std::string normalize_star_name(std::string name);

// Resolve an object identifier to coordinates via the SIMBAD sim-tap ADQL
// service (port of resolve_simbad.sl). Returns std::nullopt when the name
// resolves to no row (the S-Lang returns NULL) or every service is
// unreachable. The caller decides how to treat a miss.
std::optional<SimbadResult> resolve_simbad(const std::string& name);

}  // namespace sed
