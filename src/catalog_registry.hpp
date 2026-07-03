// The photometric catalog registry (port of photometric_catalogs() +
// is_vizier_catalogue()). Rows come from the auto-generated catalog_table.inc.
#pragma once

#include <string>
#include <vector>

namespace sed {

struct CatalogRow {
  std::string catalogue;       // e.g. "II/335/galex_ais" or "PS1_DR2"
  std::string ra_colname;      // catalog's RA column
  std::string dec_colname;     // catalog's DEC column
  std::string system;          // photometric system ("GALEX", "Johnson", ...)
  std::string filter_colname;  // magnitude/color column in the catalog
  std::string error_colname;   // uncertainty column, or "NaN" if none
  std::string passband;        // our passband label ("FUV", "BmV", ...)
  std::string type;            // "magnitude" or "color"
  double dec_min = -90.0;
  double dec_max = 90.0;
  double angular_accuracy = 0.0;  // arcsec; sets a per-catalog min search radius
};

// VizieR naming convention: a catalogue name with >=2 '/' is a VizieR table.
bool is_vizier_catalogue(const std::string& cat_name);

// The full registry (all 308 rows), loaded once.
const std::vector<CatalogRow>& catalog_registry();

// Unique catalogue names in registry order (first-appearance order), matching
// the S-Lang `unique(...)[array_sort(...)]` set but preserving determinism.
std::vector<std::string> unique_catalogues();

}  // namespace sed
