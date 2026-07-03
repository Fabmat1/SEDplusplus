#include "catalog_registry.hpp"

#include <stdexcept>

#include "util.hpp"

namespace sed {

#include "catalog_table.inc"

bool is_vizier_catalogue(const std::string& cat_name) {
  int slashes = 0;
  for (char c : cat_name)
    if (c == '/') ++slashes;
  return slashes >= 2;
}

const std::vector<CatalogRow>& catalog_registry() {
  static const std::vector<CatalogRow> rows = [] {
    std::vector<CatalogRow> out;
    const size_t n = sizeof(CATALOG_TABLE) / sizeof(CATALOG_TABLE[0]);
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      auto tok = split_ws(CATALOG_TABLE[i]);
      if (tok.size() != 11)
        throw std::runtime_error("bad catalog_table.inc row: " +
                                 std::string(CATALOG_TABLE[i]));
      CatalogRow r;
      r.catalogue = tok[0];
      r.ra_colname = tok[1];
      r.dec_colname = tok[2];
      r.system = tok[3];
      r.filter_colname = tok[4];
      r.error_colname = tok[5];
      r.passband = tok[6];
      r.type = tok[7];
      r.dec_min = std::stod(tok[8]);
      r.dec_max = std::stod(tok[9]);
      r.angular_accuracy = std::stod(tok[10]);
      out.push_back(std::move(r));
    }
    return out;
  }();
  return rows;
}

std::vector<std::string> unique_catalogues() {
  std::vector<std::string> out;
  for (const auto& r : catalog_registry()) {
    bool seen = false;
    for (const auto& c : out)
      if (c == r.catalogue) { seen = true; break; }
    if (!seen) out.push_back(r.catalogue);
  }
  return out;
}

}  // namespace sed
