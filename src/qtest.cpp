// Temporary smoke test for the query-layer foundation (http + votable).
// Issues the same VizieR CGI VOTable bulk request the S-Lang uses and prints
// what we parsed back, so we can confirm the plumbing end-to-end on live data.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "catalog_registry.hpp"
#include "http.hpp"
#include "votable.hpp"

int main(int argc, char** argv) {
  double ra = 145.977, dec = 27.7829;  // star 1
  if (argc >= 3) { ra = std::atof(argv[1]); dec = std::atof(argv[2]); }

  // Registry sanity + dec-range selection at this dec.
  const auto& reg = sed::catalog_registry();
  auto uniq = sed::unique_catalogues();
  int vizier_cats = 0, vizier_rows_in_dec = 0;
  for (const auto& c : uniq)
    if (sed::is_vizier_catalogue(c)) ++vizier_cats;
  for (const auto& r : reg)
    if (sed::is_vizier_catalogue(r.catalogue) && dec >= r.dec_min &&
        dec <= r.dec_max)
      ++vizier_rows_in_dec;
  std::printf("registry: %zu rows, %zu unique catalogues (%d VizieR)\n",
              reg.size(), uniq.size(), vizier_cats);
  std::printf("at dec=%.4f: %d VizieR registry rows in dec-range\n\n", dec,
              vizier_rows_in_dec);

  std::vector<std::string> cats = {
      "II/246/out",      // 2MASS
      "I/350/gaiaedr3",  // Gaia EDR3
      "II/335/galex_ais",// GALEX
      "II/168/ubvmeans", // Johnson
      "II/349/ps1",      // PS1 DR1
      "V/147/sdss12",    // SDSS
  };

  std::string source;
  for (size_t i = 0; i < cats.size(); ++i)
    source += (i ? "," : "") + cats[i];

  double radius_arcsec = 45.0;  // generous; per-catalog filtering comes later
  char params[512];
  std::snprintf(params, sizeof(params),
                "-source=%s&-c=%.6f%%20%+.6f&-c.rs=%.4f&-out.add=_r&-out.all",
                source.c_str(), ra, dec, radius_arcsec);
  std::string url =
      std::string("http://vizier.u-strasbg.fr/viz-bin/votable?") + params;

  std::printf("GET %s\n\n", url.c_str());
  auto resp = sed::http::get(url);
  std::printf("HTTP %ld, %zu bytes\n\n", resp.status, resp.body.size());

  auto doc = sed::votable::parse(resp.body);
  std::printf("parsed %zu tables:\n", doc.tables.size());
  for (const auto& cat : cats) {
    const auto* t = doc.find(cat);
    if (!t) { std::printf("  %-18s : NOT FOUND\n", cat.c_str()); continue; }
    std::printf("  %-18s id=%-22s rows=%zu fields=%zu\n", cat.c_str(),
                t->id.c_str(), t->rows.size(), t->fields.size());
    if (!t->rows.empty()) {
      int ir = t->field_index("_r");
      std::printf("      _r col idx=%d; first row _r=%s\n", ir,
                  ir >= 0 ? t->rows[0][ir].c_str() : "(none)");
    }
  }
  return 0;
}
