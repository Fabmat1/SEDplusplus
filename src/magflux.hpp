// magnitudes_to_flux (src/plotting/xfig_photometry.sl:3-185): convert the
// observed magnitudes/colours to flux data points relative to a model spectrum
// (l, f), for the photometry_fit_mag.txt / photometry_fit_col.txt products.
// Only the branch actually taken by photometry.sl:1347 is ported (theta +
// extinction scaling, errbar and flagged qualifiers, VizieR_catalog present).
#pragma once

#include <string>
#include <vector>

#include "passbands.hpp"
#include "photometry_table.hpp"
#include "util.hpp"

namespace sed {

// One row of mout.mag (struct-field order = print_struct column order).
struct MagFluxRow {
  double lambda_min, lambda, lambda_max;
  double flux_min, flux, flux_max;
  double diff, diff_err;
  std::string passband, system;
  int flag;
  std::string vizier_catalog;
};

// One row of mout.col.
struct ColFluxRow {
  double diff, diff_err;
  std::string passband, system;
  int flag;
  std::string vizier_catalog;
};

struct MagFluxResult {
  std::vector<MagFluxRow> mag;
  std::vector<ColFluxRow> col;
};

// Ports magnitudes_to_flux(l, f, entries; theta, E_44m55, R_55, errbar,
// flagged). `f` is the unscaled/unreddened model flux on `l`; it is scaled by
// (theta/2)^2 and reddened internally. `entries` is a copy of
// photo.photometric_entries (all rows carry a VizieR_catalog). `db` supplies
// the ZP table, filter curves and colour definitions.
MagFluxResult magnitudes_to_flux(const dvec& l, dvec f,
                                 std::vector<PhotEntry> entries,
                                 const PassbandDB& db, double theta,
                                 double E_44m55, double R_55, bool errbar,
                                 bool flagged);

}  // namespace sed
