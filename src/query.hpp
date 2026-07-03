// Photometry orchestrator (port of query_photometry.sl). Given (ra, dec),
// issues one bulk VizieR VOTable request for all in-declination catalogues,
// plus separate TAP requests for the non-VizieR surveys the reference stars
// need (MAST PS1_DR2, UKIDSS/WSA, NOAO datalab DELVE), applies the per-catalog
// quality-flag logic, and assembles a PhotometryTable.
//
// Not ported (no reference star needs them): Gaia XP / IUE / MAST box-filter
// spectra, SIMBAD name resolution (we always pass ra/dec), and the many
// surveys absent from the three validation fixtures. The per-catalog dispatch
// is structured so those can be added later.
#pragma once

#include "photometry_table.hpp"

namespace sed {

struct QueryOptions {
  double search_radius_deg = 0.01;  // qualifier search_radius
  double disk_b = 10.0;             // WISE disk-blending latitude cut (deg)
  double disk_l = 100.0;            // WISE disk-blending longitude cut (deg)
  bool force_search_radius = false;
  bool skip_tap = false;  // VizieR only (skip PS1_DR2/WSA/DELVE); for testing
};

// Build a photometry table for the object at (ra, dec) in decimal degrees
// (J2000). Reddening is NOT filled here (query_reddening is a separate step,
// mirroring the S-Lang where the in-function call is commented out).
PhotometryTable query_photometry(double ra, double dec,
                                 const QueryOptions& opt = {});

// Galactic coordinates (l, b) in degrees from equatorial (ra, dec) J2000.
void radec2lb(double ra, double dec, double& l, double& b);

// Great-circle angular separation in degrees (port of angular_separation.sl).
double angular_separation_deg(double ra1, double dec1, double ra2, double dec2);

}  // namespace sed
