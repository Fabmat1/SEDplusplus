// Photometry orchestrator (port of query_photometry.sl). Given (ra, dec),
// issues one bulk VizieR VOTable request for all in-declination catalogues,
// plus separate TAP requests for the non-VizieR surveys the reference stars
// need (MAST PS1_DR2, UKIDSS/WSA, NOAO datalab DELVE), applies the per-catalog
// quality-flag logic, appends box-filter pseudo-photometry from spectra
// (Gaia XP / IUE / MAST, see spectra.hpp), and assembles a PhotometryTable.
//
// Not ported: SIMBAD name resolution lives in simbad.hpp (callers pass
// ra/dec); a few surveys absent from all validation fixtures remain
// unimplemented TAP branches. The per-catalog dispatch is structured so
// these slot in as additional branches.
#pragma once

#include "photometry_table.hpp"
#include "votable.hpp"

namespace sed {

struct QueryOptions {
  double search_radius_deg = 0.01;  // qualifier search_radius
  double disk_b = 10.0;             // WISE disk-blending latitude cut (deg)
  double disk_l = 100.0;            // WISE disk-blending longitude cut (deg)
  bool force_search_radius = false;
  bool skip_tap = false;  // VizieR only (skip PS1_DR2/WSA/DELVE); for testing
  // Spectra -> box-filter rows. The template's standard call passes the IUE
  // and MAST qualifiers, so both default on; Gaia XP is attempted regardless
  // (as in ISIS). spectra_dir is the parent of the IUE/ and MAST/ caches
  // (ISIS uses $cwd; callers pass the output file's directory).
  bool iue = true;
  bool mast = true;
  std::string spectra_dir = ".";
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

// Single/multi-catalogue VizieR VOTable fetch (mirror fall-over); used by the
// bulk catalogue query and by the IUE spectra path (VI/110/inescat).
votable::Document vizier_fetch(const std::string& source, double ra, double dec,
                               double maxr_arcsec);

}  // namespace sed
