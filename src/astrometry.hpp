// Live Gaia DR3 astrometry query (port of query_astrometry.sl as invoked by
// templates/photometry.sl:213-215 with `quality_warnings, corrected_values,
// search_radius=4/3600, columns=designation,parallax,parallax_error,ruwe,
// ipd_gof_harmonic_amplitude,ipd_frac_multi_peak,visibility_periods_used`).
//
// The S-Lang, when `corrected_values`+DR3 and the ARI Heidelberg TAP
// (https://gaia.ari.uni-heidelberg.de/tap/) answers, queries
// `gaiadr3.gaia_source INNER JOIN gaiadr3.gaia_source_corrections` to obtain a
// server-side Lindegren-2020 zero-point-corrected `parallax_corr`, and computes
// the El-Badry/Rix/Heintz-2021 (Eq. 16) inflated `parallax_error_corr` locally.
// It then applies the parallax NaN-masking of photometry.sl:229-248 and, when
// the corrected pair is valid, uses it as the final parallax/parallax_error.
// If Heidelberg is offline it falls back to gaia.aip.de then the ESA Gaia TAP
// with no `parallax_corr` (raw values used).
//
// See PHASE2.md "Phase 3" for the full provenance and correction formulas.
#pragma once

#include <optional>
#include <string>

namespace sed {

struct AstrometryResult {
  std::string designation;
  double ra = 0, dec = 0;  // coordinates of the matched source
  double parallax = 0;         // final (corrected if available), NaN-masked
  double parallax_error = 0;   // final (corrected if available)
  std::optional<double> ruwe, ipd_gof_harmonic_amplitude,
      visibility_periods_used;
  bool good = false;  // parallax not-NaN and > 0 after masking
};

// Query the Gaia DR3 astrometry for the source closest to (ra, dec) in decimal
// degrees. Returns std::nullopt if every TAP service is unreachable or no row
// is returned. search_radius_deg defaults to 4 arcsec (the template value).
std::optional<AstrometryResult> query_astrometry(
    double ra, double dec, double search_radius_deg = 4.0 / 3600.0);

}  // namespace sed
