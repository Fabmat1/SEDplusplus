// Monte-Carlo propagation to stellar parameters (R, M, L) and the
// mode_and_HDI statistic (iterative histogram with Savitzky-Golay smoothing),
// ported from the photometry.sl template and mode_and_HDI.sl.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "util.hpp"

namespace sed {

struct ModeHDI {
  double mode, HDI_lo, HDI_hi;
  double median, quantile_lo, quantile_hi;
  double p_s, p_a;
};

ModeHDI mode_and_HDI(const dvec& d, double p = 0.6827);

// Savitzky-Golay smoothing (non-periodic edges), NR3 sect. 14.9.
dvec savitzky_golay_smoothing(const dvec& data, int nl, int nr, int p);

// Asymmetric-Gaussian MC sampler: value + N(0,1)*delta_minus (negative draws)
// or *delta_plus (positive draws), as used throughout the template.
class AsymSampler;

struct StellarMCInput {
  double logtheta, d_logtheta_minus, d_logtheta_plus;
  double teff, d_teff_minus, d_teff_plus;
  double logg, d_logg_minus, d_logg_plus;
  double sur_ratio, d_sur_minus, d_sur_plus;
  double parallax, parallax_error;  // mas
  double confidence = 0.6827;       // HDI/quantile coverage (p=confidence)
  long n_mc;
  unsigned long seed;
  // Keep the filtered R/L/M draws in the result (for the SED_results.fits
  // MC_c* extensions, photometry.sl:1234-1237). Off in bulk mode: 3x n_mc
  // doubles per component.
  bool retain_arrays = false;
};

struct StellarMCResult {
  ModeHDI R, M, L;
  // v_grav = GM/(Rc) and v_esc = sqrt(2 g R) in km/s (photometry.sl:
  // 1176-1177), computed from the same draws as R/M/L (after the
  // parallax/sur_ratio filter) and filtered >0 independently.
  ModeHDI vgrav, vesc;
  bool valid = false;
  dvec R_arr, L_arr, M_arr;  // filtered draws; only when retain_arrays
};

StellarMCResult stellar_mc(const StellarMCInput& in);

// Gaia distance MC (photometry.sl:936-945): mode_and_HDI of
// 1e3/(parallax + grand()*parallax_error), own RNG stream.
ModeHDI gaia_distance_mc(double parallax, double parallax_error,
                         double confidence, long n_mc, unsigned long seed);

}  // namespace sed
