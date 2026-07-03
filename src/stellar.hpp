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
};

struct StellarMCResult {
  ModeHDI R, M, L;
  bool valid = false;
};

StellarMCResult stellar_mc(const StellarMCInput& in);

}  // namespace sed
