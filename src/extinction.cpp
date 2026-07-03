#include "extinction.hpp"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <cmath>
#include <memory>

#include "f19_table.inc"

namespace sed {

namespace {

struct Splines {
  gsl_spline* k0;
  gsl_spline* s;
  gsl_interp_accel* acc_k0;
  gsl_interp_accel* acc_s;
  Splines() {
    k0 = gsl_spline_alloc(gsl_interp_akima, 102);
    s = gsl_spline_alloc(gsl_interp_akima, 102);
    gsl_spline_init(k0, F19_X, F19_K0, 102);
    gsl_spline_init(s, F19_X, F19_S, 102);
    acc_k0 = gsl_interp_accel_alloc();
    acc_s = gsl_interp_accel_alloc();
  }
};

const Splines& splines() {
  static Splines sp;
  return sp;
}

}  // namespace

void extinction_curve_components(const dvec& lambda, dvec& k0, dvec& s) {
  const Splines& sp = splines();
  const size_t n = lambda.size();
  k0.resize(n);
  s.resize(n);
  for (size_t i = 0; i < n; ++i) {
    double x = 1.0e4 / lambda[i];
    if (x < 8.7) {
      k0[i] = gsl_spline_eval(sp.k0, x, sp.acc_k0);
      s[i] = gsl_spline_eval(sp.s, x, sp.acc_s);
    } else {
      // linear extrapolation from the last two table points, as in the
      // S-Lang code (which extrapolates k0 + s*(R-3.02) jointly; splitting
      // into components is exact because the extrapolation is linear).
      double t = (x - 8.600) / (8.700 - 8.600);
      k0[i] = 8.264 + (8.464 - 8.264) * t;
      s[i] = -1.609 + (-1.650 - -1.609) * t;
    }
  }
}

dvec extinction_factor(const dvec& lambda, double E_44m55, double R_55) {
  dvec k0, s;
  extinction_curve_components(lambda, k0, s);
  dvec out(lambda.size());
  for (size_t i = 0; i < lambda.size(); ++i) {
    double k = k0[i] + s[i] * (R_55 - 3.02);
    out[i] = std::pow(10.0, (-0.4 * E_44m55) * (k + R_55));
  }
  return out;
}

}  // namespace sed
