// Interstellar extinction after Fitzpatrick et al. (2019, ApJ, 886, 108),
// mirroring interstellar_extinction.sl: Akima-spline interpolation of the
// Table-3 k(lambda-55) curve in x = 1e4/lambda, linear extrapolation for
// x >= 8.7, extinction factor 10^(-0.4*E44m55*(k + R55)).
#pragma once

#include "util.hpp"

namespace sed {

// Evaluate the two spline components k0(x) and s(x) on a wavelength grid
// (Angstroem). The extinction factor for given (E, R) is then
//   10^(-0.4*E*(k0 + s*(R-3.02) + R)).
void extinction_curve_components(const dvec& lambda, dvec& k0, dvec& s);

// Convenience: full factor for one wavelength grid.
dvec extinction_factor(const dvec& lambda, double E_44m55, double R_55);

}  // namespace sed
