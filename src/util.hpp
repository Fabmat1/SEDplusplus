// Small numeric helpers shared across modules. Semantics follow the ISIS
// counterparts (linear interpolation with linear extrapolation, S-Lang
// 'union' of sorted arrays, trapezoidal integration of a piecewise-linear
// function) so that fit results are comparable to the S-Lang pipeline.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace sed {

using dvec = std::vector<double>;

// ISIS interpolate_d: linear interpolation, linear extrapolation at both ends.
double interp_linear(double x, const dvec& xp, const dvec& yp);
dvec interp_linear(const dvec& x, const dvec& xp, const dvec& yp);

// Trapezoidal integral over the full range of (x, y).
double integrate_trapez(const dvec& x, const dvec& y);

// S-Lang union(): sorted array of unique values from both inputs.
dvec sorted_union(const dvec& a, const dvec& b);

// S-Lang quantile(): sorted a[int(clamp(p*(n-1)))]  (no interpolation).
double quantile(double p, dvec a);

// ISIS median(): middle element for odd n, average of the two middle for even.
double median(dvec a);

std::string trim(const std::string& s);
std::vector<std::string> split_ws(const std::string& s);

}  // namespace sed
