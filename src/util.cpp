#include "util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

namespace sed {

double interp_linear(double x, const dvec& xp, const dvec& yp) {
  const size_t n = xp.size();
  // ISIS interpolate_d: bsearch for n0 with xp[n0] <= x (clamped), then
  // linear interpolation/extrapolation using points n0, n0+1 (or n0-1 at the
  // upper end).
  size_t n0;
  if (x < xp[0]) {
    n0 = 0;
  } else if (x >= xp[n - 1]) {
    n0 = n - 1;
  } else {
    n0 = std::upper_bound(xp.begin(), xp.end(), x) - xp.begin() - 1;
  }
  double x0 = xp[n0];
  size_t n1 = n0 + 1;
  if (x == x0) return yp[n0];
  if (n1 == n) {
    if (n0 == 0) return yp[n0];
    n1 = n0 - 1;
  }
  double x1 = xp[n1];
  if (x1 == x0) return yp[n0];
  return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}

dvec interp_linear(const dvec& x, const dvec& xp, const dvec& yp) {
  dvec out(x.size());
  for (size_t i = 0; i < x.size(); ++i) out[i] = interp_linear(x[i], xp, yp);
  return out;
}

double integrate_trapez(const dvec& x, const dvec& y) {
  double s = 0.0;
  for (size_t i = 1; i < x.size(); ++i)
    s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return s;
}

dvec sorted_union(const dvec& a, const dvec& b) {
  dvec out;
  out.reserve(a.size() + b.size());
  out.insert(out.end(), a.begin(), a.end());
  out.insert(out.end(), b.begin(), b.end());
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

double quantile(double p, dvec a) {
  const size_t n = a.size();
  std::sort(a.begin(), a.end());
  double idx = std::max(0.0, std::min(p * (n - 1), double(n - 1)));
  return a[size_t(idx)];
}

double median(dvec a) {
  // ISIS find_median: element (n-1)/2 of the sorted array (lower middle).
  const size_t k = (a.size() - 1) / 2;
  std::nth_element(a.begin(), a.begin() + k, a.end());
  return a[k];
}

std::string trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

std::vector<std::string> split_ws(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream is(s);
  std::string tok;
  while (is >> tok) out.push_back(tok);
  return out;
}

}  // namespace sed
