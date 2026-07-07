#include "stellar.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace sed {

namespace {
constexpr double Const_pc_cgs = 3.0856775814913674e18;
constexpr double Const_Rsun_cgs = 6.957e10;
constexpr double Const_GMsun_cgs = 1.3271244e26;
constexpr double Const_Teffsun = 5772.0;
constexpr double Const_c = 2.99792458e10;  // cm/s

// solve (p+1)x(p+1) symmetric system by Gauss-Jordan for SG coefficients
dvec sg_coefficients(int nl, int nr, int p) {
  const int m = p + 1;
  std::vector<std::vector<double>> A(m, std::vector<double>(2 * m, 0.0));
  for (int i = 0; i < m; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = 0.0;
      for (int k = -nl; k <= nr; ++k) s += std::pow(double(k), i + j);
      A[i][j] = A[j][i] = s;
    }
  for (int i = 0; i < m; ++i) A[i][m + i] = 1.0;
  // Gauss-Jordan with partial pivoting
  for (int col = 0; col < m; ++col) {
    int piv = col;
    for (int r = col + 1; r < m; ++r)
      if (std::fabs(A[r][col]) > std::fabs(A[piv][col])) piv = r;
    std::swap(A[piv], A[col]);
    double d = A[col][col];
    for (int c = 0; c < 2 * m; ++c) A[col][c] /= d;
    for (int r = 0; r < m; ++r) {
      if (r == col) continue;
      double f = A[r][col];
      for (int c = 0; c < 2 * m; ++c) A[r][c] -= f * A[col][c];
    }
  }
  // row 0 of the inverse (derivative order 0)
  dvec c(nl + nr + 1);
  for (int i = -nl; i <= nr; ++i) {
    double s = 0.0;
    for (int j = 0; j < m; ++j) s += A[0][m + j] * std::pow(double(i), j);
    c[i + nl] = s;
  }
  return c;
}
}  // namespace

dvec savitzky_golay_smoothing(const dvec& data, int nl, int nr, int p) {
  const int len = int(data.size());
  if (nl + nr > len)
    throw std::runtime_error("savitzky_golay: window larger than data");
  dvec coeffs = sg_coefficients(nl, nr, p);
  dvec out(len, 0.0);
  // interior (periodic bulk equals plain convolution away from edges;
  // edge values are overwritten below exactly as in the S-Lang code)
  for (int i = 0; i < len; ++i) {
    double s = 0.0;
    for (int k = -nl; k <= nr; ++k) {
      int idx = i + k;
      // periodic wrap (matches 'shift'); edges get replaced below
      idx = ((idx % len) + len) % len;
      s += data[idx] * coeffs[k + nl];
    }
    out[i] = s;
  }
  // beginning
  for (int i = 0; i <= nl - 1 && i < len; ++i) {
    int nr_min = std::min(nr, len - i - 1);
    dvec c = sg_coefficients(i, nr_min, p);
    double s = 0.0;
    for (int k = 0; k <= nr_min + i; ++k) s += data[k] * c[k];
    out[i] = s;
  }
  // end
  int start = std::max(nl + 1, len - nr - 1);
  for (int i = start; i < len; ++i) {
    dvec c = sg_coefficients(nl, len - 1 - i, p);
    double s = 0.0;
    for (int k = i - nl; k <= len - 1; ++k) s += data[k] * c[k - (i - nl)];
    out[i] = s;
  }
  return out;
}

ModeHDI mode_and_HDI(const dvec& din, double p) {
  dvec d;
  d.reserve(din.size());
  for (double x : din)
    if (!std::isnan(x)) d.push_back(x);
  const size_t n = d.size();
  if (n < 4) throw std::runtime_error("mode_and_HDI: too few points");

  ModeHDI r{};
  r.p_s = p;
  r.median = median(d);
  r.quantile_lo = quantile(0.5 * (1.0 - p), d);
  r.quantile_hi = quantile(1.0 - 0.5 * (1.0 - p), d);

  double dmin = *std::min_element(d.begin(), d.end());
  double dmax = *std::max_element(d.begin(), d.end());
  if (dmin == dmax) {
    r.mode = d[0];
    r.HDI_lo = r.HDI_hi = std::numeric_limits<double>::quiet_NaN();
    return r;
  }
  const double EPS = 2.220446049250313e-16;

  dvec bin_center, h;
  const int iterations = 10;
  double summe = 0.0;
  for (int it = 0; it < iterations; ++it) {
    if (it != iterations - 1) {
      dvec bounds;
      if (it == 0) {
        size_t np = std::min<size_t>(202, n);
        bounds.resize(np);
        double lo = dmin - EPS, hi = dmax + EPS;
        for (size_t i = 0; i < np; ++i)
          bounds[i] = lo + (hi - lo) * double(i) / double(np - 1);
      } else {
        double temp =
            0.5 * std::max(r.HDI_hi - r.HDI_lo, r.quantile_hi - r.quantile_lo);
        size_t np = std::min<size_t>(200, n);
        dvec mid(np);
        double lo = r.HDI_lo - temp, hi = r.HDI_hi + temp;
        for (size_t i = 0; i < np; ++i)
          mid[i] = lo + (hi - lo) * double(i) / double(np - 1);
        bounds = sorted_union({dmin - EPS, dmax + EPS}, mid);
      }
      const size_t nb = bounds.size() - 1;
      bin_center.resize(nb);
      h.assign(nb, 0.0);
      for (size_t i = 0; i < nb; ++i)
        bin_center[i] = 0.5 * (bounds[i] + bounds[i + 1]);
      // histogram: lo[i] <= x < hi[i]
      dvec sorted_d = d;
      std::sort(sorted_d.begin(), sorted_d.end());
      for (size_t i = 0; i < nb; ++i) {
        auto lo_it = std::lower_bound(sorted_d.begin(), sorted_d.end(),
                                      bounds[i]);
        auto hi_it = (i + 1 == nb)
                         ? std::upper_bound(sorted_d.begin(), sorted_d.end(),
                                            bounds[i + 1])
                         : std::lower_bound(sorted_d.begin(), sorted_d.end(),
                                            bounds[i + 1]);
        h[i] = double(hi_it - lo_it) / double(n);
      }
      // remove overflow bins
      bin_center.erase(bin_center.begin());
      bin_center.pop_back();
      h.erase(h.begin());
      h.pop_back();
    } else {
      // final iteration: smooth and refine
      double stepsize = bin_center[1] - bin_center[0];
      int nn = int(std::lround(
          0.75 * std::min((r.mode - r.HDI_lo) / stepsize,
                          (r.HDI_hi - r.mode) / stepsize)));
      std::vector<int> ind;
      for (size_t i = 0; i < h.size(); ++i)
        if (h[i] > 0) ind.push_back(int(i));
      if (nn > 0 && int(ind.size()) > 2 * nn) {
        dvec hh(ind.size());
        for (size_t i = 0; i < ind.size(); ++i) hh[i] = h[ind[i]];
        dvec sm = savitzky_golay_smoothing(hh, nn, nn, std::min(4, nn));
        for (size_t i = 0; i < ind.size(); ++i) h[ind[i]] = sm[i];
      }
      const int factor = 100;
      size_t nf = size_t(factor) * bin_center.size();
      dvec fine(nf);
      double lo = bin_center.front(), hi = bin_center.back();
      for (size_t i = 0; i < nf; ++i)
        fine[i] = lo + (hi - lo) * double(i) / double(nf - 1);
      dvec hf = interp_linear(fine, bin_center, h);
      for (double& v : hf) v /= factor;
      bin_center = std::move(fine);
      h = std::move(hf);
    }
    // binary search for the HDI
    size_t imax = 0;
    for (size_t i = 1; i < h.size(); ++i)
      if (h[i] > h[imax]) imax = i;
    const long lenm1 = long(bin_center.size()) - 1;
    long x0_old = 0, x1_old = lenm1, x0_new = long(imax), x1_new = long(imax);
    double y0 = 0.0, y1 = h[imax];
    while (x0_old != x0_new || x1_old != x1_new) {
      double y2 = 0.5 * (y0 + y1);
      x0_old = x0_new;
      x1_old = x1_new;
      x0_new = long(imax);
      x1_new = long(imax);
      while (h[x0_new] > y2 && x0_new > 1) --x0_new;
      while (h[x1_new] > y2 && x1_new < lenm1) ++x1_new;
      summe = 0.0;
      for (long i = x0_new; i <= x1_new; ++i) summe += h[i];
      if (summe < p)
        y1 = y2;
      else
        y0 = y2;
    }
    r.mode = bin_center[imax];
    r.HDI_lo = bin_center[x0_new];
    r.HDI_hi = bin_center[x1_new];
    r.p_a = summe;
  }
  if (std::fabs(r.p_a - r.p_s) > 0.05) {
    r.mode = r.HDI_lo = r.HDI_hi = std::numeric_limits<double>::quiet_NaN();
  }
  return r;
}

StellarMCResult stellar_mc(const StellarMCInput& in) {
  std::mt19937_64 rng(in.seed);
  std::normal_distribution<double> gauss(0.0, 1.0);
  const long n = in.n_mc;

  dvec R, M, L, VG, VE;
  R.reserve(n);
  M.reserve(n);
  L.reserve(n);
  VG.reserve(n);
  VE.reserve(n);
  auto asym = [&](double value, double dminus, double dplus) {
    double g = gauss(rng);
    return value + (g < 0 ? g * dminus : g * dplus);
  };
  for (long i = 0; i < n; ++i) {
    double logtheta = asym(in.logtheta, in.d_logtheta_minus, in.d_logtheta_plus);
    double parallax = in.parallax + gauss(rng) * in.parallax_error;
    double sur = asym(in.sur_ratio, in.d_sur_minus, in.d_sur_plus);
    double logg = asym(in.logg, in.d_logg_minus, in.d_logg_plus);
    double teff = asym(in.teff, in.d_teff_minus, in.d_teff_plus);
    if (!(parallax > 1e-4 && sur > 0)) continue;
    double r = std::pow(10.0, logtheta) / (2e-3 * parallax) * Const_pc_cgs /
               Const_Rsun_cgs;
    r *= std::sqrt(sur);
    double m = std::pow(10.0, logg) * r * r * Const_Rsun_cgs * Const_Rsun_cgs /
               Const_GMsun_cgs;
    double l = r * r * std::pow(teff / Const_Teffsun, 4);
    // photometry.sl:1176-1177 (no new draws -- RNG sequence unchanged)
    double vg = m * Const_GMsun_cgs / (r * Const_Rsun_cgs * Const_c) / 1e5;
    double ve = std::sqrt(2.0 * std::pow(10.0, logg) * r * Const_Rsun_cgs) /
                1e5;
    if (r > 0) R.push_back(r);
    if (m > 0) M.push_back(m);
    if (l > 0) L.push_back(l);
    if (vg > 0) VG.push_back(vg);
    if (ve > 0) VE.push_back(ve);
  }
  StellarMCResult out;
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  out.vgrav = out.vesc = ModeHDI{NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
  if (R.size() < 4 || M.size() < 4 || L.size() < 4) return out;
  out.R = mode_and_HDI(R, in.confidence);
  out.M = mode_and_HDI(M, in.confidence);
  out.L = mode_and_HDI(L, in.confidence);
  if (VG.size() >= 4) out.vgrav = mode_and_HDI(VG, in.confidence);
  if (VE.size() >= 4) out.vesc = mode_and_HDI(VE, in.confidence);
  out.valid = true;
  if (in.retain_arrays) {
    out.R_arr = std::move(R);
    out.L_arr = std::move(L);
    out.M_arr = std::move(M);
  }
  return out;
}

ModeHDI gaia_distance_mc(double parallax, double parallax_error,
                         double confidence, long n_mc, unsigned long seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gauss(0.0, 1.0);
  // photometry.sl:939-940: no positivity filter on the parallax draws.
  dvec d(n_mc);
  for (long i = 0; i < n_mc; ++i)
    d[i] = 1e3 / (parallax + gauss(rng) * parallax_error);
  return mode_and_HDI(d, confidence);
}

}  // namespace sed
