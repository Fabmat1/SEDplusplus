#include "spectrum.hpp"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <cmath>
#include <stdexcept>

namespace sed {

namespace {

// S-Lang compute_range_multiplier (slarray.c), used by the [a:b:dx] range.
double range_multiplier(double xmin, double dx) {
  double multiplier = 1.0;
  while (multiplier < 1e9) {
    volatile double xmin1 = multiplier * xmin;
    volatile double dx1 = multiplier * dx;
    if ((xmin1 == (double)(int)xmin1) && (dx1 == (double)(int)dx1))
      return multiplier;
    multiplier *= 10.0;
  }
  return 1.0;
}

// S-Lang inline_implicit_floating_array for [xmin : xmax : dx] (dx>0).
dvec slang_range(double xmin, double xmax, double dx) {
  dvec out;
  if ((xmax <= xmin && dx >= 0.0) || (xmax >= xmin && dx <= 0.0)) return out;
  int n;
  if ((xmin + dx == (volatile double)xmin) || (xmax + dx == (volatile double)xmax))
    return out;
  n = (int)(1.5 + ((xmax - xmin) / dx));
  if (n <= 0) return out;
  double multiplier = range_multiplier(xmin, dx);
  double last = ((xmin * multiplier) + (n - 1) * (dx * multiplier)) / multiplier;
  if (dx > 0.0) {
    if (last >= xmax) n -= 1;
  } else if (last <= xmax)
    n -= 1;
  out.resize(n);
  if (multiplier != 1.0) {
    int ixmin = (int)std::floor(multiplier * xmin + 0.5);
    int idx = (int)std::floor(multiplier * dx + 0.5);
    for (int i = 0; i < n; ++i)
      out[i] = (ixmin + (double)i * idx) / multiplier;
  } else {
    for (int i = 0; i < n; ++i) out[i] = xmin + i * dx;
  }
  return out;
}

// S-Lang [xmin : xmax : #nels]
dvec slang_range_n(double xmin, double xmax, int nels) {
  dvec out;
  if (nels <= 0) return out;
  double dx = (nels == 1) ? 0.0 : (xmax - xmin) / (nels - 1);
  out.resize(nels);
  for (int i = 0; i < nels; ++i) out[i] = xmin + i * dx;
  if (nels > 1) out[nels - 1] = xmax;
  return out;
}

// Port of c_functions->convolve (convolve_single) exactly, including the
// summation order (j from length_g-1 down to 0) and the flat-edge handling.
void convolve_single(const dvec& cx, dvec& cy, const dvec& factor,
                     const dvec& fx, const dvec& fy, const dvec& gx,
                     const dvec& gy) {
  const int length_c = int(cx.size());
  const int length_f = int(fx.size());
  const int length_g = int(gx.size());
  cy.assign(length_c, 0.0);
  int k = 0;
  {  // binary_search for cx[0] + factor[0]*gx[0]
    double x = cx[0] + factor[0] * gx[0];
    int lo = 0, hi = length_f - 1;
    if (x <= fx[0])
      k = 0;
    else if (x >= fx[hi])
      k = hi;
    else {
      while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (fx[mid] <= x)
          lo = mid;
        else
          hi = mid;
      }
      k = lo;
    }
  }
  for (int i = 0; i < length_c; ++i) {
    double sum = 0.0;
    double x = cx[i] - factor[i] * gx[length_g - 1];
    while (x < fx[k] && k > 0) k--;
    for (int j = length_g - 1; j >= 0; --j) {
      x = cx[i] - factor[i] * gx[j];
      if (x < fx[0])
        sum += fy[0] * gy[j];
      else if (x > fx[length_f - 1])
        sum += fy[length_f - 1] * gy[j];
      else {
        while (x > fx[k + 1] && k + 1 < length_f - 1) k++;
        sum += (fy[k] + ((fy[k + 1] - fy[k]) * (x - fx[k])) / (fx[k + 1] - fx[k])) *
               gy[j];
      }
    }
    cy[i] = sum;
  }
}

// GSL Akima spline evaluation of (x, y) at points xi (all within [x0, xN]).
dvec akima_eval(const dvec& x, const dvec& y, const dvec& xi) {
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* sp = gsl_spline_alloc(gsl_interp_akima, x.size());
  gsl_spline_init(sp, x.data(), y.data(), x.size());
  dvec out(xi.size());
  for (size_t i = 0; i < xi.size(); ++i)
    out[i] = gsl_spline_eval(sp, xi[i], acc);
  gsl_spline_free(sp);
  gsl_interp_accel_free(acc);
  return out;
}

}  // namespace

void convolve_syn_lowres(const dvec& lambda, const dvec& flux_in,
                         double res_slope, dvec& olambda, dvec& cflux) {
  const double fwhm_to_sigma = 0.4246609001440095;
  const double max_spacing = 0.015;
  double med = median(flux_in);
  dvec flux(flux_in.size());
  for (size_t i = 0; i < flux.size(); ++i) flux[i] = flux_in[i] / med;

  const double res_slope_times2 = 2.0 * res_slope;
  const double step = 0.5 / res_slope_times2;
  const double lam0 = lambda.front(), lamN = lambda.back();
  dvec clambda = slang_range(lam0, lamN + (1.0 - 1e-5) * step, step);
  if (clambda.empty()) throw std::runtime_error("convolve: empty clambda");
  clambda.back() = lamN;

  const int clen = int(clambda.size());
  dvec sig(clen);
  double maxsig = 0.0;
  for (int i = 0; i < clen; ++i) {
    sig[i] = clambda[i] / (res_slope * clambda[i]) * fwhm_to_sigma;
    if (sig[i] > maxsig) maxsig = sig[i];
  }
  int N_gpro_needed = int(std::lround(maxsig / max_spacing));
  if (N_gpro_needed < 15) N_gpro_needed = 15;
  dvec glambda = slang_range_n(-3.0, 3.0, 6 * N_gpro_needed + 1);
  dvec gweights(glambda.size());
  double gsum = 0.0;
  for (size_t i = 0; i < glambda.size(); ++i) {
    gweights[i] = std::exp(-glambda[i] * glambda[i] / 2.0);
    gsum += gweights[i];
  }
  for (double& w : gweights) w /= gsum;

  dvec cfl;
  convolve_single(clambda, cfl, sig, lambda, flux, glambda, gweights);
  for (double& v : cfl) v *= med;

  olambda = sorted_union(lambda, clambda);
  cflux = akima_eval(clambda, cfl, olambda);
}

// Convolve one component (its own grid) and return (olambda, cflux).
static void component_spectrum(const FitFunction& fun, size_t i, double xmin,
                               double xmax, double res_slope, dvec& olambda,
                               dvec& cflux) {
  const auto& p = fun.params();
  const int base = 3 + int(i) * 9;
  const ModelGrid& g = fun.grid(i);
  const dvec& gl = g.lambda();
  const dvec& gf = g.interpolate(p[base + 3].value, p[base + 4].value,
                                 p[base + 5].value, p[base + 6].value,
                                 p[base + 7].value, 64);
  dvec lam, flx;
  for (size_t k = 0; k < gl.size(); ++k)
    if (xmin <= gl[k] && gl[k] <= xmax) {
      lam.push_back(gl[k]);
      flx.push_back(gf[k]);
    }
  convolve_syn_lowres(lam, flx, res_slope, olambda, cflux);
}

ModelSpectrum build_model_spectrum(const FitFunction& fun, double xmin,
                                   double xmax, double res_slope) {
  const size_t len = fun.n_components();
  const auto& p = fun.params();
  // black-body components not supported here (none in the reference fits)
  const size_t bb_base = 3 + len * 9;
  if ((p[bb_base].value != 0 && p[bb_base + 1].value != 0) ||
      (p[bb_base + 2].value != 0 && p[bb_base + 3].value != 0))
    throw std::runtime_error("black-body components not implemented");

  std::vector<dvec> ol(len), cf(len);
  double lmin = 0.0, lmax = 0.0;
  dvec uni;
  for (size_t i = 0; i < len; ++i) {
    component_spectrum(fun, i, xmin, xmax, res_slope, ol[i], cf[i]);
    lmin = std::max(lmin, ol[i].front());
    lmax = (i == 0) ? ol[i].back() : std::min(lmax, ol[i].back());
    uni = (i == 0) ? ol[i] : sorted_union(uni, ol[i]);
  }
  ModelSpectrum ms;
  for (double x : uni)
    if (lmin <= x && x <= lmax) ms.l.push_back(x);

  ms.f.assign(ms.l.size(), 0.0);
  if (len > 1) ms.f_comp.assign(len, dvec(ms.l.size(), 0.0));
  for (size_t i = 0; i < len; ++i) {
    const double sur = p[3 + int(i) * 9 + 8].value;
    dvec fi = interp_linear(ms.l, ol[i], cf[i]);
    for (size_t k = 0; k < ms.l.size(); ++k) {
      double v = sur * fi[k];
      ms.f[k] += v;
      if (len > 1) ms.f_comp[i][k] = v;
    }
  }
  return ms;
}

}  // namespace sed
