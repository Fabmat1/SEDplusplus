#include "fitfun.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "extinction.hpp"

namespace sed {

namespace {
// S-Lang %.Nf rounding used when building the parameter defaults
double roundf_fmt(double v, const char* fmt) {
  char buf[64];
  std::snprintf(buf, sizeof buf, fmt, v);
  return std::atof(buf);
}

bool name_match(const std::string& pattern, const std::string& name) {
  // glob with '*' only
  size_t pi = 0, ni = 0, star = std::string::npos, mark = 0;
  while (ni < name.size()) {
    if (pi < pattern.size() &&
        (pattern[pi] == name[ni] || pattern[pi] == '?')) {
      ++pi;
      ++ni;
    } else if (pi < pattern.size() && pattern[pi] == '*') {
      star = pi++;
      mark = ni;
    } else if (star != std::string::npos) {
      pi = star + 1;
      ni = ++mark;
    } else {
      return false;
    }
  }
  while (pi < pattern.size() && pattern[pi] == '*') ++pi;
  return pi == pattern.size();
}
}  // namespace

FitFunction::FitFunction(std::vector<std::shared_ptr<ModelGrid>> grids,
                         const PassbandDB& db, int n_dummy)
    : grids_(std::move(grids)), db_(db) {
  const size_t len = grids_.size();
  // NOTE: the spectrum cache key is fixed-precision (rounded parameters), so
  // the capacity is part of the observable behaviour: a larger cache reuses
  // spectra for parameter values that differ below the key precision where a
  // smaller one recomputes them. Keep the ISIS value of 6 per component.
  cache_capacity_ = len * 6;
  auto add = [&](const std::string& name, double value, int freeze,
                 double hmin, double hmax, double mn, double mx, double step,
                 double relstep) {
    Param p;
    p.name = name;
    p.value = value;
    p.freeze = freeze;
    p.hard_min = hmin;
    p.hard_max = hmax;
    p.min = mn;
    p.max = mx;
    p.step = step;
    p.relstep = relstep;
    params_.push_back(std::move(p));
  };
  add("logtheta", -9, 0, -DBL_MAX, 5, -18, -2, 0, 0.0001);
  add("E_44m55", 0, 0, 0, 100, 0, 10, 0, 0.001);
  add("R_55", 3.02, 0, 2.5, 6, 2.5, 6, 0, 0.001);
  for (size_t i = 0; i < len; ++i) {
    const GridCoverage& cov = grids_[i]->coverage();
    std::string c = "c" + std::to_string(i + 1) + "_";
    add(c + "vsini", 0, 1, 0, 1000, 0, 1000, 1, 0);
    add(c + "zeta", 0, 1, 0, 1000, 0, 1000, 1, 0);
    add(c + "vrad", 0, 1, -3e5, 3e5, -3e5, 3e5, 0.5, 0);
    auto add_grid_par = [&](const std::string& nm, const dvec& t,
                            const char* fmt, double relstep,
                            bool clamp_hard_min_zero) {
      if (t.size() == 1) {
        double v = roundf_fmt(t[0], fmt);
        add(nm, v, 1, v, v, v, v, 0, relstep);
      } else {
        double v = roundf_fmt(0.5 * (t.front() + t.back()), fmt);
        double hmin = roundf_fmt(t.front() - (t[1] - t.front()), fmt);
        if (clamp_hard_min_zero) hmin = std::max(hmin, 0.0);
        double hmax = roundf_fmt(t.back() + (t.back() - t[t.size() - 2]), fmt);
        add(nm, v, 0, hmin, hmax, roundf_fmt(t.front(), fmt),
            roundf_fmt(t.back(), fmt), 0, relstep);
      }
    };
    add_grid_par(c + "teff", cov.t, "%.0f", 0.0001, true);
    add_grid_par(c + "logg", cov.g, "%.3f", 0.001, true);
    add_grid_par(c + "xi", cov.x, "%.2f", 0.001, true);
    add_grid_par(c + "z", cov.z, "%.2f", 0.001, false);
    add_grid_par(c + "HE", cov.HHE, "%.3f", 0.001, false);
    if (i == 0)
      add(c + "sur_ratio", 1, 1, 1, 1, 1, 1, 0, 0.0001);
    else
      add(c + "sur_ratio", 1, 0, 0, 1e10, 0, 1500, 0, 0.0001);
  }
  add("bb_teff", 0, 1, 0, 1e10, 0, 1e6, 0, 0.0001);
  add("bb_sur_ratio", 0, 1, 0, 1e10, 0, 1500, 0, 0.0001);
  add("bb2_teff", 0, 1, 0, 1e10, 0, 1e6, 0, 0.0001);
  add("bb2_sur_ratio", 0, 1, 0, 1e10, 0, 1500, 0, 0.0001);
  for (int d = 1; d <= n_dummy; ++d)
    add("dummy_" + std::to_string(d), 0, 1, -DBL_MAX, DBL_MAX, 0, 0, 0, 0);
}

int FitFunction::index_of(const std::string& name) const {
  for (size_t i = 0; i < params_.size(); ++i)
    if (params_[i].name == name) return int(i);
  return -1;
}

// The three overloads replicate ISIS do_set_par exactly: min/max/freeze (and
// step/relstep) are applied first via Fit_set_param_control; the value is
// assigned only when it lies inside the (possibly just-updated) [min,max]
// interval, otherwise it is rejected (the parameter keeps its previous value)
// and set_par throws so the caller can warn -- but the min/max/freeze changes
// have already persisted. A min/max that violates the hard limits is a hard
// failure that applies nothing.
void FitFunction::set_par(const std::string& name, double value) {
  bool matched = false, value_rejected = false;
  for (auto& p : params_) {
    if (name_match(name, p.name)) {
      matched = true;
      if (value >= p.min && value <= p.max)
        p.value = value;
      else
        value_rejected = true;
    }
  }
  if (!matched) throw std::runtime_error("no parameter matches " + name);
  if (value_rejected)
    throw std::runtime_error("value " + std::to_string(value) + " for " + name +
                             " lies outside [min,max]");
}

void FitFunction::set_par(const std::string& name, double value, int freeze) {
  bool matched = false, value_rejected = false;
  for (auto& p : params_) {
    if (name_match(name, p.name)) {
      matched = true;
      if (freeze >= 0) p.freeze = freeze;
      if (value >= p.min && value <= p.max)
        p.value = value;
      else
        value_rejected = true;
    }
  }
  if (!matched) throw std::runtime_error("no parameter matches " + name);
  if (value_rejected)
    throw std::runtime_error("value " + std::to_string(value) + " for " + name +
                             " lies outside [min,max]");
}

void FitFunction::set_par(const std::string& name, double value, int freeze,
                          double mn, double mx) {
  bool matched = false, value_rejected = false;
  for (auto& p : params_) {
    if (name_match(name, p.name)) {
      matched = true;
      // min/max must lie within the hard limits; a (0,0) range means "use the
      // hard limits". A violation here is a hard failure (nothing applied to
      // this parameter).
      double nmin = mn, nmax = mx;
      if (nmin == 0.0 && nmax == 0.0) {
        nmin = p.hard_min;
        nmax = p.hard_max;
      }
      if (nmin < p.hard_min || nmin > p.hard_max || nmax < p.hard_min ||
          nmax > p.hard_max)
        throw std::runtime_error("min/max for " + p.name +
                                 " violate hard limits");
      p.min = nmin;
      p.max = nmax;
      if (freeze >= 0) p.freeze = freeze;
      if (value >= p.min && value <= p.max)
        p.value = value;
      else
        value_rejected = true;
    }
  }
  if (!matched) throw std::runtime_error("no parameter matches " + name);
  if (value_rejected)
    throw std::runtime_error("value " + std::to_string(value) + " for " + name +
                             " lies outside [min,max]");
}

void FitFunction::mark_fun(const std::string& name) {
  int i = index_of(name);
  if (i < 0) throw std::runtime_error("no parameter " + name);
  params_[i].is_fun = true;
  params_[i].min = -DBL_MAX;
  params_[i].max = DBL_MAX;
  params_[i].freeze = 1;
}

const dvec& FitFunction::model_flux(const std::vector<double>& par,
                                    dvec& scratch) const {
  const size_t len = grids_.size();
  // gather per-component spectra
  if (len == 1) {
    const int base = 3;
    const dvec& f = grids_[0]->interpolate(par[base + 3], par[base + 4],
                                           par[base + 5], par[base + 6],
                                           par[base + 7], cache_capacity_);
    check_bb(par);
    return f;
  } else {
    dvec& f_uni = scratch;
    f_uni.assign(l_uni_.size(), 0.0);
    for (size_t i = 0; i < len; ++i) {
      const int base = 3 + int(i) * 9;
      const double sur = par[base + 8];
      const dvec& f = grids_[i]->interpolate(par[base + 3], par[base + 4],
                                             par[base + 5], par[base + 6],
                                             par[base + 7], cache_capacity_);
      const Resample& rs = resample_[i];
      for (size_t k = 0; k < l_uni_.size(); ++k) {
        double fv = f[rs.j0[k]] + (f[rs.j1[k]] - f[rs.j0[k]]) * rs.t[k];
        f_uni[k] += sur * fv;
      }
    }
    check_bb(par);
    return f_uni;
  }
}

void FitFunction::check_bb(const std::vector<double>& par) const {
  // black-body components are not supported in this version
  const size_t bb_base = 3 + grids_.size() * 9;
  if ((par[bb_base] != 0 && par[bb_base + 1] != 0) ||
      (par[bb_base + 2] != 0 && par[bb_base + 3] != 0))
    throw std::runtime_error("black-body components not implemented yet");
}

void FitFunction::prepare(const std::vector<int>& mag_ind) {
  const size_t len = grids_.size();
  if (len == 1) {
    l_uni_ = grids_[0]->lambda();
  } else {
    double lmin = 0.0, lmax = 0.0;
    dvec u;
    for (size_t i = 0; i < len; ++i) {
      const dvec& l = grids_[i]->lambda();
      lmin = std::max(lmin, l.front());
      lmax = (i == 0) ? l.back() : std::min(lmax, l.back());
      u = (i == 0) ? l : sorted_union(u, l);
    }
    l_uni_.clear();
    for (double x : u)
      if (lmin <= x && x <= lmax) l_uni_.push_back(x);
    // precompute resampling of each component onto l_uni
    resample_.resize(len);
    for (size_t i = 0; i < len; ++i) {
      const dvec& l = grids_[i]->lambda();
      Resample& rs = resample_[i];
      rs.j0.resize(l_uni_.size());
      rs.j1.resize(l_uni_.size());
      rs.t.resize(l_uni_.size());
      for (size_t k = 0; k < l_uni_.size(); ++k) {
        double x = l_uni_[k];
        size_t n0;
        const size_t n = l.size();
        if (x < l[0])
          n0 = 0;
        else if (x >= l[n - 1])
          n0 = n - 1;
        else
          n0 = std::upper_bound(l.begin(), l.end(), x) - l.begin() - 1;
        size_t n1 = n0 + 1;
        if (x == l[n0] || (n1 == n && n0 == 0)) {
          rs.j0[k] = rs.j1[k] = int(n0);
          rs.t[k] = 0;
          continue;
        }
        if (n1 == n) n1 = n0 - 1;
        rs.j0[k] = int(n0);
        rs.j1[k] = int(n1);
        rs.t[k] = (x - l[n0]) / (l[n1] - l[n0]);
      }
    }
  }
  extinction_curve_components(l_uni_, ext_k0_, ext_s_);
  synth_ = std::make_unique<SynthMag>(db_, l_uni_, mag_ind);
  ext_E_ = ext_R_ = -1.0;  // l_uni changed: invalidate the extinction cache
}

void FitFunction::evaluate(const std::vector<double>& par, dvec& mags) const {
  static thread_local dvec scratch, f_red;
  const dvec& f_uni = model_flux(par, scratch);
  const double logtheta = par[0], E = par[1], R = par[2];
  if (!(E == ext_E_ && R == ext_R_)) {
    ext_fac_.resize(l_uni_.size());
    for (size_t k = 0; k < l_uni_.size(); ++k) {
      double kk = ext_k0_[k] + ext_s_[k] * (R - 3.02);
      ext_fac_[k] = std::pow(10.0, (-0.4 * E) * (kk + R));
    }
    ext_E_ = E;
    ext_R_ = R;
  }
  f_red.resize(f_uni.size());
  for (size_t k = 0; k < f_uni.size(); ++k) f_red[k] = f_uni[k] * ext_fac_[k];
  synth_->magnitudes(f_red, mags);
  const double theta_term = 1.505149978319906 - 5.0 * logtheta;
  for (double& m : mags) m += theta_term;  // NaNs stay NaN
  db_.compute_colors(mags);
  for (size_t i = 0; i < mags.size(); ++i) mags[i] += db_.entries()[i].ZP;
}

}  // namespace sed
