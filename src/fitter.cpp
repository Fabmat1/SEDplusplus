#include "fitter.hpp"

extern "C" {
#include <mpfit.h>
}

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace sed {

namespace {

struct MpfitCtx {
  Fitter* fitter;
  FitFunction* fun;
  std::vector<int>* free_idx;
  std::vector<double>* fullpar;
  const dvec* data;
  const dvec* err;
  const std::vector<int>* ind_n;
  dvec mags;  // scratch
};

int mpfit_objective(int m, int n, double* p, double* fvec, double** /*dvec*/,
                    void* vctx) {
  MpfitCtx* ctx = static_cast<MpfitCtx*>(vctx);
  for (int j = 0; j < n; ++j) (*ctx->fullpar)[(*ctx->free_idx)[j]] = p[j];
  ctx->fun->evaluate(*ctx->fullpar, ctx->mags);
  for (int k = 0; k < m; ++k) {
    double model = ctx->mags[(*ctx->ind_n)[k]];
    fvec[k] = ((*ctx->data)[k] - model) / (*ctx->err)[k];
  }
  return 0;
}

mp_config isis_mpfit_config() {
  mp_config c{};
  c.ftol = 1e-10;
  c.xtol = 1e-10;
  c.gtol = 1e-10;
  c.stepfactor = 100.0;
  c.epsfcn = MP_MACHEP0;
  c.covtol = 1e-14;
  c.maxiter = 200;
  c.maxfev = 0;
  c.nprint = 0;
  return c;
}

}  // namespace

Fitter::Fitter(FitFunction& fun, PhotometryTable& phot)
    : fun_(fun), phot_(phot) {}

int Fitter::num_free_params() const {
  int n = 0;
  for (const auto& p : fun_.params())
    if (p.freeze == 0 && !p.is_fun) ++n;
  return n;
}

void Fitter::model_vector(const std::vector<double>& fullpar, dvec& out) const {
  static thread_local dvec mags;
  fun_.evaluate(fullpar, mags);
  out.resize(ind_n_.size());
  for (size_t k = 0; k < ind_n_.size(); ++k) out[k] = mags[ind_n_[k]];
}

void Fitter::build_dataset() {
  const auto& db = fun_.db();
  data_.resize(ind_m_.size());
  err_.resize(ind_m_.size());
  for (size_t k = 0; k < ind_m_.size(); ++k) {
    const PhotEntry& e = phot_.entries[ind_m_[k]];
    double zperr = db.entries()[ind_n_[k]].ZP_err;
    data_[k] = e.magnitude;
    err_[k] = std::sqrt(e.uncertainty * e.uncertainty + zperr * zperr +
                        norm_chi_red_ * norm_chi_red_);
  }
}

void Fitter::select_data(bool verbose) {
  const auto& db = fun_.db();
  // synthetic magnitudes at the initial parameter values, all filters
  std::vector<int> all_mag(db.n_magnitudes());
  for (size_t i = 0; i < all_mag.size(); ++i) all_mag[i] = int(i);
  fun_.prepare(all_mag);
  std::vector<double> fullpar;
  for (const auto& p : fun_.params()) fullpar.push_back(p.value);
  dvec n_mags;
  fun_.evaluate(fullpar, n_mags);

  ind_m_.clear();
  ind_n_.clear();
  mag_ind_.clear();
  std::vector<int> mag_set;
  for (size_t i = 0; i < phot_.entries.size(); ++i) {
    const PhotEntry& e = phot_.entries[i];
    int j = db.find(e.system, e.passband);
    if (j < 0 || std::isnan(n_mags[j])) {
      if (verbose && e.flag == 0)
        std::fprintf(stderr,
                     "Warning: unable to compute synthetic magnitude for "
                     "%s: %s\n",
                     e.system.c_str(), e.passband.c_str());
      continue;
    }
    if (e.flag != 0) continue;
    ind_m_.push_back(int(i));
    ind_n_.push_back(j);
    if (db.entries()[j].type == "magnitude") {
      mag_set.push_back(j);
    } else {
      for (int d : db.color_deps()[j - int(db.n_magnitudes())])
        mag_set.push_back(d);
    }
  }
  if (ind_m_.empty())
    throw std::runtime_error(
        "Unable to compute synthetic magnitudes/colors for any observed "
        "magnitude/color.");
  std::sort(mag_set.begin(), mag_set.end());
  mag_set.erase(std::unique(mag_set.begin(), mag_set.end()), mag_set.end());
  mag_ind_ = mag_set;
}

double Fitter::fit_once() {
  auto& params = fun_.params();
  std::vector<int> free_idx;
  for (size_t i = 0; i < params.size(); ++i)
    if (params[i].freeze == 0 && !params[i].is_fun) free_idx.push_back(int(i));
  std::vector<double> fullpar(params.size());
  for (size_t i = 0; i < params.size(); ++i) fullpar[i] = params[i].value;

  const int nfree = int(free_idx.size());
  if (nfree == 0) return statistic_at_current();

  std::vector<double> p(nfree);
  std::vector<mp_par> mppar(nfree);
  for (int j = 0; j < nfree; ++j) {
    const Param& par = params[free_idx[j]];
    p[j] = par.value;
    mp_par& ps = mppar[j];
    ps = mp_par{};
    ps.fixed = 0;
    ps.limited[0] = 1;
    ps.limited[1] = 1;
    ps.limits[0] = par.min;
    ps.limits[1] = par.max;
    ps.step = par.step;
    ps.relstep = par.relstep;
    ps.side = 0;
  }

  MpfitCtx ctx{this, &fun_, &free_idx, &fullpar, &data_, &err_, &ind_n_, {}};
  mp_config cfg = isis_mpfit_config();
  mp_result result{};
  mpfit(mpfit_objective, int(data_.size()), nfree, p.data(), mppar.data(),
        &cfg, &ctx, &result);
  for (int j = 0; j < nfree; ++j) params[free_idx[j]].value = p[j];
  return result.bestnorm;
}

double Fitter::statistic_at_current() {
  std::vector<double> fullpar;
  for (const auto& p : fun_.params()) fullpar.push_back(p.value);
  dvec model;
  model_vector(fullpar, model);
  double s = 0.0;
  for (size_t k = 0; k < model.size(); ++k) {
    double r = (data_[k] - model[k]) / err_[k];
    s += r * r;
  }
  return s;
}

double Fitter::refit() { return fit_once(); }
double Fitter::eval_statistic() { return statistic_at_current(); }

FitResults Fitter::run(const ParSet& par, const ParSet& par_full,
                       int conf_level, double remove_outliers, bool verbose) {
  auto& params = fun_.params();
  auto warn = [&](const std::string& msg) {
    if (verbose)
      std::fprintf(stderr, "Warning in 'photometric_fitting': %s\n",
                   msg.c_str());
  };
  // ------- set parameters from qualifiers (set_par, then set_par_full)
  for (size_t i = 0; i < par.name.size(); ++i) {
    try {
      fun_.set_par(par.name[i], par.value[i], par.freeze[i]);
    } catch (const std::exception& e) {
      warn(e.what());
    }
  }
  for (size_t i = 0; i < par_full.name.size(); ++i) {
    try {
      fun_.set_par(par_full.name[i], par_full.value[i], par_full.freeze[i],
                   par_full.min[i], par_full.max[i]);
    } catch (const std::exception& e) {
      warn(e.what());
    }
  }

  // ------- determine fittable data and needed magnitudes
  select_data(verbose);

  // ------- initial guesses (E_44m55 and logtheta)
  {
    std::vector<double> fullpar;
    for (const auto& p : params) fullpar.push_back(p.value);
    dvec n_mags;
    fun_.evaluate(fullpar, n_mags);
    const double theta_term = 1.505149978319906 - 5.0 * fullpar[0];
    const auto& db = fun_.db();
    // n without the theta shift (colors are shift-free by construction)
    dvec n_pure = n_mags;
    for (size_t i = 0; i < db.n_magnitudes(); ++i) n_pure[i] -= theta_term;

    auto in_par = [&](const std::string& nm, const ParSet& ps) {
      for (const auto& s : ps.name)
        if (s == nm) return true;
      return false;
    };
    auto obs_first = [&](const char* sys, const char* pb, double& out) {
      for (const auto& e : phot_.entries)
        if (e.system == sys && e.passband == pb) {
          out = e.magnitude;
          return true;
        }
      return false;
    };
    auto fitted = [&](const char* sys, const char* pb) {
      for (int im : ind_m_) {
        const PhotEntry& e = phot_.entries[im];
        if (e.system == sys && e.passband == pb) return true;
      }
      return false;
    };

    int iE = fun_.index_of("E_44m55");
    if (!(params[iE].freeze == 1 || in_par("E_44m55", par_full) ||
          in_par("E_44m55", par))) {
      dvec est;
      double mv;
      int nBmV = db.find("Johnson", "BmV");
      if (fitted("Johnson", "BmV") && obs_first("Johnson", "BmV", mv))
        est.push_back(mv - n_pure[nBmV]);
      if (fitted("Johnson", "B") && fitted("Johnson", "V")) {
        double mB, mV;
        obs_first("Johnson", "B", mB);
        obs_first("Johnson", "V", mV);
        est.push_back(mB - mV - n_pure[nBmV]);
      }
      int nbmy = db.find("Stroemgren", "bmy");
      if (fitted("Stroemgren", "bmy") && obs_first("Stroemgren", "bmy", mv))
        est.push_back(mv - n_pure[nbmy]);
      if (fitted("Stroemgren", "b") && fitted("Stroemgren", "y")) {
        double mb, my;
        obs_first("Stroemgren", "b", mb);
        obs_first("Stroemgren", "y", my);
        est.push_back(mb - my - n_pure[nbmy]);
      }
      double E = est.empty() ? 0.0 : std::max(0.0, median(est));
      try {
        fun_.set_par("E_44m55", E);
      } catch (const std::exception& e) {
        warn(e.what());
      }
    }
    if (!(params[0].freeze == 1 || in_par("logtheta", par_full) ||
          in_par("logtheta", par))) {
      dvec est;
      for (size_t k = 0; k < ind_m_.size(); ++k) {
        if (db.entries()[ind_n_[k]].type == "magnitude") {
          double m = phot_.entries[ind_m_[k]].magnitude;
          est.push_back(std::log10(2.0) +
                        0.5 * (-0.4 * (m - n_pure[ind_n_[k]])));
        }
      }
      if (est.empty())
        throw std::runtime_error(
            "Unable to compute theta from colors alone. Specify at least one "
            "magnitude.");
      try {
        fun_.set_par("logtheta", median(est));
      } catch (const std::exception& e) {
        warn(e.what());
      }
    }
  }

  // ------- fitting
  norm_chi_red_ = 0.0;
  build_dataset();
  int len = int(ind_m_.size());
  const int nfree = num_free_params();
  if (len < nfree)
    throw std::runtime_error(
        "fewer data points than free parameters (simann fallback not "
        "implemented)");
  // debug: override start values from SEDFIT_INIT (name value lines), applied
  // after the initial guesses to probe convergence-path sensitivity.
  if (const char* initf = std::getenv("SEDFIT_INIT")) {
    std::ifstream in(initf);
    std::string nm;
    double val;
    while (in >> nm >> val) {
      int i = fun_.index_of(nm);
      if (i >= 0) params[i].value = val;
    }
  }
  fit_once();

  // ------- outlier removal
  if (remove_outliers > 0) {
    bool repeat;
    do {
      repeat = false;
      if (len > nfree + 1) {
        build_dataset();
        fit_once();
        std::vector<double> fullpar;
        for (const auto& p : params) fullpar.push_back(p.value);
        dvec model;
        model_vector(fullpar, model);
        int worst = 0;
        double worst_chi = 0.0;
        for (int k = 0; k < len; ++k) {
          double chi = std::fabs((model[k] - data_[k]) / err_[k]);
          if (chi > worst_chi) {
            worst_chi = chi;
            worst = k;
          }
        }
        int n_mag_left = 0;
        for (int j : ind_n_)
          if (fun_.db().entries()[j].type == "magnitude") ++n_mag_left;
        bool worst_is_mag =
            fun_.db().entries()[ind_n_[worst]].type == "magnitude";
        if ((!worst_is_mag || n_mag_left > 1) && worst_chi > remove_outliers) {
          repeat = true;
          const PhotEntry& e = phot_.entries[ind_m_[worst]];
          if (verbose)
            std::fprintf(stderr,
                         "Information: the data point %s: %s was identified "
                         "as outlier (flag=-2).\n",
                         e.system.c_str(), e.passband.c_str());
          phot_.entries[ind_m_[worst]].flag = -2;
          select_data(false);  // rebuild ind_m/ind_n/mag_ind without it
          len = int(ind_m_.size());
        }
      }
    } while (repeat);
    build_dataset();
    fit_once();
  }

  // ------- norm_chi_red (excess-noise normalization)
  double chisqr_red = std::numeric_limits<double>::quiet_NaN();
  if (len > nfree) {
    std::vector<double> fullpar;
    for (const auto& p : params) fullpar.push_back(p.value);
    dvec model;
    model_vector(fullpar, model);
    const double dof = len - nfree;
    auto red_chi2 = [&](const dvec& mod) {
      double s = 0.0;
      for (int k = 0; k < len; ++k) {
        double r = (mod[k] - data_[k]) / err_[k];
        s += r * r;
      }
      return s / dof;
    };
    chisqr_red = red_chi2(model);
    if (chisqr_red > 1.0) {
      // initial guess: average of the two limiting-case roots
      double c = chisqr_red - 1.0;
      double a = 0.0, b = 0.0;
      for (int k = 0; k < len; ++k) {
        double d2 = (model[k] - data_[k]) * (model[k] - data_[k]);
        a -= d2 / std::pow(err_[k], 4);
        b += d2;
      }
      a /= dof;
      b /= dof;
      double root1 = std::sqrt(-c / a);
      double root2 = std::sqrt(b);
      double norm_new = 0.5 * (root1 + root2);
      const double thres = 2e-6;
      int it = 0;
      double norm;
      do {
        norm = norm_new;
        norm_chi_red_ = norm;
        build_dataset();
        fit_once();
        for (size_t i = 0; i < params.size(); ++i)
          fullpar[i] = params[i].value;
        model_vector(fullpar, model);
        chisqr_red = red_chi2(model);
        double diff_sum = 0.0;
        for (int k = 0; k < len; ++k)
          diff_sum += (model[k] - data_[k]) * (model[k] - data_[k]) /
                      std::pow(err_[k], 4);
        double chisqr_red_diff = -2.0 * norm * diff_sum / dof;
        norm_new = norm - (chisqr_red - 1.0) / chisqr_red_diff;
        if (norm_new < thres) norm_new = 0.5 * norm;
        if (++it > 100)
          throw std::runtime_error(
              "more than 100 iterations to find norm_chi_red");
      } while (std::fabs(norm_new - norm) > thres);
    }
  } else if (verbose) {
    warn("not enough data points for norm_chi_red");
  }

  // ------- confidence limits
  std::vector<int> free_idx;
  for (size_t i = 0; i < params.size(); ++i)
    if (params[i].freeze == 0 && !params[i].is_fun) free_idx.push_back(int(i));
  dvec pmin(free_idx.size(), std::numeric_limits<double>::quiet_NaN());
  dvec pmax(free_idx.size(), std::numeric_limits<double>::quiet_NaN());
  bool have_conf = false;
  static const double DELTA_CHISQR[3] = {1.00, 2.71, 6.63};
  if (conf_level != -1 && len > nfree) {
    have_conf = true;
    const double delta = DELTA_CHISQR[conf_level];
    const double tol = 1e-3;
    // Confidence-limit driver, mirroring conf_loop.sl: when a parameter's
    // single-parameter search re-fits the others and finds a lower statistic
    // (or fails to bracket the best value), the improved fit is adopted and the
    // whole loop restarts (invalidating any limits found so far). This lets a
    // fit that stopped short on a flat/degenerate surface settle into the true
    // minimum before clean limits are computed. conf() leaves the improved
    // parameters in the table on an improvement and restores the best-fit
    // values otherwise.
    auto conf_all = [&]() {
      double best_stat = statistic_at_current();
      bool restart = true;
      int guard = 0;
      while (restart) {
        restart = false;
        if (++guard > 1000) break;
        for (size_t j = 0; j < free_idx.size(); ++j) {
          std::vector<double> save;
          for (const auto& p : params) save.push_back(p.value);
          double pbest_val = save[free_idx[j]];
          double lo, hi;
          conf(free_idx[j], delta, tol, lo, hi);
          double stat_after = statistic_at_current();
          bool not_bracketed = (pbest_val < lo || hi < pbest_val);
          if (stat_after < best_stat - 1e-9) {
            // adopt the improved fit and settle it, then restart the loop
            fit_once();
            best_stat = statistic_at_current();
            restart = true;
            break;
          }
          // no improvement: restore the best fit and record the limits
          for (size_t i = 0; i < params.size(); ++i)
            params[i].value = save[i];
          pmin[j] = lo;
          pmax[j] = hi;
          (void)not_bracketed;
        }
      }
    };
    if (norm_chi_red_ != 0.0) {
      double chisqr_red_old = statistic_at_current() / double(len - nfree);
      bool repeat;
      do {
        repeat = false;
        conf_all();
        double cr = statistic_at_current() / double(len - nfree);
        if (cr < 0.99 * chisqr_red_old) {
          repeat = true;
          // renormalize with Newton, starting from the previous norm
          double norm_new = norm_chi_red_;
          double norm;
          int it = 0;
          std::vector<double> fullpar;
          dvec model;
          const double dof = len - nfree;
          do {
            norm = norm_new;
            norm_chi_red_ = norm;
            build_dataset();
            fit_once();
            fullpar.clear();
            for (const auto& p : params) fullpar.push_back(p.value);
            model_vector(fullpar, model);
            double s = 0.0, diff_sum = 0.0;
            for (int k = 0; k < len; ++k) {
              double r = (model[k] - data_[k]) / err_[k];
              s += r * r;
              diff_sum += (model[k] - data_[k]) * (model[k] - data_[k]) /
                          std::pow(err_[k], 4);
            }
            double cr2 = s / dof;
            double crdiff = -2.0 * norm * diff_sum / dof;
            norm_new = norm - (cr2 - 1.0) / crdiff;
            if (norm_new < 1e-6) norm_new = 0.5 * norm;
            if (++it > 100)
              throw std::runtime_error(
                  "more than 100 iterations to find norm_chi_red");
          } while (std::fabs(norm_new - norm) > 1e-6);
          chisqr_red_old = statistic_at_current() / double(len - nfree);
        }
      } while (repeat);
    } else {
      conf_all();
    }
  } else if (conf_level != -1 && verbose) {
    warn("cannot compute confidence limits (too few data points)");
  }

  // ------- final statistics
  {
    double s = statistic_at_current();
    if (len > nfree) chisqr_red = s / double(len - nfree);
  }

  // ------- assemble results
  FitResults r;
  r.chisqr_red = chisqr_red;
  r.norm_chi_red = norm_chi_red_;
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  for (size_t i = 0; i < params.size(); ++i) {
    r.index.push_back(int(i) + 1);
    r.name.push_back(params[i].name);
    r.value.push_back(params[i].value);
    r.freeze.push_back(params[i].is_fun ? 1 : params[i].freeze);
    r.min.push_back(params[i].min);
    r.max.push_back(params[i].max);
    auto it = std::find(free_idx.begin(), free_idx.end(), int(i));
    if (have_conf && it != free_idx.end()) {
      size_t j = it - free_idx.begin();
      r.conf_min.push_back(pmin[j]);
      r.conf_max.push_back(pmax[j]);
      r.buf_below.push_back((pmin[j] - params[i].min) /
                            (params[i].max - params[i].min));
      r.buf_above.push_back((params[i].max - pmax[j]) /
                            (params[i].max - params[i].min));
      r.tex.push_back("");  // TeX string filled by the caller (texval)
    } else {
      r.conf_min.push_back(NaN);
      r.conf_max.push_back(NaN);
      r.buf_below.push_back(NaN);
      r.buf_above.push_back(NaN);
      r.tex.push_back("\\ldots");
    }
  }
  return r;
}

// ---------------------------------------------------------------------------
// Confidence-limit search, ported from ISIS fit-chisqrconf.c
// ---------------------------------------------------------------------------

namespace {
enum {
  EVAL_ERROR = -3,
  EVAL_FAILED = -2,
  EVAL_INVALID = -1,
  EVAL_OK = 0,
  EVAL_IMPROVED = 1
};
enum { USE_QUADRATIC = 0, USE_LINEAR = 1, USE_BISECT = 2 };
constexpr int MAX_BRACKET_TRIES = 30;
constexpr int MAX_CONVERGE_TRIES = 30;

struct Sample {
  double v = 0.0;
  double chisqr = 0.0;
  std::vector<double> par;  // variable params (excluding the frozen one)
};
}  // namespace

int Fitter::find_limit(int ipar, double ptest, double prange, double pbest,
                       double min_chisqr, double delta, double tol,
                       double& conf_limit) {
  auto& params = fun_.params();
  std::vector<int> varidx;
  for (size_t i = 0; i < params.size(); ++i)
    if (params[i].freeze == 0 && !params[i].is_fun && int(i) != ipar)
      varidx.push_back(int(i));
  const bool find_best = !varidx.empty();
  const double improve_tol = 0.5 * tol * delta;

  auto get_vars = [&]() {
    std::vector<double> v(varidx.size());
    for (size_t j = 0; j < varidx.size(); ++j) v[j] = params[varidx[j]].value;
    return v;
  };
  auto set_vars = [&](const std::vector<double>& v) {
    for (size_t j = 0; j < varidx.size(); ++j) params[varidx[j]].value = v[j];
  };

  // examine_fit_statistic: freeze ipar at p, re-minimize others, return chi2
  auto examine = [&](double p, double& chisqr) -> int {
    params[ipar].value = p;
    if (find_best) {
      // temporarily freeze ipar for the sub-fit
      int oldfreeze = params[ipar].freeze;
      params[ipar].freeze = 1;
      chisqr = fit_once();
      params[ipar].freeze = oldfreeze;
    } else {
      chisqr = statistic_at_current();
    }
    double dchisqr = chisqr - min_chisqr;
    if (dchisqr < 0.0 && -dchisqr > improve_tol) return EVAL_IMPROVED;
    return EVAL_OK;
  };

  conf_limit = prange;
  Sample a, b;
  double chisqr = 0.0, target_chisqr = min_chisqr + delta;

  // ---- bracket_target
  b.v = pbest;
  b.chisqr = min_chisqr;
  b.par = get_vars();
  int status = examine(ptest, chisqr);
  if (status != EVAL_OK) {
    conf_limit = ptest;
    return status;
  }
  int k = 0;
  while (chisqr < target_chisqr) {
    ++k;
    if (k > MAX_BRACKET_TRIES || ptest == prange) {
      conf_limit = prange;
      return EVAL_INVALID;
    }
    b.par = get_vars();
    ptest = 0.5 * (prange + ptest);
    status = examine(ptest, chisqr);
    if (status != EVAL_OK) {
      conf_limit = ptest;
      return status;
    }
    if (chisqr >= target_chisqr) break;
    b.v = ptest;
    b.chisqr = chisqr;
  }
  a.v = ptest;
  a.chisqr = chisqr;
  a.par = get_vars();

  // ---- converge on the target chi2
  const double accept_tol = tol * delta;
  int mode = USE_QUADRATIC;
  int count = 0;
  double test = 0.0;
  k = 0;
  std::vector<double> curpar(varidx.size());
  do {
    if (std::fabs(a.chisqr - b.chisqr) < accept_tol) {
      conf_limit = ptest;
      return EVAL_INVALID;
    }
    // interp_par: propose new ptest and starting values for the sub-fit
    switch (mode) {
      case USE_BISECT: {
        ptest = 0.5 * (a.v + b.v);
        for (size_t j = 0; j < curpar.size(); ++j)
          curpar[j] = 0.5 * (a.par[j] + b.par[j]);
        mode = USE_LINEAR;
        break;
      }
      case USE_LINEAR: {
        double rr = (target_chisqr - b.chisqr) / (a.chisqr - b.chisqr);
        ptest = rr * a.v + (1.0 - rr) * b.v;
        for (size_t j = 0; j < curpar.size(); ++j)
          curpar[j] = rr * a.par[j] + (1.0 - rr) * b.par[j];
        break;
      }
      case USE_QUADRATIC: {
        double factor =
            std::sqrt((target_chisqr - min_chisqr) / (a.chisqr - min_chisqr));
        ptest = b.v + (a.v - b.v) * factor;
        for (size_t j = 0; j < curpar.size(); ++j)
          curpar[j] = b.par[j] + (a.par[j] - b.par[j]) * factor;
        mode = USE_LINEAR;
        break;
      }
    }
    set_vars(curpar);
    status = examine(ptest, chisqr);
    if (status != EVAL_OK) {
      conf_limit = ptest;
      return status;
    }
    test = chisqr - target_chisqr;
    Sample* update = &a;
    if (test > 0.0) {
      ++count;
      update = &a;
    } else if (test < 0.0) {
      --count;
      update = &b;
    }
    update->v = ptest;
    update->chisqr = chisqr;
    update->par = get_vars();
    if (std::abs(count) > 2 || std::fabs(test) > 3) {
      count = 0;
      mode = USE_BISECT;
    }
    ++k;
  } while (std::fabs(test) > accept_tol && k < MAX_CONVERGE_TRIES);

  conf_limit = ptest;
  if (std::fabs(test) > accept_tol) {
    if (k == MAX_CONVERGE_TRIES && a.chisqr <= target_chisqr &&
        target_chisqr <= b.chisqr)
      return EVAL_OK;
    return EVAL_INVALID;
  }
  return EVAL_OK;
}

int Fitter::conf(int ipar, double delta_stat, double tol, double& lo,
                 double& hi) {
  auto& params = fun_.params();
  double curr_value = params[ipar].value;
  lo = hi = curr_value;

  // save the full initial state
  std::vector<double> initial_values;
  for (const auto& p : params) initial_values.push_back(p.value);
  int initial_freeze = params[ipar].freeze;

  params[ipar].freeze = 1;
  double pbest = curr_value;
  double min_chisqr = statistic_at_current();

  // lower limit
  double pstart = 0.5 * (pbest + params[ipar].min);
  double conf_min = curr_value, conf_max = curr_value;
  int ret = find_limit(ipar, pstart, params[ipar].min, pbest, min_chisqr,
                       delta_stat, tol, conf_min);
  if (ret == EVAL_IMPROVED) {
    // keep the improved fit: parameters already hold it
    params[ipar].freeze = initial_freeze;
    lo = hi = conf_min;
    return EVAL_IMPROVED;
  }
  if (ret == EVAL_ERROR || ret == EVAL_FAILED) {
    for (size_t i = 0; i < params.size(); ++i)
      params[i].value = initial_values[i];
    params[ipar].freeze = initial_freeze;
    lo = conf_min;
    hi = conf_max;
    return ret;
  }
  // restore, then search the upper limit
  for (size_t i = 0; i < params.size(); ++i)
    params[i].value = initial_values[i];
  pstart = pbest + 2.0 * (pbest - conf_min);
  if (pstart >= params[ipar].max) pstart = 0.5 * (pbest + params[ipar].max);
  ret = find_limit(ipar, pstart, params[ipar].max, pbest, min_chisqr,
                   delta_stat, tol, conf_max);
  if (ret == EVAL_IMPROVED) {
    params[ipar].freeze = initial_freeze;
    lo = hi = conf_max;
    return EVAL_IMPROVED;
  }
  for (size_t i = 0; i < params.size(); ++i)
    params[i].value = initial_values[i];
  params[ipar].freeze = initial_freeze;
  lo = conf_min;
  hi = conf_max;
  return 0;
}

}  // namespace sed
