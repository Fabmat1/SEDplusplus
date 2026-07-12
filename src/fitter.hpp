// The photometric_fitting flow: chi-square minimization with mpfit (same
// engine and configuration as ISIS), initial guesses, outlier removal,
// excess-noise normalization (norm_chi_red), and single-parameter confidence
// limits ported from ISIS fit-chisqrconf.c (delta chi^2 = 1.00/2.71/6.63,
// tolerance 1e-3).
#pragma once

#include <limits>
#include <string>
#include <vector>

#include "fitfun.hpp"
#include "photometry_table.hpp"
#include "util.hpp"

namespace sed {

struct FitResults {
  std::vector<int> index;  // 1-based parameter indices (ISIS convention)
  std::vector<std::string> name;
  dvec value;
  std::vector<int> freeze;
  dvec min, max;
  dvec conf_min, conf_max;
  dvec buf_below, buf_above;
  std::vector<std::string> tex;
  double chisqr_red = std::numeric_limits<double>::quiet_NaN();
  double norm_chi_red = 0.0;
};

struct ParSet {
  std::vector<std::string> name;
  dvec value;
  std::vector<int> freeze;
  dvec min, max;  // only used by par_full
};

class Fitter {
 public:
  Fitter(FitFunction& fun, PhotometryTable& phot);

  // The full photometric_fitting equivalent. conf_level: -1,0,1,2.
  // remove_outliers: sigma threshold (0 disables). Updates phot flags (-2)
  // for removed outliers.
  FitResults run(const ParSet& par, const ParSet& par_full, int conf_level,
                 double remove_outliers, bool verbose);

  // Re-minimize with the current dataset (used by the prescribed-parameter
  // error propagation); returns best statistic.
  double refit();
  // Evaluate statistic at current parameter values.
  double eval_statistic();

  FitFunction& fun() { return fun_; }
  const std::vector<int>& ind_m() const { return ind_m_; }
  const std::vector<int>& ind_n() const { return ind_n_; }
  int n_data() const { return int(ind_m_.size()); }
  int num_free_params() const;

  double norm_chi_red() const { return norm_chi_red_; }

  // Cap on conf-loop restart-on-improvement cycles (and renormalize repeats).
  void set_max_conf_restarts(int n) { max_conf_restarts_ = n; }
  // Conf-limit search tolerance (ISIS conf() tol, default 1e-3 = parity).
  void set_conf_tol(double t) { conf_tol_ = t; }
  // Error estimation: false = ISIS conf-limit search (parity), true =
  // mpfit covariance errors scaled to the requested delta-chi^2 level
  // (no re-fits at all; the dominant cost of a bulk fit disappears).
  void set_covar_errors(bool on) { covar_errors_ = on; }

 private:
  void build_dataset();  // data/err arrays from phot + ZP_err (+norm)
  void select_data(bool verbose);
  double fit_once();  // mpfit over free params, returns statistic
  double statistic_at_current() ;
  void model_vector(const std::vector<double>& fullpar, dvec& out) const;
  // conf-limit search for free parameter ipar; returns (lo, hi)
  int conf(int ipar, double delta_stat, double tol, double& lo, double& hi);
  int find_limit(int ipar, double ptest, double prange, double pbest,
                 double min_chisqr, double delta, double tol,
                 double& conf_limit);

  FitFunction& fun_;
  PhotometryTable& phot_;
  std::vector<int> ind_m_;    // indices into phot entries (flag==0, computable)
  std::vector<int> ind_n_;    // indices into db entries
  std::vector<int> mag_ind_;  // magnitude rows needed
  dvec data_, err_;           // observed magnitudes and total uncertainties
  double norm_chi_red_ = 0.0;
  int max_conf_restarts_ = 1000;
  double conf_tol_ = 1e-3;
  bool covar_errors_ = false;
  // 1-sigma covariance errors of the free parameters from the most recent
  // fit_once() (mpfit xerror; free-parameter order). 0 for pegged parameters.
  dvec xerror_;
};

}  // namespace sed
