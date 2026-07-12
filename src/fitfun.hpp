// The fit_theta_interext_atmparams model: parameter table with ISIS default
// values/limits and evaluation of synthetic magnitudes for a parameter
// vector. Parameter layout (per initialize_grid_fit_photometry):
//   0 logtheta, 1 E_44m55, 2 R_55,
//   per component i: c{i}_vsini, zeta, vrad, teff, logg, xi, z, HE, sur_ratio,
//   bb_teff, bb_sur_ratio, bb2_teff, bb2_sur_ratio, dummy_1..dummy_n.
#pragma once

#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "modelgrid.hpp"
#include "passbands.hpp"
#include "synthmag.hpp"
#include "util.hpp"

namespace sed {

struct Param {
  std::string name;
  double value = 0.0;
  int freeze = 1;
  double hard_min = 0.0, hard_max = 0.0;
  double min = 0.0, max = 0.0;
  double step = 0.0, relstep = 0.0;
  bool is_fun = false;  // set via set_par_fun (dummy_2..4); excluded from fit
};

class FitFunction {
 public:
  FitFunction(std::vector<std::shared_ptr<ModelGrid>> grids,
              const PassbandDB& db, int n_dummy);

  std::vector<Param>& params() { return params_; }
  const std::vector<Param>& params() const { return params_; }

  int index_of(const std::string& name) const;  // -1 if absent
  // set_par equivalents; wildcard '*' in name supported. Throws on
  // out-of-range values like ISIS (caller catches and warns).
  void set_par(const std::string& name, double value);
  void set_par(const std::string& name, double value, int freeze);
  void set_par(const std::string& name, double value, int freeze, double mn,
               double mx);
  void mark_fun(const std::string& name);  // set_par_fun: min/max -> +-DBL_MAX

  size_t n_components() const { return grids_.size(); }
  const ModelGrid& grid(size_t i) const { return *grids_[i]; }

  // Combined model spectrum (l_uni, f_uni) at the given parameter values,
  // without extinction/theta scaling. l_uni is fixed after prepare().
  // Returns a reference to the grid-cache spectrum (single component) or to
  // `scratch` (multiple components); valid until the next call.
  const dvec& l_uni() const { return l_uni_; }
  const dvec& model_flux(const std::vector<double>& par, dvec& scratch) const;

  // Prepare fixed wavelength grid and per-filter integration weights for the
  // magnitude rows in mag_ind (indices into db entries).
  void prepare(const std::vector<int>& mag_ind);

  // Full model evaluation: synthetic observed-system magnitudes for all db
  // rows (AB + theta term + colors + ZP). Requires prepare().
  void evaluate(const std::vector<double>& par, dvec& mags) const;

  const PassbandDB& db() const { return db_; }

 private:
  void check_bb(const std::vector<double>& par) const;

  std::vector<std::shared_ptr<ModelGrid>> grids_;
  const PassbandDB& db_;
  std::vector<Param> params_;
  dvec l_uni_;
  std::unique_ptr<SynthMag> synth_;
  // extinction spline components evaluated on l_uni
  dvec ext_k0_, ext_s_;
  // cached extinction factor 10^(-0.4*E*(k0+s*(R-3.02)+R)) for the last (E,R);
  // the fit varies one parameter at a time, so most evaluations reuse it
  mutable dvec ext_fac_;
  mutable double ext_E_ = -1.0, ext_R_ = -1.0;  // E>=0, so -1 = invalid

  // ---- subset fit path (single component only): the synthetic magnitudes
  // read the model flux only at the wavelengths inside the prepared filters,
  // so the fit interpolates/extincts just those points. Values are
  // bit-identical to the full path. The band integrals are additionally
  // memoized on (spectrum-cache generation, E, R) -- they are a pure function
  // of those, so replayed evaluations (e.g. conf-limit restarts re-running
  // identical fit trajectories) skip the gather entirely.
  bool sub_active_ = false;
  std::vector<int> sub_idx_;
  dvec ext_k0_sub_, ext_s_sub_;
  mutable dvec ext_fac_sub_;
  mutable double ext_Es_ = -1.0, ext_Rs_ = -1.0;
  struct MemoKey {
    std::uint64_t gen;
    double E, R;
    bool operator==(const MemoKey& o) const {
      return gen == o.gen && E == o.E && R == o.R;
    }
  };
  struct MemoHash {
    size_t operator()(const MemoKey& k) const {
      auto mix = [](std::uint64_t h, std::uint64_t v) {
        h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
      };
      std::uint64_t e, r;
      std::memcpy(&e, &k.E, 8);
      std::memcpy(&r, &k.R, 8);
      return size_t(mix(mix(k.gen, e), r));
    }
  };
  mutable std::unordered_map<MemoKey, dvec, MemoHash> memo_;

 public:
  // instrumentation (printed by main under SEDFIT_STATS)
  mutable long n_eval_ = 0, n_memo_hit_ = 0, n_ext_recompute_ = 0;
  size_t subset_size() const { return sub_idx_.size(); }

  // config fast_ext: compute the extinction factor with fast_exp10 instead of
  // pow(10,x) on the fit path (~1e-15 relative difference, not byte-exact)
  void set_fast_ext(bool on) { fast_ext_ = on; }

  // config use_mag_grid: prefer the precomputed per-band mag grid
  // (<griddir>/mag/, see maggrid.hpp) over flux interpolation + integration.
  // prepare() falls back to the flux path when the grid is absent/stale or
  // lacks a band the flux path could compute.
  void set_use_mag_grid(bool on, std::uint64_t filters_hash) {
    use_mag_ = on;
    filters_hash_ = filters_hash;
  }
  bool mag_active() const { return mag_active_; }

 private:
  bool fast_ext_ = false;
  bool use_mag_ = false;
  std::uint64_t filters_hash_ = 0;
  bool mag_active_ = false;
  struct MagRow {
    int row;     // db entry
    int band;    // index into the grid's mag band table
    bool hbeta;  // skip the AB reference constant
  };
  std::vector<MagRow> mag_rows_;
  mutable std::vector<ModelGrid::MagCorner> mag_corners_;
  // per-component interpolation of grid lambda -> l_uni (for multi-grid)
  struct Resample {
    std::vector<int> j0, j1;
    dvec t;
  };
  std::vector<Resample> resample_;
  size_t cache_capacity_;
};

}  // namespace sed
