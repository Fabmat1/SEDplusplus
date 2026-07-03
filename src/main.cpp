// sedfit: C++ port of the stellar_isisscripts photometry.sl SED-fitting
// pipeline (fit engine only; photometry/astrometry are read from cached
// files produced by the original pipeline, remote querying comes later).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <set>
#include <string>

#include "config.hpp"
#include "fitfun.hpp"
#include "fitter.hpp"
#include "modelgrid.hpp"
#include "photometry_table.hpp"
#include "query.hpp"
#include "reddening.hpp"
#include "results.hpp"
#include "stellar.hpp"
#include "texval.hpp"

using namespace sed;

namespace {

struct Astrometry {
  double parallax = std::numeric_limits<double>::quiet_NaN();
  double parallax_error = std::numeric_limits<double>::quiet_NaN();
  double ra = 0, dec = 0;
  std::optional<double> ruwe, ipd_gof_harmonic_amplitude,
      visibility_periods_used;
};

std::string fmt(const char* f, double v) {
  char buf[128];
  std::snprintf(buf, sizeof buf, f, v);
  return buf;
}

// Build photometry.dat from (ra, dec): remote query + IRSA reddening, then
// write in photometric_table.sl format. Mirrors query_photometry followed by
// query_reddening in the template.
void build_photometry(double ra, double dec, const std::string& outfile) {
  PhotometryTable phot = query_photometry(ra, dec);
  try {
    phot.reddening = query_reddening(ra, dec);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "Warning: reddening query failed: %s\n", e.what());
  }
  phot.write(outfile);
}

}  // namespace

int main(int argc, char** argv) {
  auto t0 = std::chrono::steady_clock::now();
  // Standalone query mode: sedfit --query <ra> <dec> <out.dat>
  if (argc == 5 && std::string(argv[1]) == "--query") {
    build_photometry(std::atof(argv[2]), std::atof(argv[3]), argv[4]);
    std::fprintf(stderr, "- wrote %s\n", argv[4]);
    return 0;
  }
  if (argc != 2) {
    std::fprintf(stderr,
                 "usage: sedfit <config.json>\n"
                 "       sedfit --query <ra> <dec> <out.dat>\n");
    return 1;
  }
  Config cfg = Config::load(argv[1]);
  const std::string base_in = cfg.workdir + "/" + cfg.basename;
  const std::string base_out = cfg.outdir + "/" + cfg.basename;
  std::filesystem::create_directories(cfg.outdir);

  // ---------------- astrometry (from cached results header or fix_distance)
  std::optional<Astrometry> astrometry;
  bool good_astrometry = false;
  auto hdr = ResultsHeader::read(base_in + "photometry_results.txt");
  if (hdr && hdr->has("parallax") && hdr->has("parallax_error") &&
      hdr->has("RA") && hdr->has("DEC")) {
    Astrometry a;
    a.parallax = hdr->num("parallax");
    a.parallax_error = hdr->num("parallax_error");
    a.ra = hdr->num("RA");
    a.dec = hdr->num("DEC");
    if (hdr->has("ruwe")) a.ruwe = hdr->num("ruwe");
    if (hdr->has("ipd_gof_harmonic_amplitude"))
      a.ipd_gof_harmonic_amplitude = hdr->num("ipd_gof_harmonic_amplitude");
    if (hdr->has("visibility_periods_used"))
      a.visibility_periods_used = hdr->num("visibility_periods_used");
    astrometry = a;
    good_astrometry = !std::isnan(a.parallax) && a.parallax > 0;
    cfg.ra = a.ra;
    cfg.dec = a.dec;
  }
  if (!astrometry && cfg.fix_distance) {
    Astrometry a;
    a.parallax = 1.0 / *cfg.fix_distance;
    a.parallax_error = a.parallax * *cfg.fix_distance_err / *cfg.fix_distance;
    a.ruwe = 1.0;
    a.ipd_gof_harmonic_amplitude = 0.0;
    a.visibility_periods_used = 1;
    a.ra = cfg.ra.value_or(0);
    a.dec = cfg.dec.value_or(0);
    astrometry = a;
    good_astrometry = true;
  }
  if (!astrometry)
    throw std::runtime_error(
        "no cached astrometry found (" + base_in +
        "photometry_results.txt); remote Gaia queries are not implemented in "
        "this version");
  double gdist = 0.0;
  if (good_astrometry) {
    Astrometry& a = *astrometry;
    if (!(a.parallax > 0))
      a.parallax = std::numeric_limits<double>::quiet_NaN();
    if (!(a.parallax / a.parallax_error >= 1.0))
      a.parallax = std::numeric_limits<double>::quiet_NaN();
    good_astrometry = !std::isnan(a.parallax) && a.parallax > 0;
    if (good_astrometry) gdist = 1.0e3 / a.parallax;
  }

  // ---------------- photometry
  const std::string photfile = base_in + "photometry.dat";
  if (!std::filesystem::exists(photfile)) {
    // Mirror the template's stat_file(...)==NULL gate: query and cache it.
    if (!cfg.ra || !cfg.dec)
      throw std::runtime_error(
          photfile +
          " not found and no ra/dec available to query it (need cached "
          "photometry_results.txt or ra/dec in the config).");
    std::fprintf(stderr, "- querying photometry for RA=%.6f DEC=%.6f ...\n",
                 *cfg.ra, *cfg.dec);
    build_photometry(*cfg.ra, *cfg.dec, photfile);
  }
  PhotometryTable phot = PhotometryTable::read(photfile);
  if (phot.entries.empty())
    throw std::runtime_error("photometry.dat contains no entries");
  {
    bool any_unflagged = false;
    for (const auto& e : phot.entries) any_unflagged |= (e.flag == 0);
    if (!any_unflagged) throw std::runtime_error("all magnitudes are flagged");
  }
  // generic uncertainty for unrealistically small/NaN uncertainties
  for (auto& e : phot.entries)
    if (e.flag == 0 && (e.uncertainty < 1e-5 || std::isnan(e.uncertainty)))
      e.uncertainty = 0.025;

  // ---------------- first guess for E(44-55) from reddening maps
  ParSet par_full = cfg.par_full;
  {
    bool have_E = false;
    for (const auto& nm : par_full.name)
      if (nm == "E_44m55") have_E = true;
    std::optional<double> E3d = phot.reddening.meanStilism;
    std::optional<double> E2d = phot.reddening.meanSandF;
    if (!E2d) E2d = phot.reddening.meanSFD;
    if (!have_E && (E3d || E2d)) {
      bool outside_stilism = false;
      const Astrometry& a = *astrometry;
      if (a.parallax / a.parallax_error > 2 && a.parallax <= 0.5)
        outside_stilism = true;
      std::optional<double> E;
      if (!outside_stilism && E3d)
        E = *E3d / 0.91;
      else if (E2d)
        E = *E2d / 0.91;
      if (E) {
        if (*E > 2) E = 2.0;
        par_full.name.insert(par_full.name.begin(), "E_44m55");
        par_full.value.insert(par_full.value.begin(), *E);
        par_full.freeze.insert(par_full.freeze.begin(), 0);
        par_full.min.insert(par_full.min.begin(), 0.0);
        par_full.max.insert(par_full.max.begin(), 5.0);
      }
    }
  }

  // ---------------- grids and fit function
  auto griddirs = search_grid_dirs(cfg.bpaths, cfg.griddirectories, "grid.fits");
  std::vector<std::shared_ptr<ModelGrid>> grids;
  for (const auto& d : griddirs) grids.push_back(std::make_shared<ModelGrid>(d));

  std::vector<std::string> boxes;
  {
    std::set<std::string> seen;
    for (const auto& e : phot.entries)
      if (e.system == "box" && seen.insert(e.passband).second)
        boxes.push_back(e.passband);
  }
  PassbandDB db(cfg.refdata, cfg.apply_ZPO_corr, boxes);

  int ndummy = good_astrometry ? 4 : 0;
  FitFunction fun(grids, db, ndummy);
  if (good_astrometry) {
    fun.set_par("dummy_1", astrometry->parallax, -1, 0.0, 200.0);
    fun.mark_fun("dummy_2");
    fun.mark_fun("dummy_3");
    fun.mark_fun("dummy_4");
  }

  // ---------------- fit
  Fitter fitter(fun, phot);
  auto tfit0 = std::chrono::steady_clock::now();
  FitResults res = fitter.run(cfg.par, par_full, cfg.conf_level,
                              cfg.remove_outliers, true);
  auto tfit1 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "- fit done in %.1fs\n",
               std::chrono::duration<double>(tfit1 - tfit0).count());

  // debug: evaluate reduced chi^2 at externally supplied parameter values
  // (SEDFIT_EVAL=<file> of "name value" lines) using the final dataset, to
  // distinguish convergence deficiencies from genuine degeneracies. Prints
  // and exits without writing outputs.
  if (const char* evalf = std::getenv("SEDFIT_EVAL")) {
    std::ifstream in(evalf);
    std::string nm;
    double val;
    auto& params = fun.params();
    while (in >> nm >> val) {
      int i = fun.index_of(nm);
      if (i >= 0) params[i].value = val;
    }
    double stat = fitter.eval_statistic();
    int dof = fitter.n_data() - fitter.num_free_params();
    std::fprintf(stderr, "[SEDFIT_EVAL] stat=%.10g reduced=%.10g (dof=%d)\n",
                 stat, stat / dof, dof);
    std::fprintf(stderr, "[SEDFIT_EVAL] my best reduced=%.10g\n", res.chisqr_red);
    return 0;
  }

  // TeX strings for parameters with confidence limits
  for (size_t i = 0; i < res.name.size(); ++i) {
    if (!std::isnan(res.conf_min[i]))
      res.tex[i] = tex_value_pm_error(res.value[i], res.conf_min[i],
                                      res.conf_max[i], res.min[i], res.max[i],
                                      6);
  }

  // ---------------- propagate prescribed-parameter ranges into logtheta
  const int ind_logtheta = 0;
  const int ind_c1_teff = fun.index_of("c1_teff");
  if (res.tex[ind_logtheta] != "\\ldots") {
    auto& params = fun.params();
    std::vector<double> snapshot;
    for (const auto& p : params) snapshot.push_back(p.value);

    struct Asym {
      double value, minu, plus;
    };
    Asym lt{res.value[ind_logtheta],
            res.value[ind_logtheta] - res.conf_min[ind_logtheta],
            res.conf_max[ind_logtheta] - res.value[ind_logtheta]};
    Asym t1{res.value[ind_c1_teff],
            res.value[ind_c1_teff] - res.conf_min[ind_c1_teff],
            res.conf_max[ind_c1_teff] - res.value[ind_c1_teff]};

    static const char* fields[] = {"teff", "logg", "xi",     "z",
                                   "HE",   "sur_ratio", "E_44m55"};
    for (size_t comp = 1; comp <= grids.size(); ++comp) {
      for (const char* field : fields) {
        // locate the prescribed entry in par_full
        int t_full = -1;
        std::string pname;
        if (std::string(field) == "E_44m55") {
          if (comp != 1) continue;
          for (int k = int(par_full.name.size()) - 1; k >= 0; --k)
            if (par_full.name[k] == "E_44m55") {
              t_full = k;
              break;
            }
          pname = "E_44m55";
        } else {
          pname = "c" + std::to_string(comp) + "_" + field;
          for (size_t k = 0; k < par_full.name.size(); ++k)
            if (par_full.name[k] == pname) t_full = int(k);
        }
        if (t_full < 0 || par_full.freeze[t_full] != 1) continue;
        int ip = fun.index_of(pname);
        if (ip < 0) continue;
        try {
          auto apply = [&](double v) {
            params[ip].value = v;
            fitter.refit();
            double tdiff = lt.value - params[ind_logtheta].value;
            if (tdiff > 0)
              lt.minu = std::sqrt(lt.minu * lt.minu + tdiff * tdiff);
            else
              lt.plus = std::sqrt(lt.plus * lt.plus + tdiff * tdiff);
            double teff_diff = t1.value - params[ind_c1_teff].value;
            if (!std::isnan(t1.minu) && !std::isnan(t1.plus)) {
              if (teff_diff > 0)
                t1.minu = std::sqrt(t1.minu * t1.minu + teff_diff * teff_diff);
              else
                t1.plus = std::sqrt(t1.plus * t1.plus + teff_diff * teff_diff);
            } else {
              if (teff_diff > 0)
                t1.minu = std::fabs(teff_diff);
              else
                t1.plus = std::fabs(teff_diff);
            }
          };
          apply(par_full.min[t_full]);
          apply(par_full.max[t_full]);
          for (size_t i = 0; i < params.size(); ++i)
            params[i].value = snapshot[i];
        } catch (const std::exception& e) {
          std::fprintf(stderr,
                       "Warning: prescribed error propagation failed for %s\n",
                       pname.c_str());
          for (size_t i = 0; i < params.size(); ++i)
            params[i].value = snapshot[i];
        }
      }
    }
    res.conf_min[ind_logtheta] = lt.value - lt.minu;
    res.conf_max[ind_logtheta] = lt.value + lt.plus;
    res.buf_below[ind_logtheta] =
        (res.conf_min[ind_logtheta] - res.min[ind_logtheta]) /
        (res.max[ind_logtheta] - res.min[ind_logtheta]);
    res.buf_above[ind_logtheta] =
        (res.max[ind_logtheta] - res.conf_max[ind_logtheta]) /
        (res.max[ind_logtheta] - res.min[ind_logtheta]);
    res.tex[ind_logtheta] =
        tex_value_pm_error(lt.value, res.conf_min[ind_logtheta],
                           res.conf_max[ind_logtheta], res.min[ind_logtheta],
                           res.max[ind_logtheta], 6);
    res.conf_min[ind_c1_teff] = t1.value - t1.minu;
    res.conf_max[ind_c1_teff] = t1.value + t1.plus;
    res.buf_below[ind_c1_teff] =
        (res.conf_min[ind_c1_teff] - res.min[ind_c1_teff]) /
        (res.max[ind_c1_teff] - res.min[ind_c1_teff]);
    res.buf_above[ind_c1_teff] =
        (res.max[ind_c1_teff] - res.conf_max[ind_c1_teff]) /
        (res.max[ind_c1_teff] - res.min[ind_c1_teff]);
    if (res.freeze[ind_c1_teff] == 0)
      res.tex[ind_c1_teff] =
          tex_value_pm_error(t1.value, res.conf_min[ind_c1_teff],
                             res.conf_max[ind_c1_teff], res.min[ind_c1_teff],
                             res.max[ind_c1_teff], 6);
  }

  // ---------------- derived dummy parameters
  if (good_astrometry) {
    auto& params = fun.params();
    double logtheta = params[ind_logtheta].value;
    double parallax = params[fun.index_of("dummy_1")].value;
    double c1_teff = params[ind_c1_teff].value;
    double c1_logg = params[fun.index_of("c1_logg")].value;
    double R = std::pow(10.0, logtheta) / (2e-3 * parallax) *
               4.4353565926280975e7;
    double L = R * R * std::pow(c1_teff / 5772.0, 4);
    double M = std::pow(10.0, c1_logg) * R * R * 3.6469715273112305e-5;
    auto patch = [&](const char* nm, double v) {
      int i = fun.index_of(nm);
      params[i].value = v;
      res.value[i] = v;
    };
    patch("dummy_2", R);
    patch("dummy_3", L);
    patch("dummy_4", M);
    res.value[fun.index_of("dummy_1")] = parallax;
  }

  // ---------------- write photometry_results.txt
  int nmag_good = 0;
  for (const auto& e : phot.entries)
    if (e.flag == 0) ++nmag_good;
  auto grid_short = [](const std::string& dir) {
    // last two path components before the trailing slash
    std::vector<std::string> parts;
    std::string cur;
    for (char c : dir) {
      if (c == '/') {
        if (!cur.empty()) parts.push_back(cur);
        cur.clear();
      } else
        cur += c;
    }
    if (!cur.empty()) parts.push_back(cur);
    size_t n = parts.size();
    return parts[n - 2] + "/" + parts[n - 1];
  };
  std::string header;
  {
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "# RA = %.8f; DEC = %.8f; norm_chi_red = %.4f; "
                  "chisqr_reduced = %.3f;\n",
                  phot.ra.value_or(0.0), phot.dec.value_or(0.0),
                  res.norm_chi_red, res.chisqr_red);
    header = buf;
    std::snprintf(buf, sizeof buf, "# nmag_good = %d; grid = %s;\n", nmag_good,
                  grid_short(griddirs[0]).c_str());
    header += buf;
    const Astrometry& a = *astrometry;
    header += "# ";
    if (a.ruwe) header += "ruwe = " + fmt("%.3g", *a.ruwe) + "; ";
    if (a.ipd_gof_harmonic_amplitude)
      header += "ipd_gof_harmonic_amplitude = " +
                fmt("%.3g", *a.ipd_gof_harmonic_amplitude) + "; ";
    if (a.visibility_periods_used)
      header += "visibility_periods_used = " +
                fmt("%.3g", *a.visibility_periods_used) + "; ";
    header += "parallax = " + fmt("%.3g", a.parallax) + "; ";
    header += "parallax_error = " + fmt("%.3g", a.parallax_error) + "; ";
    header += "\n";
    std::string rl;
    if (phot.reddening.meanSFD)
      rl += "meanSFD = " + fmt("%.3g", *phot.reddening.meanSFD) + "; ";
    if (phot.reddening.meanSandF)
      rl += "meanSandF = " + fmt("%.3g", *phot.reddening.meanSandF) + "; ";
    if (phot.reddening.meanStilism)
      rl += "meanStilism = " + fmt("%.3g", *phot.reddening.meanStilism) + "; ";
    if (!rl.empty()) header += "# " + rl + "\n";
  }
  write_results_txt(base_out + "photometry_results.txt", header, res);

  // ---------------- MC stellar parameters per component
  double confidence = 0.68268;
  double perr_scale = 1.0;
  if (cfg.conf_level == 1) {
    confidence = 0.9;
    perr_scale = 1.645;
  } else if (cfg.conf_level == 2) {
    confidence = 0.99;
    perr_scale = 2.576;
  }
  if (good_astrometry) {
    for (size_t comp = 1; comp <= grids.size(); ++comp) {
      std::string cstr = "c" + std::to_string(comp);
      int i_teff = fun.index_of(cstr + "_teff");
      int i_logg = fun.index_of(cstr + "_logg");
      int i_sur = fun.index_of(cstr + "_sur_ratio");
      // Asymmetric MC deltas, matching the template exactly: a free parameter
      // uses its confidence interval; a frozen parameter uses its prescribed
      // [min,max] range ONLY when par_full holds the *literal* name
      // c{i}_{field} (the template's where(par_full.name==sprintf("c%d_%s",..))
      // does not match wildcard entries such as "c*_teff", so those contribute
      // zero spread).
      auto delta = [&](int idx, const char* suffix, double& dm, double& dp) {
        dm = dp = 0.0;
        if (cfg.conf_level == -1) return;
        std::string lit = cstr + "_" + suffix;
        bool prescribed = false;
        for (size_t k = 0; k < par_full.name.size(); ++k)
          if (par_full.name[k] == lit && par_full.freeze[k] == 1)
            prescribed = true;
        if (res.freeze[idx] == 0) {  // free parameter
          dp = res.conf_max[idx] - res.value[idx];
          dm = res.value[idx] - res.conf_min[idx];
        } else if (prescribed) {  // prescribed via a literal-named par_full entry
          dp = res.max[idx] - res.value[idx];
          dm = res.value[idx] - res.min[idx];
        }
      };
      StellarMCInput in{};
      in.logtheta = res.value[ind_logtheta];
      if (cfg.conf_level != -1 && res.freeze[ind_logtheta] == 0 &&
          !std::isnan(res.conf_min[ind_logtheta])) {
        in.d_logtheta_plus =
            res.conf_max[ind_logtheta] - res.value[ind_logtheta];
        in.d_logtheta_minus =
            res.value[ind_logtheta] - res.conf_min[ind_logtheta];
      }
      in.teff = res.value[i_teff];
      in.logg = res.value[i_logg];
      in.sur_ratio = res.value[i_sur];
      delta(i_teff, "teff", in.d_teff_minus, in.d_teff_plus);
      delta(i_logg, "logg", in.d_logg_minus, in.d_logg_plus);
      delta(i_sur, "sur_ratio", in.d_sur_minus, in.d_sur_plus);
      in.parallax = astrometry->parallax;
      in.parallax_error = astrometry->parallax_error * perr_scale;
      in.confidence = confidence;
      in.n_mc = cfg.nMC;
      in.seed = cfg.mc_seed + comp;

      StellarMCResult mc = stellar_mc(in);
      // mode_and_HDI at the requested confidence
      std::vector<StellarRow> rows;
      const double NaN = std::numeric_limits<double>::quiet_NaN();
      auto push = [&](const std::string& nm, const ModeHDI* s) {
        StellarRow row{nm, NaN, NaN, NaN, NaN, NaN, NaN};
        if (s) {
          row.value = s->mode;
          row.conf_min = s->HDI_lo;
          row.conf_max = s->HDI_hi;
          row.median_value = s->median;
          row.median_conf_min = s->quantile_lo;
          row.median_conf_max = s->quantile_hi;
        }
        rows.push_back(row);
      };
      push(cstr + "_R", mc.valid ? &mc.R : nullptr);
      push(cstr + "_M", mc.valid ? &mc.M : nullptr);
      push(cstr + "_L", mc.valid ? &mc.L : nullptr);
      write_stellar_txt(
          base_out + "photometry_results_stellar_" + cstr + ".txt", rows);
    }
  }


  auto t1 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "- script completed in %.1fs\n",
               std::chrono::duration<double>(t1 - t0).count());
  return 0;
}
