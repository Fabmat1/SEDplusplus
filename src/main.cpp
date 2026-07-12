// sedfit: C++ port of the stellar_isisscripts photometry.sl SED-fitting
// pipeline. Photometry, reddening, and Gaia DR3 astrometry are queried live
// when no cached files from the original pipeline are present.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>

#include "astrometry.hpp"
#include "config.hpp"
#include "simbad.hpp"
#include "extinction.hpp"
#include "fitfun.hpp"
#include "fitter.hpp"
#include "magflux.hpp"
#include "maggrid.hpp"
#include "modelgrid.hpp"
#include "photometry_table.hpp"
#include "query.hpp"
#include "reddening.hpp"
#include "results.hpp"
#include "spectrum.hpp"
#include "stellar.hpp"
#include "fitsout.hpp"
#include "texout.hpp"
#include "texval.hpp"
#include "util.hpp"

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

// Strict double parse: true only if the whole token is a number.
bool parse_double(const char* s, double& out) {
  char* end = nullptr;
  double v = std::strtod(s, &end);
  if (end == s || *end != '\0') return false;
  out = v;
  return true;
}

// Resolve an object name to coordinates via SIMBAD, or exit(1) with the same
// DataError-style message the template throws on an unresolvable object.
void resolve_or_die(const std::string& name, double& ra, double& dec) {
  std::string star = normalize_star_name(name);
  auto s = resolve_simbad(star);
  if (!s) {
    std::fprintf(stderr,
                 "Data error: Object '%s' could not be resolved and no "
                 "coordinates were provided either.\n",
                 star.c_str());
    std::exit(1);
  }
  ra = s->ra;
  dec = s->dec;
}

// Build photometry.dat from (ra, dec): remote query (IUE+MAST spectra boxes
// enabled, mirroring the template's `query_photometry(...; IUE, MAST, ...)`)
// + IRSA reddening, then write in photometric_table.sl format. Downloaded
// spectra are cached in IUE/ and MAST/ next to the output file (ISIS caches
// them in $cwd, which is where the template writes photometry.dat).
void build_photometry(double ra, double dec, const std::string& outfile) {
  QueryOptions qopt;
  std::string dir = std::filesystem::path(outfile).parent_path().string();
  qopt.spectra_dir = dir.empty() ? "." : dir;
  PhotometryTable phot = query_photometry(ra, dec, qopt);
  try {
    phot.reddening = query_reddening(ra, dec);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "Warning: reddening query failed: %s\n", e.what());
  }
  phot.write(outfile);
}

// Shared per-process model-grid pool for --multi: node-file caches persist
// across stars (exact data), the rounded-key interpolation cache is cleared
// per star so each star behaves exactly as in a single-config run.
struct GridPool {
  std::map<std::string, std::shared_ptr<ModelGrid>> grids;
  std::shared_ptr<ModelGrid> get(const std::string& dir) {
    auto it = grids.find(dir);
    if (it == grids.end())
      it = grids.emplace(dir, std::make_shared<ModelGrid>(dir)).first;
    it->second->clear_interp_cache();
    return it->second;
  }
};

void run_star(Config cfg, GridPool& pool);

}  // namespace

int main(int argc, char** argv) {
  // Standalone query mode: sedfit --query <ra> <dec> <out.dat>
  //                    or:  sedfit --query <name> <out.dat>
  if (argc >= 4 && std::string(argv[1]) == "--query") {
    double ra, dec;
    const char* out = nullptr;
    if (argc == 5 && parse_double(argv[2], ra) && parse_double(argv[3], dec)) {
      out = argv[4];
    } else if (argc == 4) {
      resolve_or_die(argv[2], ra, dec);
      out = argv[3];
    } else {
      std::fprintf(stderr,
                   "usage: sedfit --query <ra> <dec> <out.dat>\n"
                   "       sedfit --query <name> <out.dat>\n");
      return 1;
    }
    build_photometry(ra, dec, out);
    std::fprintf(stderr, "- wrote %s\n", out);
    return 0;
  }
  // Standalone astrometry mode: sedfit --astrometry <ra> <dec>
  //                        or:  sedfit --astrometry <name>
  // Prints the header-line-3 fields (%.3g, same as photometry_results.txt).
  if (argc >= 3 && std::string(argv[1]) == "--astrometry") {
    double ra, dec;
    if (argc == 4 && parse_double(argv[2], ra) && parse_double(argv[3], dec)) {
      // numeric coordinates (unchanged path)
    } else if (argc == 3) {
      resolve_or_die(argv[2], ra, dec);
    } else {
      std::fprintf(stderr,
                   "usage: sedfit --astrometry <ra> <dec>\n"
                   "       sedfit --astrometry <name>\n");
      return 1;
    }
    auto q = query_astrometry(ra, dec);
    if (!q) {
      std::fprintf(stderr, "no astrometry found\n");
      return 1;
    }
    std::printf(
        "designation = %s; good = %d;\n"
        "ruwe = %.3g; ipd_gof_harmonic_amplitude = %.3g; "
        "visibility_periods_used = %.3g; parallax = %.3g; "
        "parallax_error = %.3g; \n",
        q->designation.c_str(), q->good ? 1 : 0,
        q->ruwe.value_or(std::nan("")),
        q->ipd_gof_harmonic_amplitude.value_or(std::nan("")),
        q->visibility_periods_used.value_or(std::nan("")), q->parallax,
        q->parallax_error);
    return 0;
  }
  // Preprocessing mode: sedfit --premag <griddir> <refdata> [--force]
  // [box ...] -- build/refresh the precomputed mag grid (<griddir>/mag/).
  if (argc >= 3 && std::string(argv[1]) == "--premag")
    return premag_main(argc, argv);
  // Batch mode: sedfit --multi <file with one config path per line>.
  // Grids / filter archives / node FITS files are loaded once per process.
  if (argc == 3 && std::string(argv[1]) == "--multi") {
    std::ifstream in(argv[2]);
    if (!in) {
      std::fprintf(stderr, "cannot open config list %s\n", argv[2]);
      return 1;
    }
    GridPool pool;
    std::string line;
    int n = 0, nfail = 0;
    while (std::getline(in, line)) {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      ++n;
      std::fprintf(stderr, "=== [%d] %s\n", n, line.c_str());
      try {
        run_star(Config::load(line), pool);
      } catch (const std::exception& e) {
        ++nfail;
        std::fprintf(stderr, "ERROR [%s]: %s\n", line.c_str(), e.what());
      }
    }
    std::fprintf(stderr, "--multi: %d configs processed, %d failed\n", n,
                 nfail);
    return nfail == 0 ? 0 : 1;
  }
  if (argc != 2) {
    std::fprintf(stderr,
                 "usage: sedfit <config.json>\n"
                 "       sedfit --multi <configlist.txt>\n"
                 "       sedfit --query <ra> <dec> <out.dat>\n"
                 "       sedfit --query <name> <out.dat>\n"
                 "       sedfit --astrometry <ra> <dec>\n"
                 "       sedfit --astrometry <name>\n");
    return 1;
  }
  GridPool pool;
  run_star(Config::load(argv[1]), pool);
  return 0;
}

namespace {

void run_star(Config cfg, GridPool& pool) {
  auto t0 = std::chrono::steady_clock::now();
  const std::string base_in = cfg.workdir + "/" + cfg.basename;
  const std::string base_out = cfg.outdir + "/" + cfg.basename;
  std::filesystem::create_directories(cfg.outdir);

  // ---------------- astrometry (cached results header, fix_distance, or live
  // Gaia DR3 query -- same precedence as the template)
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
  // SIMBAD name resolution (photometry.sl:197-208): only when coordinates are
  // still unknown after the cached-header / explicit-`coordinates` paths.
  if ((!cfg.ra || !cfg.dec) && !cfg.star.empty()) {
    std::string star = normalize_star_name(cfg.star);
    std::fprintf(stderr, "- resolving '%s' via SIMBAD ...\n", star.c_str());
    auto s = resolve_simbad(star);
    if (!s)
      throw std::runtime_error(
          "Data error: Object '" + star +
          "' could not be resolved and no coordinates were provided either.");
    cfg.ra = s->ra;
    cfg.dec = s->dec;
    std::fprintf(stderr, "  -> RA=%.6f DEC=%.6f (%s)\n", s->ra, s->dec,
                 s->main_id.c_str());
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
  if (!astrometry && cfg.ra && cfg.dec) {
    // Live Gaia DR3 TAP query (port of query_astrometry as called in
    // photometry.sl:213-215). Corrections + NaN-masking happen inside.
    std::fprintf(stderr, "- querying astrometry for RA=%.6f DEC=%.6f ...\n",
                 *cfg.ra, *cfg.dec);
    auto q = query_astrometry(*cfg.ra, *cfg.dec);
    if (q) {
      Astrometry a;
      a.parallax = q->parallax;
      a.parallax_error = q->parallax_error;
      a.ruwe = q->ruwe;
      a.ipd_gof_harmonic_amplitude = q->ipd_gof_harmonic_amplitude;
      a.visibility_periods_used = q->visibility_periods_used;
      // The template keeps using the input coordinates (the header RA/DEC are
      // the requested position, not the Gaia source position).
      a.ra = *cfg.ra;
      a.dec = *cfg.dec;
      astrometry = a;
      good_astrometry = q->good;
    }
  }
  if (!astrometry)
    throw std::runtime_error(
        "no astrometry available: no cached " + base_in +
        "photometry_results.txt header, no fix_distance in the config, and "
        "the live Gaia query found nothing (or no ra/dec in the config)");
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
  for (const auto& d : griddirs) grids.push_back(pool.get(d));

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
  fun.set_fast_ext(cfg.fast_ext);
  if (cfg.use_mag_grid)
    fun.set_use_mag_grid(true, mag_filters_hash(cfg.refdata));
  if (good_astrometry) {
    fun.set_par("dummy_1", astrometry->parallax, -1, 0.0, 200.0);
    fun.mark_fun("dummy_2");
    fun.mark_fun("dummy_3");
    fun.mark_fun("dummy_4");
  }

  // ---------------- fit
  Fitter fitter(fun, phot);
  fitter.set_max_conf_restarts(cfg.max_conf_restarts);
  fitter.set_conf_tol(cfg.conf_tol);
  fitter.set_covar_errors(cfg.error_mode == "covar");
  auto tfit0 = std::chrono::steady_clock::now();
  FitResults res = fitter.run(cfg.par, par_full, cfg.conf_level,
                              cfg.remove_outliers, true);
  auto tfit1 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "- fit done in %.1fs\n",
               std::chrono::duration<double>(tfit1 - tfit0).count());
  if (std::getenv("SEDFIT_STATS"))
    std::fprintf(stderr,
                 "- stats: %ld evals (%ld memo hits, %ld ext recomputes), "
                 "subset %zu points\n",
                 fun.n_eval_, fun.n_memo_hit_, fun.n_ext_recompute_,
                 fun.subset_size());

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
    return;
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
  std::vector<StellarMCResult> mc_all;  // per component, for the tex table
  std::vector<std::vector<StellarRow>> stellar_all;  // for SED_results.fits
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
      in.retain_arrays = cfg.save_MC;  // stage 3: MC_c* FITS extensions

      StellarMCResult mc;
      if (cfg.nMC > 0) {
        mc = stellar_mc(in);
      } else {
        // nMC = 0: first-order asymmetric error propagation instead of the
        // MC (bulk fast path). Same central values as the derived dummies;
        // the deltas combine in quadrature in ln-space, correlations between
        // fit parameters are ignored (the MC ignores them too -- it samples
        // each delta independently). vgrav/vesc stay invalid (not written).
        const double plx = in.parallax, splx = in.parallax_error;
        mc.valid = std::isfinite(in.logtheta) && plx > 0 &&
                   std::isfinite(splx) && in.sur_ratio > 0;
        if (mc.valid) {
          const double R = std::pow(10.0, in.logtheta) / (2e-3 * plx) *
                           4.4353565926280975e7 * std::sqrt(in.sur_ratio);
          const double L = R * R * std::pow(in.teff / 5772.0, 4);
          const double M =
              std::pow(10.0, in.logg) * R * R * 3.6469715273112305e-5;
          auto q = [](double a, double b, double c) {
            return std::sqrt(a * a + b * b + c * c);
          };
          const double sp = splx / plx;
          const double lnR_m = q(M_LN10 * in.d_logtheta_minus,
                                 0.5 * in.d_sur_minus / in.sur_ratio, sp);
          const double lnR_p = q(M_LN10 * in.d_logtheta_plus,
                                 0.5 * in.d_sur_plus / in.sur_ratio, sp);
          const double lnL_m = q(2 * lnR_m, 4 * in.d_teff_minus / in.teff, 0);
          const double lnL_p = q(2 * lnR_p, 4 * in.d_teff_plus / in.teff, 0);
          const double lnM_m = q(2 * lnR_m, M_LN10 * in.d_logg_minus, 0);
          const double lnM_p = q(2 * lnR_p, M_LN10 * in.d_logg_plus, 0);
          auto fill = [](ModeHDI& s, double v, double rm, double rp) {
            s.mode = s.median = v;
            s.HDI_lo = s.quantile_lo = std::max(0.0, v * (1.0 - rm));
            s.HDI_hi = s.quantile_hi = v * (1.0 + rp);
            s.p_s = s.p_a = 0.0;
          };
          fill(mc.R, R, lnR_m, lnR_p);
          fill(mc.L, L, lnL_m, lnL_p);
          fill(mc.M, M, lnM_m, lnM_p);
        }
      }
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
      stellar_all.push_back(std::move(rows));
      mc_all.push_back(std::move(mc));
    }
  }

  // ---------------- photometry_results.tex (photometry.sl:845-1304)
  if (cfg.write_tex) {
    TexInput ti;
    ti.star = cfg.star;
    ti.conf_level = cfg.conf_level;
    ti.reddening = phot.reddening;
    ti.good_astrometry = good_astrometry;
    ti.parallax = astrometry->parallax;
    // conf_level 1/2 scale the parallax error ONCE (photometry.sl:858-866);
    // the stellar MC above already used the same scaled value.
    ti.parallax_error = astrometry->parallax_error * perr_scale;
    ti.ruwe = astrometry->ruwe.value_or(
        std::numeric_limits<double>::quiet_NaN());
    if (good_astrometry && astrometry->parallax > 0)
      ti.distance = gaia_distance_mc(ti.parallax, ti.parallax_error,
                                     confidence, cfg.nMC, cfg.mc_seed);
    ti.res = &res;
    ti.par_full = &par_full;
    ti.ncomp = grids.size();
    ti.mc = &mc_all;
    ti.norm_chi_red = res.norm_chi_red;
    ti.chisqr_red = res.chisqr_red;
    write_results_tex(base_out + "photometry_results.tex", ti);
    run_pdflatex(cfg.outdir, cfg.basename + "photometry_results.tex");
  }

  // ---------------- model spectrum + write_model text products
  // (photometry.sl:1306-1476). Gated on write_model (plot implies write_model).
  std::optional<MagFluxResult> mout_opt;  // kept for SED_results.fits
  dvec spec_l, spec_f;
  std::vector<dvec> spec_fcomp;
  if (cfg.write_model || cfg.plot) {
    auto& params = fun.params();
    const double res_slope = 1.0 / 6.0;  // FWHM = 6 A (low-res IUE-like)
    const double xmin = 1100, xmax = 59000;
    const double logtheta = params[fun.index_of("logtheta")].value;
    const double E_44m55 = params[fun.index_of("E_44m55")].value;
    const double R_55 = params[fun.index_of("R_55")].value;
    const double theta = std::pow(10.0, logtheta);

    // (l,f) = syn_spec(...; res_offset=0, res_slope=1/6, xmin-10..xmax+10)
    ModelSpectrum ms =
        build_model_spectrum(fun, xmin - 10, xmax + 10, res_slope);

    // Uncertainty inflation for the figure/tables only (photometry.sl:1344).
    std::vector<PhotEntry> entries = phot.entries;
    if (res.norm_chi_red != 0)
      for (auto& e : entries)
        e.uncertainty = std::sqrt(e.uncertainty * e.uncertainty +
                                  res.norm_chi_red * res.norm_chi_red);

    MagFluxResult mout = magnitudes_to_flux(ms.l, ms.f, entries, db, theta,
                                            E_44m55, R_55, true, true);
    write_mag_txt(base_out + "photometry_fit_mag.txt", mout.mag);
    if (!mout.col.empty())
      write_col_txt(base_out + "photometry_fit_col.txt", mout.col);

    // sout = {l, f*sf [, f_c{i}*sf]}; sf = extinction * (theta/2)^2.
    dvec sf = extinction_factor(ms.l, E_44m55, R_55);
    const double s2 = 0.5 * theta * (0.5 * theta);
    dvec fscaled(ms.l.size());
    for (size_t i = 0; i < ms.l.size(); ++i) fscaled[i] = ms.f[i] * sf[i] * s2;
    std::vector<dvec> fcomp;
    for (const auto& fc : ms.f_comp) {
      dvec c(fc.size());
      for (size_t i = 0; i < fc.size(); ++i) c[i] = fc[i] * sf[i] * s2;
      fcomp.push_back(std::move(c));
    }
    write_spectrum_txt(base_out + "photometry_fit.txt", ms.l, fscaled, fcomp);
    mout_opt = std::move(mout);
    spec_l = std::move(ms.l);
    spec_f = std::move(fscaled);
    spec_fcomp = std::move(fcomp);
  }

  // ---------------- SED_results.fits (photometry.sl:1585-1664)
  const std::string fits_path = base_out + "SED_results.fits";
  if (cfg.write_fits) {
    FitsOutInput fo;
    fo.ra = phot.ra.value_or(0.0);
    fo.dec = phot.dec.value_or(0.0);
    fo.norm_chi_red = res.norm_chi_red;
    fo.chisqr_red = res.chisqr_red;
    fo.nmag_good = nmag_good;
    fo.grid_short = grid_short(griddirs[0]);
    {
      const Astrometry& a = *astrometry;
      fo.ruwe = a.ruwe;
      fo.ipd_gof_harmonic_amplitude = a.ipd_gof_harmonic_amplitude;
      fo.visibility_periods_used = a.visibility_periods_used;
      if (!std::isnan(a.parallax)) fo.parallax = a.parallax;
      // the ISIS struct carries the conf_level-scaled error by this point
      if (!std::isnan(a.parallax_error))
        fo.parallax_error = a.parallax_error * perr_scale;
    }
    fo.have_reddening = phot.reddening.meanSFD || phot.reddening.meanSandF ||
                        phot.reddening.meanStilism;
    fo.reddening = phot.reddening;
    fo.res = &res;
    fo.stellar = stellar_all;
    for (const auto& d : griddirs) fo.stellar_grids.push_back(grid_short(d));
    if (mout_opt) {
      fo.mout = &*mout_opt;
      fo.spec_l = &spec_l;
      fo.spec_f = &spec_f;
      fo.spec_fcomp = &spec_fcomp;
    }
    select_iue_spectrum(phot.entries, cfg.workdir, fo);
    if (cfg.save_MC)
      for (const auto& mc : mc_all) fo.mc.push_back(&mc);
    write_sed_results_fits(fits_path, fo);
  }

  // ---------------- SED plot via scripts/plot_sed.py (replaces the xfig
  // rendering, photometry.sl:1380-1525)
  if (cfg.plot) {
    std::string script = cfg.plot_script;
    if (script.empty()) {
      std::error_code ec;
      auto exe = std::filesystem::canonical("/proc/self/exe", ec);
      if (!ec) {
        auto cand = exe.parent_path().parent_path() / "scripts" / "plot_sed.py";
        if (std::filesystem::exists(cand)) script = cand.string();
      }
    }
    if (script.empty() || !std::filesystem::exists(script)) {
      std::fprintf(stderr,
                   "Warning: plot_sed.py not found (set \"plot_script\" in "
                   "the config); skipping the SED plot\n");
    } else {
      std::string cmd = "python3 '" + script + "' '" + fits_path + "' '" +
                        base_out + "photometry_SED.pdf'";
      if (std::system(cmd.c_str()) != 0)
        std::fprintf(stderr, "Warning: SED plot failed: %s\n", cmd.c_str());
    }
  }

  auto t1 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "- script completed in %.1fs\n",
               std::chrono::duration<double>(t1 - t0).count());
}

}  // namespace
