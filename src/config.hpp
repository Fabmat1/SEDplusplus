// JSON run configuration mirroring the header block of photometry.sl
// (everything above "variable predict_mag = 0;").
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "fitter.hpp"

namespace sed {

struct Config {
  std::string star;
  std::string basename;  // derived from star if requested
  std::optional<double> ra, dec;
  std::optional<double> fix_distance, fix_distance_err;  // kpc

  ParSet par;       // name/value/freeze
  ParSet par_full;  // name/value/freeze/min/max

  std::vector<std::string> griddirectories;
  std::vector<std::string> bpaths;

  int conf_level = 0;
  // Phase-3 output-product toggles (all OFF by default; ISIS defaults differ).
  bool write_model = false;  // photometry_fit*.txt + model FITS extensions
  bool write_fits = false;   // SED_results.fits
  bool write_tex = false;    // photometry_results.tex (+ pdflatex)
  bool plot = false;         // photometry_SED.pdf (implies write_fits+write_model)
  bool save_MC = false;      // MC_c* extensions in SED_results.fits
  bool apply_ZPO_corr = true;
  double remove_outliers = 5;
  long nMC = 2000000;
  int stilism_distance_simple = 1;
  int stilism_ebmv_simple = 1;
  int stilism_ebmv_rerun = 1;
  double mass_can = 0;
  double delta_mass_can = 0.05;
  int derive_logg = 0;
  int hb_distance = 0;
  int derive_logg_c2 = 0;
  double z_c2 = -0.9;
  int derive_sr = 0;
  double sdOB_radius = 0.2;
  double R1 = 0;
  double R1_err = 0.01;

  // C++-specific
  std::string refdata;  // stellar_isisscripts refdata directory
  std::string workdir = ".";  // where photometry.dat etc. are read from
  std::string outdir;         // where outputs go (default: workdir)
  std::string plot_script;    // override for scripts/plot_sed.py
  unsigned long mc_seed = 42;

  static Config load(const std::string& json_path);
};

}  // namespace sed
