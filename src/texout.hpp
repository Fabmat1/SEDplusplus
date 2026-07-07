// photometry_results.tex writer, port of the TeX block photometry.sl:845-1304
// for this pipeline's scope: no rstilism row, no parallax_offset, no
// mass_can/R1/hb_distance branches, blackbody gates always false.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "fitter.hpp"
#include "photometry_table.hpp"
#include "stellar.hpp"

namespace sed {

struct TexInput {
  std::string star;
  int conf_level = 0;
  Reddening reddening;
  bool good_astrometry = false;
  double parallax = 0, parallax_error = 0;  // parallax_error conf-scaled once
  double ruwe = 0;
  std::optional<ModeHDI> distance;  // Gaia distance MC (only when parallax>0)
  const FitResults* res = nullptr;
  const ParSet* par_full = nullptr;
  size_t ncomp = 1;
  // per-component stellar MC (photometry.sl:1176-1222); empty or invalid
  // entries skip the R/M/L/vgrav/vesc rows (same as !good_astrometry in ISIS)
  const std::vector<StellarMCResult>* mc = nullptr;
  double norm_chi_red = 0, chisqr_red = 0;
};

void write_results_tex(const std::string& path, const TexInput& in);

// pdflatex -halt-on-error -file-line-error, run inside `outdir` so the .pdf
// lands there (photometry.sl:1303), then rm aux/log/out. Returns false (with
// a warning) when pdflatex is not in PATH.
bool run_pdflatex(const std::string& outdir, const std::string& texname);

}  // namespace sed
