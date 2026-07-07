// Model spectrum reconstruction for the output products (SED plot,
// spectrum_fit table, filters/colours flux points). Mirrors the template's
//   (l,f) = syn_spec(...; res_offset=0, res_slope=1/6, xmin=xmin-10, xmax=xmax+10)
// followed by the per-component overplot loop. Only the branch actually taken
// by photometry.sl (no macroscopic broadening: vsini=zeta=0; res_offset=0
// longslit) is implemented; other branches throw.
#pragma once

#include <vector>

#include "fitfun.hpp"
#include "util.hpp"

namespace sed {

struct ModelSpectrum {
  dvec l;                    // combined wavelength grid l_uni (Angstroem)
  dvec f;                    // combined model flux on l (unscaled, unreddened)
  std::vector<dvec> f_comp;  // per-component flux on l (only when >1 component)
};

// convolve_syn(lambda, flux, {vsini=0}; res_offset=0, res_slope) for the
// longslit / no-macro-broadening branch. Returns (olambda, cflux).
void convolve_syn_lowres(const dvec& lambda, const dvec& flux, double res_slope,
                         dvec& olambda, dvec& cflux);

// Build the plotted/output model spectrum from the fit function at its current
// parameter values, over the window [xmin, xmax] (Angstroem; the template
// passes xmin-10 .. xmax+10). res_slope = 1/6 as in photometry.sl.
ModelSpectrum build_model_spectrum(const FitFunction& fun, double xmin,
                                   double xmax, double res_slope);

}  // namespace sed
