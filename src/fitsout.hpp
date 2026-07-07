// SED_results.fits writer, port of the FITS block photometry.sl:1585-1664:
// primary header (fit summary + Gaia + reddening keywords), params_fit,
// stellar_c{i} (+GRID), filters/colours/spectrum_fit (write_model), IUE
// (when an IUE spectrum backs a box magnitude), MC_c{i} (save_MC).
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "fitter.hpp"
#include "magflux.hpp"
#include "photometry_table.hpp"
#include "results.hpp"
#include "stellar.hpp"

namespace sed {

struct FitsOutInput {
  double ra = 0, dec = 0;
  double norm_chi_red = 0, chisqr_red = 0;
  int nmag_good = 0;
  std::string grid_short;  // of griddirectories[0]
  // Gaia fields (photometry.sl:1593-1600): written only when set and finite.
  // parallax_error carries the conf_level scaling, like the ISIS struct.
  std::optional<double> ruwe, ipd_gof_harmonic_amplitude,
      visibility_periods_used, parallax, parallax_error;
  bool have_reddening = false;
  Reddening reddening;  // three cards, unset fields -> undefined keywords

  const FitResults* res = nullptr;
  // per-component stellar rows + grid_short of that component's grid
  std::vector<std::vector<StellarRow>> stellar;
  std::vector<std::string> stellar_grids;

  // write_model products (null when write_model is off)
  const MagFluxResult* mout = nullptr;
  const dvec* spec_l = nullptr;
  const dvec* spec_f = nullptr;  // f*sf (extincted + theta-scaled)
  const std::vector<dvec>* spec_fcomp = nullptr;

  // IUE spectrum backing a box magnitude (photometry.sl:1487-1518,1634-1642)
  std::string iue_file;
  dvec iue_wavelength, iue_flux, iue_error;

  // save_MC: retained filtered draws per component (ragged; padded with NaN
  // in the table -- ISIS semantics unverified, no fixture has these)
  std::vector<const StellarMCResult*> mc;
};

void write_sed_results_fits(const std::string& path, const FitsOutInput& in);

// Port of the IUE-spectrum selection in the plot section
// (photometry.sl:1487-1518): for each box passband with flag==0, resolve the
// backing IUE file ("*_VI/110/inescat" -> IUE/<name>.FITS.gz relative to
// `workdir`; "Average" -> the flag==2 row with the closest magnitude), read
// it (columns wavelength/flux/sigma/quality) and keep rows with flux>0 &&
// quality==0. The LAST successfully read spectrum wins (template loop
// semantics). Fills iue_* in `out`; leaves them empty when nothing matches.
void select_iue_spectrum(const std::vector<PhotEntry>& entries,
                         const std::string& workdir, FitsOutInput& out);

}  // namespace sed
