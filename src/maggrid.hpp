// Precomputed per-band magnitude grids ("mag/" next to a grid's ATLAS12/
// flux tree), built by `sedfit --premag` and consumed by the fit's mag path.
//
// For every grid node and passband the preprocessing stores the unreddened
// band integral I0 = sum_i w_i f_i (w = SynthMag's fixed trapezoid-x-filter
// weights) together with flux-weighted moments of the extinction curve
// k(lambda, R) = k0 + s*(R - 3.02) + R across the band:
//   K   = <k0>,  S   = <s>            (p_i = w_i f_i / I0)
//   V00 = Var(k0), V01 = Cov(k0, s), V11 = Var(s)
// The reddened band integral is then reconstructed as
//   I(E, R) = I0 * exp(-beta*E*x + 0.5*(beta*E)^2*V),   beta = 0.4*ln 10,
//   x = K + S*r + R,  V = V00 + 2*V01*r + V11*r^2,      r = R - 3.02,
// a second-order cumulant expansion that is exact for narrow bands and
// accurate to O((beta*E)^3 * skewness) for broad ones (sub-mmag at E <~ 0.3,
// ~0.01 mag for the broadest bands at E ~ 1). Because band integration is
// linear in flux, interpolating I(E, R) across nodes with the flux path's
// multilinear weights is exactly equivalent to integrating the interpolated
// spectrum under this expansion.
//
// The grid depends only on the filter curves (+ analytic boxes named in the
// band table) and the node fluxes: zero-point offsets are applied at fit
// time. mag/manifest.txt records MAG_GRID_VERSION and an FNV-1a hash of the
// filter archive; a mismatch at fit time falls back to the flux path.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "util.hpp"

namespace sed {

// AB reference constant, same value as synthmag.cpp (Stroemgren Hbeta bands
// skip it).
constexpr double MAG_AB_REF = -2.407948242680184;
constexpr int MAG_GRID_VERSION = 1;

struct MagBandInfo {
  std::string system, passband;
  bool hbeta = false;
  double leff = 0.0;  // weight-averaged wavelength (Angstroem)
};

// Per-node arrays, one entry per band (bands.fits row order). NaN entries
// mark bands that could not be integrated for this node.
struct MagNodeData {
  dvec I0, K, S, V00, V01, V11;
};

// FNV-1a hash of the filter archive bytes (refdata/filter_passbands.fits[.gz])
// mixed with MAG_GRID_VERSION; cached per refdata path within the process.
std::uint64_t mag_filters_hash(const std::string& refdata);

// mag/bands.fits + mag/manifest.txt I/O (shared by --premag and the fit path)
std::vector<MagBandInfo> read_mag_bands(const std::string& bands_fits);
bool read_mag_manifest(const std::string& manifest_txt, int& version,
                       std::uint64_t& filters_hash);

// `sedfit --premag <griddir> <refdata> [--force] [box ...]` driver.
int premag_main(int argc, char** argv);

}  // namespace sed
