// Photometric passband database: the zero-point table from
// photometric_magnitudes.sl (with optional empirical ZPO corrections), filter
// transmission curves from refdata/filter_passbands.fits.gz, analytic box
// filters, and the fixed list of colors from photometric_colors.sl.
#pragma once

#include <map>
#include <string>
#include <vector>

#include "util.hpp"

namespace sed {

struct BandEntry {
  std::string system;
  std::string passband;
  double ZP = 0.0;
  double ZP_err = 0.0;
  std::string type;  // "magnitude" or "color"
};

struct FilterCurve {
  dvec l, f;
};

class PassbandDB {
 public:
  // boxes: passband names of the form "%F_%F" for additional box filters.
  PassbandDB(const std::string& refdata_dir, bool apply_ZPO_corr,
             const std::vector<std::string>& boxes);

  // Row layout: n_magnitudes() magnitude rows followed by the 21 color rows.
  const std::vector<BandEntry>& entries() const { return entries_; }
  size_t n_magnitudes() const { return n_mag_; }

  // index of system+passband, or -1
  int find(const std::string& system, const std::string& passband) const;

  // filter curve for a magnitude row (filter_passbands equivalent)
  const FilterCurve& filter(size_t idx) const;

  // For each color row: indices of the magnitude rows it depends on.
  const std::vector<std::vector<int>>& color_deps() const { return deps_; }

  // Compute the 21 color magnitudes in place: mags has size entries().size(),
  // color slots are overwritten from the magnitude slots.
  void compute_colors(dvec& mags) const;

 private:
  void append_colors();
  // read every extension of filter_passbands.fits.gz in one pass (cfitsio
  // decompresses the whole archive per open, so one open for all curves)
  void load_archive() const;

  std::string refdata_;
  std::vector<BandEntry> entries_;
  size_t n_mag_ = 0;
  std::vector<std::vector<int>> deps_;
  mutable std::map<std::string, FilterCurve> curves_;
  mutable std::map<std::string, FilterCurve> archive_;  // EXTNAME -> curve
  mutable bool archive_loaded_ = false;
  // magnitude rows feeding compute_colors, resolved once in append_colors
  struct ColorIdx {
    int U, B, V, R, I;
    int sb, sy, sv, su, hn, hw;
    int GU, GB, GV, GG, GB1, GB2, GV1;
    int AF148, AF169;
  } cidx_{};
};

}  // namespace sed
