// photometry.dat reading/writing (photometric_table.sl format) and the
// "# key = value;" header of photometry_results.txt (astrometry cache).
#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "util.hpp"

namespace sed {

struct PhotEntry {
  // Bitmask, 0 = pristine (only flag==0 rows are fitted). Bits are
  // independent and can be combined:
  //   1  survey_bad    failed the survey's own quality flags / SNR / saturation
  //   2  superseded    a better measurement of this band exists (redundant,
  //                    would double-count in an SED)
  //   4  internal_bad  flagged unreliable by our own testing/curation
  //                    (e.g. all APASS, blended identifications)
  //   8  sed_outlier   inconsistent with the object's other bands (set by the
  //                    post-SED outlier pass, see Fitter remove_outliers)
  // Note: the legacy VizieR remote-query path (query.cpp/spectra.cpp) still
  // uses the historic S-Lang flag values (-1, -3, ..., 2) for parity with the
  // original tool; the bulk pipeline never mixes with that path.
  int flag = 0;
  std::string system;
  std::string passband;
  double magnitude = 0.0;
  double uncertainty = 0.0;
  std::string type;  // "magnitude" or "color"
  double angu_dist_arcsec = 0.0;
  std::string vizier_catalog;
};

struct Reddening {
  std::optional<double> meanSFD, stdSFD;
  std::optional<double> meanSandF, stdSandF;
  std::optional<double> meanStilism, stdStilism;
};

struct PhotometryTable {
  std::optional<double> ra, dec;
  Reddening reddening;
  std::vector<PhotEntry> entries;

  static PhotometryTable read(const std::string& filename);

  // Write in photometric_table.sl format: "# RA/DEC" + reddening header lines,
  // then the print_struct table (right-aligned, "%S"-formatted numbers).
  void write(const std::string& filename) const;
};

// Header of photometry_results.txt: "# key = value; key = value;" lines.
// String values are kept as strings; numeric values parsed to double.
struct ResultsHeader {
  std::map<std::string, std::string> raw;
  bool has(const std::string& k) const { return raw.count(k) > 0; }
  double num(const std::string& k) const;

  static std::optional<ResultsHeader> read(const std::string& filename);
};

}  // namespace sed
