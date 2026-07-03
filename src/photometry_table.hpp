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
