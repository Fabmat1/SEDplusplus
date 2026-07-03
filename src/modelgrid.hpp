// Model-grid loading and interpolation, mirroring initialize_grid_fit_photometry
// and interpol_syn (linear interpolation between the two enclosing nodes per
// dimension, recursing over Z / HE / X / G / T with exact-node short-circuit).
#pragma once

#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "util.hpp"

namespace sed {

struct GridCoverage {
  dvec t, g, x, z, HHE;
};

class ModelGrid {
 public:
  // location: grid directory ending in '/', containing grid.fits and
  // ATLAS12/lambda.fits plus the node tree ATLAS12/Z*/HE*/X*/G*/T*.fits.
  explicit ModelGrid(const std::string& location);

  const std::string& location() const { return location_; }
  const GridCoverage& coverage() const { return cov_; }
  const dvec& lambda() const { return lambda_; }

  // interpol_syn("ATLAS12", T, G, X, Z, HE, ...): calibrated flux on lambda().
  // Results are cached with the same fixed-precision keys as the S-Lang code
  // so that indistinguishable parameter values reuse the identical spectrum.
  const dvec& interpolate(double T, double G, double X, double Z, double HE,
                          size_t cache_capacity) const;

 private:
  const dvec& node_flux(const std::string& relpath) const;
  dvec interpol_recursive(const std::string& str, const char* const* fmts,
                          const double* values, const dvec* const* gridpoints,
                          int depth) const;

  std::string location_;
  GridCoverage cov_;
  dvec lambda_;

  // node-file cache (path -> flux) and interpolated-spectrum LRU cache
  mutable std::unordered_map<std::string, dvec> nodes_;
  mutable std::unordered_map<std::string, std::pair<dvec, std::list<std::string>::iterator>> interp_cache_;
  mutable std::list<std::string> lru_;
};

// search_grid_fit_photometry: resolve grid directories against base paths.
std::vector<std::string> search_grid_dirs(
    const std::vector<std::string>& bpaths,
    std::vector<std::string> griddirectories, const std::string& sfile);

// Enclosing gridpoints for x in ascending array (interpol_syn_get_min_max).
void enclosing_gridpoints(double x, const dvec& arr, dvec& out);

}  // namespace sed
