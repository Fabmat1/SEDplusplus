// Model-grid loading and interpolation, mirroring initialize_grid_fit_photometry
// and interpol_syn (linear interpolation between the two enclosing nodes per
// dimension, recursing over Z / HE / X / G / T with exact-node short-circuit).
#pragma once

#include <cstdint>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "maggrid.hpp"
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

  // Restrict the fit-path interpolation to a subset of lambda() indices
  // (ascending). Clears the subset caches; called once per star from
  // FitFunction::prepare().
  void set_subset(const std::vector<int>& idx) const;

  // Like interpolate(), but computes flux only at the set_subset() indices.
  // Uses the same rounded keys, capacity and LRU policy as interpolate(), so
  // the hit/miss history of the fit path is unchanged, and the values equal
  // interpolate(...)[idx] point by point. `generation` identifies the cached
  // entry: it changes whenever an entry is (re)computed, so equal generations
  // guarantee bit-identical spectra (used as a memoization key upstream).
  const dvec& interpolate_sub(double T, double G, double X, double Z, double HE,
                              size_t cache_capacity,
                              std::uint64_t& generation) const;

  // ---- precomputed mag grid (mag/ next to ATLAS12/, see maggrid.hpp) ----
  // Load bands + manifest once; false when mag/ is absent or its manifest
  // does not match (version, filter-archive hash). Warns once on mismatch.
  bool mag_available(std::uint64_t filters_hash) const;
  const std::vector<MagBandInfo>& mag_bands() const { return mag_bands_; }
  struct MagCorner {
    const MagNodeData* node;
    double w;
  };
  // Corner nodes + multilinear weights, replicating interpolate()'s
  // enclosing-gridpoint / exact-node semantics (so the mag path visits the
  // same nodes the flux path would).
  void mag_corners(double T, double G, double X, double Z, double HE,
                   std::vector<MagCorner>& out) const;

  // Drop the interpolated-spectrum caches (NOT the node-file cache). Called
  // between stars in --multi mode: the interpolation cache uses rounded keys
  // and bounded capacity, so its content is part of the per-star behaviour,
  // while node fluxes are exact and safe to share.
  void clear_interp_cache() const {
    interp_cache_.clear();
    lru_.clear();
    interp_cache_sub_.clear();
    lru_sub_.clear();
    nodes_sub_.clear();
  }

 private:
  const dvec& node_flux(const std::string& relpath) const;
  const dvec& node_flux_sub(const std::string& relpath) const;
  dvec interpol_recursive(const std::string& str, const char* const* fmts,
                          const double* values, const dvec* const* gridpoints,
                          int depth) const;
  dvec interpol_recursive_sub(const std::string& str, const char* const* fmts,
                              const double* values,
                              const dvec* const* gridpoints, int depth) const;
  dvec interpol_fast_sub(const char* const* fmts, const double* values,
                         const dvec* const* gridpoints) const;

  std::string location_;
  GridCoverage cov_;
  dvec lambda_;

  // node-file cache (path -> flux) and interpolated-spectrum LRU cache
  mutable std::unordered_map<std::string, dvec> nodes_;
  mutable std::unordered_map<std::string, std::pair<dvec, std::list<std::string>::iterator>> interp_cache_;
  mutable std::list<std::string> lru_;

  // subset fit path: lambda indices, node fluxes gathered onto the subset,
  // and the subset interpolation LRU cache (value also carries the entry's
  // generation id)
  mutable std::vector<int> subset_;
  mutable std::unordered_map<std::string, dvec> nodes_sub_;
  struct SubEntry {
    dvec f;
    std::uint64_t gen;
    std::list<std::string>::iterator it;
  };
  mutable std::unordered_map<std::string, SubEntry> interp_cache_sub_;
  mutable std::list<std::string> lru_sub_;
  mutable std::uint64_t gen_counter_ = 0;

  // mag grid: state (0 unchecked, 1 ok, -1 unavailable), band table, and the
  // node cache (exact data, persists across stars like nodes_)
  const MagNodeData& mag_node(const std::string& relpath) const;
  mutable int mag_state_ = 0;
  mutable std::uint64_t mag_hash_ = 0;
  mutable std::vector<MagBandInfo> mag_bands_;
  mutable std::unordered_map<std::string, MagNodeData> mag_nodes_;
};

// search_grid_fit_photometry: resolve grid directories against base paths.
std::vector<std::string> search_grid_dirs(
    const std::vector<std::string>& bpaths,
    std::vector<std::string> griddirectories, const std::string& sfile);

// Enclosing gridpoints for x in ascending array (interpol_syn_get_min_max).
void enclosing_gridpoints(double x, const dvec& arr, dvec& out);

}  // namespace sed
