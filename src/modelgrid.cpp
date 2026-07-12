#include "modelgrid.hpp"

#include <sys/stat.h>

#include <cstdio>
#include <stdexcept>

#include "fitsio_util.hpp"

namespace sed {

namespace {
bool file_exists(const std::string& p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0;
}
}  // namespace

std::vector<std::string> search_grid_dirs(
    const std::vector<std::string>& bpaths,
    std::vector<std::string> griddirectories, const std::string& sfile) {
  for (auto& dir : griddirectories) {
    if (dir.empty() || dir.back() != '/') dir += "/";
    for (auto bp : bpaths) {
      if (!bp.empty() && bp.back() != '/') bp += "/";
      if (file_exists(bp + dir + sfile)) {
        dir = bp + dir;
        break;
      }
    }
    if (!file_exists(dir + sfile))
      throw std::runtime_error("Requested grid " + dir + " does not exist!");
  }
  return griddirectories;
}

ModelGrid::ModelGrid(const std::string& location) : location_(location) {
  const std::string gridfits = location_ + "grid.fits";
  cov_.t = fits_read_double_column(gridfits + "[Grid]", "t");
  cov_.g = fits_read_double_column(gridfits + "[Grid]", "g");
  cov_.x = fits_read_double_column(gridfits + "[Grid]", "x");
  cov_.z = fits_read_double_column(gridfits + "[Grid]", "z");
  cov_.HHE = fits_read_double_column(gridfits + "[Grid]", "HHE");
  lambda_ = fits_read_double_column(location_ + "ATLAS12/lambda.fits[1]", "l");
}

void enclosing_gridpoints(double x, const dvec& arr, dvec& out) {
  out.clear();
  const size_t len = arr.size();
  if (len > 1) {
    // wherefirst(x < arr)
    size_t i = 0;
    while (i < len && !(x < arr[i])) ++i;
    if (i == 0)
      i = 1;
    else if (i == len)
      i = len - 1;
    out.push_back(arr[i - 1]);
    out.push_back(arr[i]);
  } else {
    out.push_back(arr[0]);
  }
}

const dvec& ModelGrid::node_flux(const std::string& relpath) const {
  auto it = nodes_.find(relpath);
  if (it != nodes_.end()) return it->second;
  dvec f = fits_read_double_column(location_ + relpath + ".fits[1]", "c");
  return nodes_.emplace(relpath, std::move(f)).first->second;
}

const dvec& ModelGrid::node_flux_sub(const std::string& relpath) const {
  auto it = nodes_sub_.find(relpath);
  if (it != nodes_sub_.end()) return it->second;
  const dvec& full = node_flux(relpath);
  dvec f(subset_.size());
  for (size_t k = 0; k < subset_.size(); ++k) f[k] = full[subset_[k]];
  return nodes_sub_.emplace(relpath, std::move(f)).first->second;
}

void ModelGrid::set_subset(const std::vector<int>& idx) const {
  // Keep the caches when the subset is unchanged (re-prepare after outlier
  // removal): the full-path cache also persisted across select_data, and the
  // fit-path results are only reproduced bit-exactly if cached spectra
  // survive the re-prepare the same way.
  if (idx == subset_) return;
  subset_ = idx;
  nodes_sub_.clear();
  interp_cache_sub_.clear();
  lru_sub_.clear();
}

dvec ModelGrid::interpol_recursive(const std::string& str,
                                   const char* const* fmts,
                                   const double* values,
                                   const dvec* const* gridpoints,
                                   int depth) const {
  if (depth == 0) return node_flux(str);
  // exact node match short-circuits the interpolation for this dimension
  const dvec& gp = *gridpoints[0];
  char buf[64];
  int nmatch = 0;
  for (double g : gp)
    if (values[0] == g) ++nmatch;
  if (nmatch == 1) {
    std::snprintf(buf, sizeof buf, fmts[0], values[0]);
    return interpol_recursive(str + buf, fmts + 1, values + 1, gridpoints + 1,
                              depth - 1);
  }
  // Lagrange interpolation over the (two) gridpoints, matching
  // interpol_polynomial arithmetic:  sum_j y_j * prod_{k!=j} (x-x_k)/(x_j-x_k)
  const size_t np = gp.size();
  std::vector<dvec> temp(np);
  for (size_t j = 0; j < np; ++j) {
    std::snprintf(buf, sizeof buf, fmts[0], gp[j]);
    temp[j] = interpol_recursive(str + buf, fmts + 1, values + 1,
                                 gridpoints + 1, depth - 1);
  }
  dvec out(temp[0].size(), 0.0);
  const double x = values[0];
  for (size_t j = 0; j < np; ++j) {
    double w = 1.0;
    for (size_t k = 0; k < np; ++k)
      if (k != j) w *= (x - gp[k]) / (gp[j] - gp[k]);
    const dvec& yj = temp[j];
    if (j == 0)
      for (size_t i = 0; i < out.size(); ++i) out[i] = yj[i] * w;
    else
      for (size_t i = 0; i < out.size(); ++i) out[i] += yj[i] * w;
  }
  return out;
}

dvec ModelGrid::interpol_recursive_sub(const std::string& str,
                                       const char* const* fmts,
                                       const double* values,
                                       const dvec* const* gridpoints,
                                       int depth) const {
  // identical recursion/arithmetic to interpol_recursive, gathered onto
  // subset_: every output value equals interpol_recursive(...)[subset_[i]]
  if (depth == 0) return node_flux_sub(str);
  const dvec& gp = *gridpoints[0];
  char buf[64];
  int nmatch = 0;
  for (double g : gp)
    if (values[0] == g) ++nmatch;
  if (nmatch == 1) {
    std::snprintf(buf, sizeof buf, fmts[0], values[0]);
    return interpol_recursive_sub(str + buf, fmts + 1, values + 1,
                                  gridpoints + 1, depth - 1);
  }
  const size_t np = gp.size();
  std::vector<dvec> temp(np);
  for (size_t j = 0; j < np; ++j) {
    std::snprintf(buf, sizeof buf, fmts[0], gp[j]);
    temp[j] = interpol_recursive_sub(str + buf, fmts + 1, values + 1,
                                     gridpoints + 1, depth - 1);
  }
  dvec out(temp[0].size(), 0.0);
  const double x = values[0];
  for (size_t j = 0; j < np; ++j) {
    double w = 1.0;
    for (size_t k = 0; k < np; ++k)
      if (k != j) w *= (x - gp[k]) / (gp[j] - gp[k]);
    const dvec& yj = temp[j];
    if (j == 0)
      for (size_t i = 0; i < out.size(); ++i) out[i] = yj[i] * w;
    else
      for (size_t i = 0; i < out.size(); ++i) out[i] += yj[i] * w;
  }
  return out;
}

dvec ModelGrid::interpol_fast_sub(const char* const* fmts,
                                  const double* values,
                                  const dvec* const* gridpoints) const {
  // Single-pass equivalent of interpol_recursive_sub for up to three
  // interpolating dimensions (the practical maximum: xi/z are frozen on all
  // grids in use). Per element it performs the identical sequence of rounded
  // multiplies/adds as the nested recursion, so the result is bit-identical;
  // it just avoids the per-call string building, temporaries and extra array
  // passes. Falls back to the recursion for >3 interpolating dimensions.
  char buf[64];
  // classify dims (exact-node match short-circuits like the recursion)
  bool exact[5];
  int d = 0;
  for (int i = 0; i < 5; ++i) {
    const dvec& gp = *gridpoints[i];
    int nmatch = 0;
    for (double g : gp)
      if (values[i] == g) ++nmatch;
    exact[i] = (nmatch == 1);
    if (!exact[i]) ++d;
  }
  if (d > 3)
    return interpol_recursive_sub("ATLAS12", fmts, values, gridpoints, 5);

  // free dims in hierarchy order, with Lagrange weights of their gridpoints
  int fdim[3];
  double fw0[3], fw1[3];
  const dvec* fgp[3];
  int nf = 0;
  for (int i = 0; i < 5; ++i) {
    if (exact[i]) continue;
    const dvec& gp = *gridpoints[i];
    fdim[nf] = i;
    fgp[nf] = &gp;
    if (gp.size() == 1) {
      fw0[nf] = 1.0;  // single gridpoint: weight 1 (empty Lagrange product)
      fw1[nf] = 0.0;
    } else {
      const double x = values[i];
      fw0[nf] = (x - gp[1]) / (gp[0] - gp[1]);
      fw1[nf] = (x - gp[0]) / (gp[1] - gp[0]);
    }
    ++nf;
  }

  // corner node paths: expand dimension by dimension in hierarchy order;
  // outermost free dim varies slowest (dim0-major corner order)
  const dvec* corner[8];
  int ncorner = 1;
  {
    std::string paths[8];
    paths[0] = "ATLAS12";
    for (int i = 0; i < 5; ++i) {
      if (exact[i]) {
        std::snprintf(buf, sizeof buf, fmts[i], values[i]);
        for (int c = 0; c < ncorner; ++c) paths[c] += buf;
        continue;
      }
      const dvec& gp = *gridpoints[i];
      const size_t np = gp.size();
      int newn = 0;
      std::string next[8];
      for (int c = 0; c < ncorner; ++c) {
        for (size_t k = 0; k < np; ++k) {
          std::snprintf(buf, sizeof buf, fmts[i], gp[k]);
          next[newn++] = paths[c] + buf;
        }
      }
      for (int c = 0; c < newn; ++c) paths[c] = std::move(next[c]);
      ncorner = newn;
    }
    for (int c = 0; c < ncorner; ++c) corner[c] = &node_flux_sub(paths[c]);
  }
  const size_t m = subset_.size();
  dvec out(m);
  if (d == 0) {
    out = *corner[0];
  } else if (d == 1) {
    const double* a0 = corner[0]->data();
    const double* a1 = ncorner > 1 ? corner[1]->data() : nullptr;
    const double w0 = fw0[0], w1 = fw1[0];
    if (!a1)
      for (size_t i = 0; i < m; ++i) out[i] = a0[i] * w0;
    else
      for (size_t i = 0; i < m; ++i) out[i] = a0[i] * w0 + a1[i] * w1;
  } else if (d == 2) {
    // ncorner == 4 unless an inner/outer dim has a single gridpoint
    const bool s0 = fgp[0]->size() == 1, s1 = fgp[1]->size() == 1;
    const double wa0 = fw0[0], wa1 = fw1[0], wb0 = fw0[1], wb1 = fw1[1];
    if (!s0 && !s1) {
      const double *n00 = corner[0]->data(), *n01 = corner[1]->data(),
                   *n10 = corner[2]->data(), *n11 = corner[3]->data();
      for (size_t i = 0; i < m; ++i) {
        double t0 = n00[i] * wb0 + n01[i] * wb1;
        double t1 = n10[i] * wb0 + n11[i] * wb1;
        out[i] = t0 * wa0 + t1 * wa1;
      }
    } else if (s0 && s1) {  // both single point
      const double* n0 = corner[0]->data();
      for (size_t i = 0; i < m; ++i) out[i] = (n0[i] * wb0) * wa0;
    } else if (s1) {  // inner dim single point: t = n*wb0
      const double *n0 = corner[0]->data(), *n1 = corner[1]->data();
      for (size_t i = 0; i < m; ++i)
        out[i] = (n0[i] * wb0) * wa0 + (n1[i] * wb0) * wa1;
    } else {  // outer dim single point
      const double *n0 = corner[0]->data(), *n1 = corner[1]->data();
      for (size_t i = 0; i < m; ++i)
        out[i] = (n0[i] * wb0 + n1[i] * wb1) * wa0;
    }
  } else {  // d == 3: general via per-element nesting over the corner table
    const bool s0 = fgp[0]->size() == 1, s1 = fgp[1]->size() == 1,
               s2 = fgp[2]->size() == 1;
    const int n2 = s2 ? 1 : 2, n1 = s1 ? 1 : 2;
    const double w0[2] = {fw0[0], fw1[0]}, w1[2] = {fw0[1], fw1[1]},
                 w2[2] = {fw0[2], fw1[2]};
    const int n0c = s0 ? 1 : 2;
    for (size_t i = 0; i < m; ++i) {
      double acc0 = 0.0;
      int c = 0;
      for (int a = 0; a < n0c; ++a) {
        double acc1 = 0.0;
        for (int b = 0; b < n1; ++b) {
          double acc2 = 0.0;
          for (int e = 0; e < n2; ++e) {
            const double v = (*corner[c])[i] * w2[e];
            acc2 = (e == 0) ? v : acc2 + v;
            ++c;
          }
          const double v = acc2 * w1[b];
          acc1 = (b == 0) ? v : acc1 + v;
        }
        const double v = acc1 * w0[a];
        acc0 = (a == 0) ? v : acc0 + v;
      }
      out[i] = acc0;
    }
  }
  return out;
}

const dvec& ModelGrid::interpolate_sub(double T, double G, double X, double Z,
                                       double HE, size_t cache_capacity,
                                       std::uint64_t& generation) const {
  char key[256];
  std::snprintf(key, sizeof key, "ATLAS12Z%.4fHE%.5fX%.4fG%.5fT%.2fc[1]", Z, HE,
                X, G, T);
  auto it = interp_cache_sub_.find(key);
  if (it != interp_cache_sub_.end()) {
    lru_sub_.splice(lru_sub_.end(), lru_sub_, it->second.it);
    generation = it->second.gen;
    return it->second.f;
  }
  static const char* fmts[5] = {"/Z%.2f", "/HE%.3f", "/X%.2f", "/G%.3f",
                                "/T%.0f"};
  double values[5] = {Z, HE, X, G, T};
  dvec gpz, gphe, gpx, gpg, gpt;
  enclosing_gridpoints(Z, cov_.z, gpz);
  enclosing_gridpoints(HE, cov_.HHE, gphe);
  enclosing_gridpoints(X, cov_.x, gpx);
  enclosing_gridpoints(G, cov_.g, gpg);
  enclosing_gridpoints(T, cov_.t, gpt);
  const dvec* gps[5] = {&gpz, &gphe, &gpx, &gpg, &gpt};
  dvec f = interpol_fast_sub(fmts, values, gps);

  lru_sub_.push_back(key);
  generation = ++gen_counter_;
  auto ins = interp_cache_sub_.emplace(
      key, SubEntry{std::move(f), generation, std::prev(lru_sub_.end())});
  while (lru_sub_.size() > cache_capacity) {
    interp_cache_sub_.erase(lru_sub_.front());
    lru_sub_.pop_front();
  }
  return ins.first->second.f;
}

// ---------------------------------------------------------------------------
// precomputed mag grid
// ---------------------------------------------------------------------------

bool ModelGrid::mag_available(std::uint64_t filters_hash) const {
  if (mag_state_ != 0 && mag_hash_ == filters_hash) return mag_state_ == 1;
  mag_state_ = -1;
  mag_hash_ = filters_hash;
  const std::string mdir = location_ + "mag/";
  if (!file_exists(mdir + "manifest.txt")) return false;
  int ver = 0;
  std::uint64_t h = 0;
  if (!read_mag_manifest(mdir + "manifest.txt", ver, h) ||
      ver != MAG_GRID_VERSION || h != filters_hash) {
    std::fprintf(stderr,
                 "Warning: %smanifest.txt is stale (rerun sedfit --premag); "
                 "using the flux path\n",
                 mdir.c_str());
    return false;
  }
  try {
    mag_bands_ = read_mag_bands(mdir + "bands.fits");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "Warning: %s -- using the flux path\n", e.what());
    return false;
  }
  mag_nodes_.clear();
  mag_state_ = 1;
  return true;
}

const MagNodeData& ModelGrid::mag_node(const std::string& relpath) const {
  auto it = mag_nodes_.find(relpath);
  if (it != mag_nodes_.end()) return it->second;
  const std::string p = location_ + relpath + ".fits[1]";
  MagNodeData d;
  d.I0 = fits_read_double_column(p, "I0");
  d.K = fits_read_double_column(p, "K");
  d.S = fits_read_double_column(p, "S");
  d.V00 = fits_read_double_column(p, "V00");
  d.V01 = fits_read_double_column(p, "V01");
  d.V11 = fits_read_double_column(p, "V11");
  if (d.I0.size() != mag_bands_.size())
    throw std::runtime_error("mag node " + relpath +
                             " does not match bands.fits");
  return mag_nodes_.emplace(relpath, std::move(d)).first->second;
}

void ModelGrid::mag_corners(double T, double G, double X, double Z, double HE,
                            std::vector<MagCorner>& out) const {
  static const char* fmts[5] = {"/Z%.2f", "/HE%.3f", "/X%.2f", "/G%.3f",
                                "/T%.0f"};
  const double values[5] = {Z, HE, X, G, T};
  const dvec* covs[5] = {&cov_.z, &cov_.HHE, &cov_.x, &cov_.g, &cov_.t};
  // expand corner paths and multilinear weights dimension by dimension,
  // mirroring interpol_recursive (exact-node short-circuit, Lagrange pair
  // weights, single-gridpoint weight 1)
  std::string paths[32];
  double ws[32];
  int n = 1;
  paths[0] = "mag";
  ws[0] = 1.0;
  char buf[64];
  dvec gp;
  for (int i = 0; i < 5; ++i) {
    enclosing_gridpoints(values[i], *covs[i], gp);
    int nmatch = 0;
    for (double g : gp)
      if (values[i] == g) ++nmatch;
    if (nmatch == 1 || gp.size() == 1) {
      const double v = nmatch == 1 ? values[i] : gp[0];
      std::snprintf(buf, sizeof buf, fmts[i], v);
      for (int c = 0; c < n; ++c) paths[c] += buf;
      continue;
    }
    const double x = values[i];
    const double w0 = (x - gp[1]) / (gp[0] - gp[1]);
    const double w1 = (x - gp[0]) / (gp[1] - gp[0]);
    for (int c = n - 1; c >= 0; --c) {
      std::snprintf(buf, sizeof buf, fmts[i], gp[0]);
      std::string p0 = paths[c] + buf;
      std::snprintf(buf, sizeof buf, fmts[i], gp[1]);
      paths[2 * c + 1] = paths[c] + buf;
      ws[2 * c + 1] = ws[c] * w1;
      paths[2 * c] = std::move(p0);
      ws[2 * c] = ws[c] * w0;
    }
    n *= 2;
  }
  out.clear();
  for (int c = 0; c < n; ++c) out.push_back({&mag_node(paths[c]), ws[c]});
}

const dvec& ModelGrid::interpolate(double T, double G, double X, double Z,
                                   double HE, size_t cache_capacity) const {
  char key[256];
  // same fixed-precision cache key as interpol_syn ("%s%sZ%.4fHE%.5fX%.4fG%.5fT%.2f" + "c[1]")
  std::snprintf(key, sizeof key, "ATLAS12Z%.4fHE%.5fX%.4fG%.5fT%.2fc[1]", Z, HE,
                X, G, T);
  auto it = interp_cache_.find(key);
  if (it != interp_cache_.end()) {
    // move-to-back without reallocating the list node
    lru_.splice(lru_.end(), lru_, it->second.second);
    return it->second.first;
  }
  static const char* fmts[5] = {"/Z%.2f", "/HE%.3f", "/X%.2f", "/G%.3f",
                                "/T%.0f"};
  double values[5] = {Z, HE, X, G, T};
  dvec gpz, gphe, gpx, gpg, gpt;
  enclosing_gridpoints(Z, cov_.z, gpz);
  enclosing_gridpoints(HE, cov_.HHE, gphe);
  enclosing_gridpoints(X, cov_.x, gpx);
  enclosing_gridpoints(G, cov_.g, gpg);
  enclosing_gridpoints(T, cov_.t, gpt);
  const dvec* gps[5] = {&gpz, &gphe, &gpx, &gpg, &gpt};
  dvec f = interpol_recursive("ATLAS12", fmts, values, gps, 5);

  lru_.push_back(key);
  auto ins = interp_cache_.emplace(
      key, std::make_pair(std::move(f), std::prev(lru_.end())));
  while (lru_.size() > cache_capacity) {
    interp_cache_.erase(lru_.front());
    lru_.pop_front();
  }
  return ins.first->second.first;
}

}  // namespace sed
