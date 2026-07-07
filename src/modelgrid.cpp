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
    for (const auto& bp : bpaths) {
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
