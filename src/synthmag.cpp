#include "synthmag.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sed {

namespace {
constexpr double AB_REF = -2.407948242680184;  // 2.5*log10(c*1e8*10^(-48.6/2.5))

// bracketing indices + weight for ISIS interpolate_d on subgrid xs
void interp_nodes(double x, const dvec& xs, int& j0, int& j1, double& t) {
  const int n = int(xs.size());
  int n0;
  if (x < xs[0])
    n0 = 0;
  else if (x >= xs[n - 1])
    n0 = n - 1;
  else {
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      int mid = (lo + hi) / 2;
      if (xs[mid] <= x)
        lo = mid;
      else
        hi = mid;
    }
    n0 = lo;
  }
  if (x == xs[n0]) {
    j0 = j1 = n0;
    t = 0.0;
    return;
  }
  int n1 = n0 + 1;
  if (n1 == n) {
    if (n0 == 0) {
      j0 = j1 = n0;
      t = 0.0;
      return;
    }
    n1 = n0 - 1;
  }
  if (xs[n1] == xs[n0]) {
    j0 = j1 = n0;
    t = 0.0;
    return;
  }
  j0 = n0;
  j1 = n1;
  t = (x - xs[n0]) / (xs[n1] - xs[n0]);
}
}  // namespace

SynthMag::SynthMag(const PassbandDB& db, const dvec& l,
                   const std::vector<int>& mag_indices)
    : db_(db), n_entries_(db.entries().size()) {
  for (int row : mag_indices) {
    const FilterCurve* fcp;
    try {
      fcp = &db_.filter(row);
    } catch (const std::exception&) {
      continue;  // no filter curve for this ZP row -> synthetic mag stays NaN
    }
    const FilterCurve& fc = *fcp;
    if (!(fc.l.front() >= l.front() && fc.l.back() <= l.back())) continue;

    // ind = where(l_filter[0] <= l <= l_filter[-1]), then prepend one point
    // to the left if possible (the S-Lang code's append branch is inert).
    std::vector<int> ind;
    for (int i = 0; i < int(l.size()); ++i)
      if (fc.l.front() <= l[i] && l[i] <= fc.l.back()) ind.push_back(i);
    if (ind.empty()) continue;
    if (ind.front() > 0) ind.insert(ind.begin(), ind.front() - 1);

    dvec u;
    {
      dvec lsub(ind.size());
      for (size_t k = 0; k < ind.size(); ++k) lsub[k] = l[ind[k]];
      u = sorted_union(fc.l, lsub);
    }
    dvec fu_filter = interp_linear(u, fc.l, fc.f);

    // subgrid used for flux interpolation
    dvec lsub(ind.size());
    for (size_t k = 0; k < ind.size(); ++k) lsub[k] = l[ind[k]];

    // trapezoid weight of each u node times the filter value there
    const size_t nu = u.size();
    FilterWeights fw;
    fw.row = row;
    const BandEntry& be = db_.entries()[row];
    fw.hbeta = (be.system == "Stroemgren" &&
                (be.passband == "Hbeta_wide" || be.passband == "Hbeta_narrow"));

    // accumulate weights on the global flux indices
    std::vector<double> acc(l.size(), 0.0);
    for (size_t k = 0; k < nu; ++k) {
      double wtrap = 0.0;
      if (k > 0) wtrap += 0.5 * (u[k] - u[k - 1]);
      if (k + 1 < nu) wtrap += 0.5 * (u[k + 1] - u[k]);
      double w = wtrap * fu_filter[k];
      if (w == 0.0) continue;
      int j0, j1;
      double t;
      interp_nodes(u[k], lsub, j0, j1, t);
      acc[ind[j0]] += w * (1.0 - t);
      acc[ind[j1]] += w * t;
    }
    for (size_t i = 0; i < acc.size(); ++i) {
      if (acc[i] != 0.0) {
        fw.idx.push_back(int(i));
        fw.coeff.push_back(acc[i]);
      }
    }
    filters_.push_back(std::move(fw));
  }
}

std::vector<int> SynthMag::used_indices() const {
  std::vector<int> u;
  for (const auto& fw : filters_) u.insert(u.end(), fw.idx.begin(), fw.idx.end());
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return u;
}

void SynthMag::remap_to_subset(const std::vector<int>& sub) {
  for (auto& fw : filters_) {
    fw.idx_sub.resize(fw.idx.size());
    for (size_t k = 0; k < fw.idx.size(); ++k) {
      auto it = std::lower_bound(sub.begin(), sub.end(), fw.idx[k]);
      if (it == sub.end() || *it != fw.idx[k])
        throw std::runtime_error("SynthMag subset does not cover filter index");
      fw.idx_sub[k] = int(it - sub.begin());
    }
  }
}

void SynthMag::integrals_sub(const dvec& f_sub, const dvec& ext_sub,
                             dvec& out) const {
  out.resize(filters_.size());
  for (size_t j = 0; j < filters_.size(); ++j) {
    const auto& fw = filters_[j];
    double integral = 0.0;
    for (size_t k = 0; k < fw.idx_sub.size(); ++k) {
      const int i = fw.idx_sub[k];
      integral += fw.coeff[k] * (f_sub[i] * ext_sub[i]);
    }
    out[j] = integral;
  }
}

void SynthMag::mags_from_integrals(const dvec& integrals, dvec& out) const {
  out.assign(n_entries_, std::numeric_limits<double>::quiet_NaN());
  for (size_t j = 0; j < filters_.size(); ++j) {
    const auto& fw = filters_[j];
    double mag = -2.5 * std::log10(integrals[j]);
    if (!fw.hbeta) mag += AB_REF;
    out[fw.row] = mag;
  }
}

void SynthMag::magnitudes(const dvec& f, dvec& out) const {
  out.assign(n_entries_, std::numeric_limits<double>::quiet_NaN());
  for (const auto& fw : filters_) {
    double integral = 0.0;
    for (size_t k = 0; k < fw.idx.size(); ++k)
      integral += fw.coeff[k] * f[fw.idx[k]];
    double mag = -2.5 * std::log10(integral);
    if (!fw.hbeta) mag += AB_REF;
    out[fw.row] = mag;
  }
}

}  // namespace sed
