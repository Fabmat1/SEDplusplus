// Synthetic AB magnitudes on a fixed wavelength grid, mirroring
// photometric_magnitudes.sl: for each filter, the flux is linearly
// interpolated onto u = union(l_filter, l[ind]) and int(f*SE, l) is computed
// by trapezoidal integration; mag = -2.5*log10(int) - 2.407948242680184
// (except Stroemgren Hbeta filters, which skip the AB reference constant).
// Because the wavelength grid is fixed during a fit, the interpolation and
// integration collapse into fixed sparse weights on the flux array.
#pragma once

#include <vector>

#include "passbands.hpp"
#include "util.hpp"

namespace sed {

class SynthMag {
 public:
  // l: fixed wavelength grid; mag_indices: rows of db to prepare (magnitude
  // rows only). Filters not fully covered by l stay NaN.
  SynthMag(const PassbandDB& db, const dvec& l,
           const std::vector<int>& mag_indices);

  // Compute magnitudes for the prepared filters from flux f (same grid as l);
  // all other entries (including colors) are NaN. Output size =
  // db.entries().size().
  void magnitudes(const dvec& f, dvec& out) const;

  const PassbandDB& db() const { return db_; }

 private:
  struct FilterWeights {
    int row;              // index into db entries
    bool hbeta;           // Stroemgren Hbeta special case
    std::vector<int> idx; // flux indices
    dvec coeff;           // weights: int = sum coeff*f[idx]
  };
  const PassbandDB& db_;
  std::vector<FilterWeights> filters_;
  size_t n_entries_;
};

}  // namespace sed
