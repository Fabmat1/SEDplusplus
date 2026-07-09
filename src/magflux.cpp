#include "magflux.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "extinction.hpp"

namespace sed {

namespace {

constexpr double AB_REF = -2.407948242680184;

// Synthetic AB magnitude of (l, f) through one filter curve, an exact port of
// photometric_magnitudes.sl:755-814. NOTE the integration grid u is built from
// the UN-extended ind (l points inside the filter); ind is extended by one
// point to the left only for the flux interpolation (to avoid extrapolating f
// at u[0]). This differs from src/synthmag.cpp (the fit engine), which folds
// the prepend point into u -- that is why the magnitudes are recomputed here
// rather than reused. The S-Lang append-right branch is inert (dead code).
double synth_ab_mag(const dvec& l, const dvec& f, const FilterCurve& fc,
                    bool hbeta) {
  if (!(fc.l.front() >= l.front() && fc.l.back() <= l.back()))
    return std::numeric_limits<double>::quiet_NaN();
  std::vector<int> ind0;
  for (int i = 0; i < int(l.size()); ++i)
    if (fc.l.front() <= l[i] && l[i] <= fc.l.back()) ind0.push_back(i);
  if (ind0.empty()) return std::numeric_limits<double>::quiet_NaN();

  dvec lsub0(ind0.size());
  for (size_t k = 0; k < ind0.size(); ++k) lsub0[k] = l[ind0[k]];
  dvec u = sorted_union(fc.l, lsub0);
  dvec fu_filter = interp_linear(u, fc.l, fc.f);

  std::vector<int> ind = ind0;
  if (ind.front() > 0) ind.insert(ind.begin(), ind.front() - 1);
  dvec lsub(ind.size()), fsub(ind.size());
  for (size_t k = 0; k < ind.size(); ++k) {
    lsub[k] = l[ind[k]];
    fsub[k] = f[ind[k]];
  }
  dvec fu = interp_linear(u, lsub, fsub);
  dvec integrand(u.size());
  for (size_t k = 0; k < u.size(); ++k) integrand[k] = fu[k] * fu_filter[k];
  double integral = integrate_trapez(u, integrand);
  double mag = -2.5 * std::log10(integral);
  if (!hbeta) mag += AB_REF;
  return mag;
}

// ISIS complement(a, b) (setfuns.sl): indices into a of the elements not in b,
// returned in ascending value order of a. Ported faithfully because ISIS
// intersection() is built from it and its ordering leaks into the output row
// order (see intersection() below).
std::vector<int> complement(const std::vector<std::string>& a,
                            const std::vector<std::string>& b) {
  const int lena = int(a.size()), lenb = int(b.size());
  std::vector<int> c;
  if (lena == 0 || lenb == 0) {
    c.resize(lena);
    std::iota(c.begin(), c.end(), 0);
    return c;
  }
  std::vector<int> sia(lena), sib(lenb);
  std::iota(sia.begin(), sia.end(), 0);
  std::iota(sib.begin(), sib.end(), 0);
  std::stable_sort(sia.begin(), sia.end(),
                   [&](int i, int j) { return a[i] < a[j]; });
  std::stable_sort(sib.begin(), sib.end(),
                   [&](int i, int j) { return b[i] < b[j]; });
  c.reserve(lena);
  int ia = 0, ib = 0;
  std::string xb = b[sib[0]];
  while (ia < lena) {
    int k = sia[ia];
    const std::string& xa = a[k];
    if (xa < xb) {
      c.push_back(k);
      ia++;
      continue;
    }
    if (xb == xa) {
      ia++;
      continue;
    }
    do {
      ib++;
    } while (ib < lenb && xa > b[sib[ib]]);
    if (ib == lenb) {
      for (int t = ia; t < lena; ++t) c.push_back(sia[t]);
      break;
    }
    xb = b[sib[ib]];
    if (xa == xb) ia++;
  }
  return c;
}

// ISIS intersection(a, b): indices into a of the elements common to a and b.
std::vector<int> intersection(const std::vector<std::string>& a,
                              const std::vector<std::string>& b) {
  std::vector<int> i1 = complement(a, b);
  std::vector<std::string> ai;
  ai.reserve(i1.size());
  for (int x : i1) ai.push_back(a[x]);
  return complement(a, ai);
}

}  // namespace

MagFluxResult magnitudes_to_flux(const dvec& l, dvec f,
                                 std::vector<PhotEntry> m, const PassbandDB& db,
                                 double theta, double E_44m55, double R_55,
                                 bool errbar, bool flagged) {
  // theta + interstellar extinction scaling (xfig_photometry.sl:53-56).
  const double tscale = 0.5 * theta;
  const double tsq = tscale * tscale;
  dvec ext = extinction_factor(l, E_44m55, R_55);
  for (size_t i = 0; i < f.size(); ++i) f[i] *= tsq * ext[i];

  // struct_filter on the observed entries (xfig_photometry.sl:59-61). NaN
  // comparisons are false, matching S-Lang.
  std::vector<PhotEntry> mf;
  for (const auto& e : m) {
    if (!(-2.0 < e.magnitude && e.magnitude < 35.0)) continue;
    if (!(0.0 <= e.uncertainty && e.uncertainty < 10.0)) continue;
    if (!(e.uncertainty < 0.3 || e.flag == 0)) continue;
    mf.push_back(e);
  }

  // Synthetic magnitudes and colours of the (scaled, reddened) model, then
  // convert AB->observed system by adding the ZP (xfig_photometry.sl:69-71).
  const auto& entries = db.entries();
  dvec nmag(entries.size(), std::numeric_limits<double>::quiet_NaN());
  for (size_t r = 0; r < db.n_magnitudes(); ++r) {
    const BandEntry& be = entries[r];
    const bool hbeta = (be.system == "Stroemgren" &&
                        (be.passband == "Hbeta_wide" ||
                         be.passband == "Hbeta_narrow"));
    try {
      nmag[r] = synth_ab_mag(l, f, db.filter(r), hbeta);
    } catch (const std::exception&) {
      // ZP-table row without a filter curve (e.g. a per-catalogue
      // zero-point override like S82calib): not plottable, keep NaN
    }
  }
  db.compute_colors(nmag);     // fills the 21 colour rows
  for (size_t i = 0; i < nmag.size(); ++i) nmag[i] += entries[i].ZP;

  // keys and the measurable subset (non-NaN synthetic magnitude/colour)
  std::vector<std::string> m_keys(mf.size());
  for (size_t i = 0; i < mf.size(); ++i) m_keys[i] = mf[i].system + mf[i].passband;
  std::vector<std::string> meas_keys;
  for (size_t i = 0; i < entries.size(); ++i)
    if (!std::isnan(nmag[i]))
      meas_keys.push_back(entries[i].system + entries[i].passband);

  std::vector<int> ind_m = intersection(m_keys, meas_keys);

  MagFluxResult out;
  for (size_t i = 0; i < ind_m.size(); ++i) {
    const int mi = ind_m[i];
    const PhotEntry& obs = mf[mi];
    const int ni = db.find(obs.system, obs.passband);
    if (ni < 0) continue;  // guarded by intersection, but be safe
    const BandEntry& ne = entries[ni];
    const double unc = obs.uncertainty;
    const double diff = nmag[ni] - obs.magnitude;
    const double diff_err = std::sqrt(unc * unc + ne.ZP_err * ne.ZP_err);
    // VizieR_catalog is indexed by the loop counter in the S-Lang, not by
    // ind_m[i] (xfig_photometry.sl:99) -- replicate the exact indexing.
    const std::string vizier = mf[i].vizier_catalog;

    if (ne.type == "magnitude") {
      const FilterCurve& fc = db.filter(ni);
      const dvec& lf = fc.l;
      const dvec& ff = fc.f;
      double fmax = ff[0];
      for (double v : ff) fmax = std::max(fmax, v);
      const double thr = 0.1 * fmax;  // full width at tenth maximum
      int i0 = 0, i1 = int(ff.size()) - 1;
      for (int k = 0; k < int(ff.size()); ++k)
        if (ff[k] >= thr) {
          i0 = k;
          break;
        }
      for (int k = int(ff.size()) - 1; k >= 0; --k)
        if (ff[k] >= thr) {
          i1 = k;
          break;
        }
      const double lambda_min = lf[i0];
      const double lambda_max = lf[i1];

      // u = union(l[where(lf[0] <= l <= lf[-1])], lf)
      dvec lin;
      for (double x : l)
        if (lf.front() <= x && x <= lf.back()) lin.push_back(x);
      dvec u = sorted_union(lin, lf);
      dvec ff_u = interp_linear(u, lf, ff);          // filter onto u
      dvec fm_u = interp_linear(u, l, f);            // model onto u
      dvec num(u.size()), num_l(u.size());
      for (size_t k = 0; k < u.size(); ++k) {
        num[k] = fm_u[k] * ff_u[k];
        num_l[k] = num[k] * u[k];
      }
      const double denom_flux = integrate_trapez(u, num);
      const double flux = denom_flux / integrate_trapez(u, ff_u);
      const double lambda_eff = integrate_trapez(u, num_l) / denom_flux;

      const double flux_obs =
          std::pow(10.0, -0.4 * (obs.magnitude - nmag[ni])) * flux;
      MagFluxRow row;
      row.lambda_min = lambda_min;
      row.lambda = lambda_eff;
      row.lambda_max = lambda_max;
      row.flux = flux_obs;
      row.diff = diff;
      row.diff_err = diff_err;
      row.passband = ne.passband;
      row.system = ne.system;
      row.flag = obs.flag;
      row.vizier_catalog = vizier;
      if (errbar) {
        row.flux_min =
            std::pow(10.0, -0.4 * (obs.magnitude + diff_err - nmag[ni])) * flux;
        row.flux_max =
            std::pow(10.0, -0.4 * (obs.magnitude - diff_err - nmag[ni])) * flux;
      } else {
        row.flux_min = flux_obs;
        row.flux_max = flux_obs;
      }
      out.mag.push_back(std::move(row));
    } else {  // colour
      ColFluxRow row;
      row.diff = diff;
      row.diff_err = diff_err;
      row.passband = ne.passband;
      row.system = ne.system;
      row.flag = obs.flag;
      row.vizier_catalog = vizier;
      out.col.push_back(std::move(row));
    }
  }

  if (!flagged) {
    out.mag.erase(std::remove_if(out.mag.begin(), out.mag.end(),
                                 [](const MagFluxRow& r) { return r.flag != 0; }),
                  out.mag.end());
    out.col.erase(std::remove_if(out.col.begin(), out.col.end(),
                                 [](const ColFluxRow& r) { return r.flag != 0; }),
                  out.col.end());
  }
  return out;
}

}  // namespace sed
