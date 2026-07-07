#include "fitsout.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <set>

#include "fitsio_util.hpp"

namespace sed {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

FitsCol dcol(std::string name, dvec v) {
  FitsCol c;
  c.name = std::move(name);
  c.kind = 'D';
  c.d = std::move(v);
  return c;
}

FitsCol icol(std::string name, std::vector<int> v) {
  FitsCol c;
  c.name = std::move(name);
  c.kind = 'J';
  c.i = std::move(v);
  return c;
}

FitsCol scol(std::string name, std::vector<std::string> v) {
  FitsCol c;
  c.name = std::move(name);
  c.kind = 'A';
  c.s = std::move(v);
  return c;
}

}  // namespace

void write_sed_results_fits(const std::string& path, const FitsOutInput& in) {
  FitsWriter w(path);
  w.primary_hdu();
  // header struct photometry.sl:1588-1610; field order = keyword order
  w.key_double("RA", in.ra);
  w.key_double("DEC", in.dec);
  // norm_chi_red stays S-Lang Integer_Type 0 when no renormalization ran
  if (in.norm_chi_red == 0.0)
    w.key_int("NORM_CHI_RED", 0);
  else
    w.key_double("NORM_CHI_RED", in.norm_chi_red);
  w.key_double("CHISQR_RED", in.chisqr_red);
  w.key_int("NMAG_GOOD", in.nmag_good);
  w.key_string("GRID_SHORT", in.grid_short);
  auto gaia = [&](const char* name, const std::optional<double>& v) {
    if (v && !std::isnan(*v)) w.key_double(name, *v);
  };
  gaia("RUWE", in.ruwe);
  gaia("IPD_GOF_HARMONIC_AMPLITUDE", in.ipd_gof_harmonic_amplitude);
  // vpu is Integer_Type in the ISIS astrometry struct (fixture: "= 19")
  if (in.visibility_periods_used && !std::isnan(*in.visibility_periods_used)) {
    double v = *in.visibility_periods_used;
    if (v == std::floor(v))
      w.key_int("VISIBILITY_PERIODS_USED", int(v));
    else
      w.key_double("VISIBILITY_PERIODS_USED", v);
  }
  gaia("PARALLAX", in.parallax);
  gaia("PARALLAX_ERROR", in.parallax_error);
  if (in.have_reddening) {
    // all three hfields_redd cards exist; NULL fields become undefined
    // keywords (photometry.sl:1604-1609 sets only existing fields)
    auto redd = [&](const char* name, const std::optional<double>& v) {
      if (v)
        w.key_double(name, *v);
      else
        w.key_undefined(name);
    };
    redd("MEANSFD", in.reddening.meanSFD);
    redd("MEANSANDF", in.reddening.meanSandF);
    redd("MEANSTILISM", in.reddening.meanStilism);
  }

  const FitResults& r = *in.res;
  w.binary_table(
      "params_fit",
      {icol("index", r.index), scol("name", r.name), dcol("value", r.value),
       icol("freeze", r.freeze), dcol("min", r.min), dcol("max", r.max),
       dcol("conf_min", r.conf_min), dcol("conf_max", r.conf_max),
       dcol("buf_below", r.buf_below), dcol("buf_above", r.buf_above),
       scol("tex", r.tex)});

  for (size_t i = 0; i < in.stellar.size(); ++i) {
    const auto& rows = in.stellar[i];
    std::vector<std::string> name;
    dvec value, cmin, cmax, mval, mcmin, mcmax;
    for (const auto& row : rows) {
      name.push_back(row.name);
      value.push_back(row.value);
      cmin.push_back(row.conf_min);
      cmax.push_back(row.conf_max);
      mval.push_back(row.median_value);
      mcmin.push_back(row.median_conf_min);
      mcmax.push_back(row.median_conf_max);
    }
    w.binary_table("stellar_c" + std::to_string(i + 1),
                   {scol("name", name), dcol("value", value),
                    dcol("conf_min", cmin), dcol("conf_max", cmax),
                    dcol("median_value", mval), dcol("median_conf_min", mcmin),
                    dcol("median_conf_max", mcmax)});
    w.key_string("GRID", in.stellar_grids[i]);
  }

  if (in.mout) {
    const auto& mag = in.mout->mag;
    dvec lmin, lam, lmax, fmin, flx, fmax, dif, derr;
    std::vector<std::string> pb, sys, cat;
    std::vector<int> flg;
    for (const auto& m : mag) {
      lmin.push_back(m.lambda_min);
      lam.push_back(m.lambda);
      lmax.push_back(m.lambda_max);
      fmin.push_back(m.flux_min);
      flx.push_back(m.flux);
      fmax.push_back(m.flux_max);
      dif.push_back(m.diff);
      derr.push_back(m.diff_err);
      pb.push_back(m.passband);
      sys.push_back(m.system);
      flg.push_back(m.flag);
      cat.push_back(m.vizier_catalog);
    }
    w.binary_table(
        "filters",
        {dcol("lambda_min", lmin), dcol("lambda", lam),
         dcol("lambda_max", lmax), dcol("flux_min", fmin), dcol("flux", flx),
         dcol("flux_max", fmax), dcol("diff", dif), dcol("diff_err", derr),
         scol("passband", pb), scol("system", sys), icol("flag", flg),
         scol("VizieR_catalog", cat)});
    if (!in.mout->col.empty()) {
      dvec cdif, cderr;
      std::vector<std::string> cpb, csys, ccat;
      std::vector<int> cflg;
      for (const auto& c : in.mout->col) {
        cdif.push_back(c.diff);
        cderr.push_back(c.diff_err);
        cpb.push_back(c.passband);
        csys.push_back(c.system);
        cflg.push_back(c.flag);
        ccat.push_back(c.vizier_catalog);
      }
      w.binary_table("colours",
                     {dcol("diff", cdif), dcol("diff_err", cderr),
                      scol("passband", cpb), scol("system", csys),
                      icol("flag", cflg), scol("VizieR_catalog", ccat)});
    }
    std::vector<FitsCol> spec = {dcol("l", *in.spec_l), dcol("f", *in.spec_f)};
    for (size_t i = 0; i < in.spec_fcomp->size(); ++i)
      spec.push_back(dcol("f_c" + std::to_string(i + 1), (*in.spec_fcomp)[i]));
    w.binary_table("spectrum_fit", spec);
  }

  if (!in.iue_wavelength.empty()) {
    w.binary_table("IUE", {dcol("wavelength", in.iue_wavelength),
                           dcol("flux", in.iue_flux),
                           dcol("error", in.iue_error)});
    w.key_string("FILE", in.iue_file);
  }

  // MC_c* keep the struct-field order R_Rsun, M_Msun, L_Lsun (to_keep applied
  // via struct_drop_fields preserves the MC struct layout). The per-quantity
  // positivity filters make the arrays ragged; short columns are padded with
  // NaN (ISIS behaviour unverified -- save_MC=0 in every fixture).
  for (size_t i = 0; i < in.mc.size() && i < 2; ++i) {
    const StellarMCResult* mc = in.mc[i];
    if (!mc || mc->R_arr.empty()) continue;
    dvec R = mc->R_arr, M = mc->M_arr, L = mc->L_arr;
    size_t n = std::max({R.size(), M.size(), L.size()});
    R.resize(n, NaN);
    M.resize(n, NaN);
    L.resize(n, NaN);
    w.binary_table("MC_c" + std::to_string(i + 1),
                   {dcol("R_Rsun", std::move(R)), dcol("M_Msun", std::move(M)),
                    dcol("L_Lsun", std::move(L))});
  }
  w.close();
}

void select_iue_spectrum(const std::vector<PhotEntry>& entries,
                         const std::string& workdir, FitsOutInput& out) {
  // union() of box passbands with flag==0 (S-Lang union sorts unique values)
  std::set<std::string> boxes;
  for (const auto& e : entries)
    if (e.system == "box" && e.flag == 0) boxes.insert(e.passband);
  auto file_from = [](const std::string& cat) {
    const std::string tag = "_VI/110/inescat";
    auto pos = cat.find(tag);
    if (pos == std::string::npos) return std::string();
    return "IUE/" + cat.substr(0, pos) + ".FITS.gz";
  };
  for (const auto& box : boxes) {
    for (size_t i = 0; i < entries.size(); ++i) {
      const auto& e = entries[i];
      if (e.passband != box || e.flag != 0) continue;
      std::string iue_file;
      if (e.vizier_catalog.find("VI/110/inescat") != std::string::npos) {
        iue_file = file_from(e.vizier_catalog);
      } else if (e.vizier_catalog == "Average") {
        // the flag==2 member whose magnitude is closest to the average
        double best = std::numeric_limits<double>::infinity();
        for (const auto& e2 : entries)
          if (e2.passband == box && e2.flag == 2 &&
              std::fabs(e.magnitude - e2.magnitude) < best) {
            best = std::fabs(e.magnitude - e2.magnitude);
            iue_file = file_from(e2.vizier_catalog);
          }
      }
      if (iue_file.empty()) continue;
      std::string p = workdir + "/" + iue_file;
      if (!std::filesystem::exists(p)) continue;
      try {
        dvec wl = fits_read_double_column(p + "[1]", "wavelength");
        dvec fl = fits_read_double_column(p + "[1]", "flux");
        dvec qu = fits_read_double_column(p + "[1]", "quality");
        dvec si = fits_read_double_column(p + "[1]", "sigma");
        dvec owl, ofl, osi;
        // bad-pixel removal (photometry.sl:1514)
        for (size_t k = 0; k < wl.size(); ++k)
          if (fl[k] > 0 && qu[k] == 0) {
            owl.push_back(wl[k]);
            ofl.push_back(fl[k]);
            osi.push_back(si[k]);
          }
        if (!owl.empty()) {
          out.iue_file = iue_file;
          out.iue_wavelength = std::move(owl);
          out.iue_flux = std::move(ofl);
          out.iue_error = std::move(osi);
        }
      } catch (const std::exception& ex) {
        std::fprintf(stderr, "Warning: cannot read IUE spectrum %s: %s\n",
                     p.c_str(), ex.what());
      }
    }
  }
}

}  // namespace sed
