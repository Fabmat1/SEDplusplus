#include "texout.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#include "texval.hpp"

namespace sed {

namespace {

int idx_of(const FitResults& r, const std::string& nm) {
  for (size_t i = 0; i < r.name.size(); ++i)
    if (r.name[i] == nm) return int(i);
  return -1;
}

// One parameter row (photometry.sl:903-922 top block, 986-1014 component
// fields): the tex column when set; otherwise "(prescribed)" when par_full
// holds the *literal* name frozen (wildcards like "c*_teff" do not match,
// same as where(par_full.name==sprintf(...))); otherwise free ($%g$, red
// when at a bound) or "(fixed) & $%g$".
void param_row(std::FILE* fp, const FitResults& r, const ParSet* par_full,
               bool check_prescribed, const std::string& pname,
               const char* label, const char* unit) {
  int i = idx_of(r, pname);
  if (i < 0) return;
  if (r.tex[i] != "\\ldots") {
    std::fprintf(fp, "%s & %s%s \\\\\n", label, r.tex[i].c_str(), unit);
    return;
  }
  if (check_prescribed && par_full) {
    int t = -1, nmatch = 0;
    for (size_t k = 0; k < par_full->name.size(); ++k)
      if (par_full->name[k] == pname) {
        if (t < 0) t = int(k);
        ++nmatch;
      }
    if (nmatch == 1 && par_full->freeze[t] == 1) {
      std::fprintf(fp, "%s (prescribed) & %s%s \\\\\n", label,
                   tex_value_pm_error(par_full->value[t], par_full->min[t],
                                      par_full->max[t], 6)
                       .c_str(),
                   unit);
      return;
    }
  }
  if (r.freeze[i] == 0) {
    if (r.value[i] == r.min[i] || r.value[i] == r.max[i])
      std::fprintf(fp, "%s & \\color{red}{$%g$}%s \\\\\n", label, r.value[i],
                   unit);
    else
      std::fprintf(fp, "%s & $%g$%s \\\\\n", label, r.value[i], unit);
  } else {
    std::fprintf(fp, "%s (fixed) & $%g$%s \\\\\n", label, r.value[i], unit);
  }
}

std::string tex_mode(const ModeHDI& m, int sci) {
  return tex_value_pm_error(m.mode, m.HDI_lo, m.HDI_hi, sci);
}
std::string tex_median(const ModeHDI& m, int sci) {
  return tex_value_pm_error(m.median, m.quantile_lo, m.quantile_hi, sci);
}

}  // namespace

void write_results_tex(const std::string& path, const TexInput& in) {
  const FitResults& r = *in.res;
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot open " + path);

  std::fputs(
      "\\documentclass{standalone}\n"
      "\\usepackage{amsmath,txfonts,color}\n"
      "\\usepackage[colorlinks=true,urlcolor=black]{hyperref}\n"
      "\\begin{document}\n"
      "\\renewcommand{\\arraystretch}{1.2}\n",
      fp);
  std::string headline = "Object: " + in.star + " & ";
  if (in.conf_level == -1 || in.conf_level == 0)
    headline += "68\\% confidence interval";
  else if (in.conf_level == 1)
    headline += "90\\% confidence interval";
  else if (in.conf_level == 2)
    headline += "99\\% confidence interval";
  else
    headline += " Value";
  std::fprintf(fp, "\\begin{tabular}{lr}\n\\hline\\hline\n%s\\\\\n\\hline\n",
               headline.c_str());

  // reddening rows (photometry.sl:873-890); mean +- std with default sci=2
  auto redd_row = [&](const char* fmt, const std::optional<double>& mean,
                      const std::optional<double>& std_) {
    if (!mean) return;
    double s = std_.value_or(0.0);
    std::fprintf(fp, fmt,
                 tex_value_pm_error(*mean, *mean - s, *mean + s, 2).c_str());
  };
  redd_row(
      "Color excess $E(B-V)$ from "
      "\\href{https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract}"
      "{SFD (1998)} & %s\\,mag \\\\\n",
      in.reddening.meanSFD, in.reddening.stdSFD);
  redd_row(
      "Color excess $E(B-V)$ from "
      "\\href{https://ui.adsabs.harvard.edu/abs/2011ApJ...737..103S/abstract}"
      "{S\\&F (2011)} & %s\\,mag \\\\\n",
      in.reddening.meanSandF, in.reddening.stdSandF);
  redd_row(
      "Color excess $E(B-V)$ from \\href{https://stilism.obspm.fr/}{Stilism} "
      "\\href{https://ui.adsabs.harvard.edu/abs/2017A%%26A...606A..65C/"
      "abstract}{(Capitanio+ 2017)} & %s\\,mag \\\\\n",
      in.reddening.meanStilism, in.reddening.stdStilism);
  std::fputs("\\hline\n", fp);

  param_row(fp, r, nullptr, false, "E_44m55", "Color excess $E(44-55)$",
            "\\,mag");
  param_row(fp, r, nullptr, false, "R_55", "Extinction parameter $R(55)$", "");
  param_row(fp, r, nullptr, false, "logtheta",
            "Angular diameter $\\log(\\Theta\\,\\mathrm{(rad)})$", "");

  if (in.good_astrometry) {
    // no parallax_offset field in our Astrometry (photometry.sl:931-934)
    std::fprintf(
        fp,
        "Parallax $\\varpi$ ({\\it Gaia}, $\\text{RUWE}=%.2f$) & %s\\,mas "
        "\\\\\n",
        in.ruwe,
        tex_value_pm_error(in.parallax, in.parallax - in.parallax_error,
                           in.parallax + in.parallax_error, 2)
            .c_str());
    if (in.distance) {
      std::fprintf(fp, "Distance $d$ ({\\it Gaia}, mode) & %s\\,pc \\\\\n",
                   tex_mode(*in.distance, 3).c_str());
      std::fprintf(fp, "Distance $d$ ({\\it Gaia}, median) & %s\\,pc \\\\\n",
                   tex_median(*in.distance, 3).c_str());
    }
  }

  // per-component blocks (photometry.sl:973-1238); our fit function is
  // always the fit_theta_interext_atmparams flavour
  struct Field {
    const char* name;
    const char* label;
    const char* unit;
  };
  static const Field fields[] = {
      {"teff", "Effective temperature $T_{\\mathrm{eff}}$", "\\,K"},
      {"logg", "Surface gravity $\\log (g\\,\\mathrm{(cm\\,s^{-2})})$", ""},
      {"xi", "Microturbulence $\\xi$", "\\,km\\,s$^{-1}$"},
      {"z", "Metallicity $z$", "\\,dex"},
      {"HE", "Helium abundance $\\log(n(\\textnormal{He}))$", ""},
      {"sur_ratio", "Surface ratio $A_\\mathrm{eff}/A_{\\mathrm{eff,}1}$", ""},
  };
  for (int comp = 1; comp <= int(in.ncomp); ++comp) {
    if (in.ncomp > 1)
      std::fprintf(fp, "\\hline\n\\multicolumn{2}{l}{Component %d:} \\\\\n",
                   comp);
    std::string cstr = "c" + std::to_string(comp);
    for (const Field& f : fields) {
      // surface ratio of component 1 is fixed to 1, never printed
      if (comp == 1 && std::string(f.name) == "sur_ratio") continue;
      param_row(fp, r, in.par_full, true, cstr + "_" + f.name, f.label,
                f.unit);
    }
    const StellarMCResult* mc =
        (in.mc && comp <= int(in.mc->size())) ? &(*in.mc)[comp - 1] : nullptr;
    if (in.good_astrometry && mc && mc->valid) {
      const char* rlab =
          comp == 1 ? "" : "(A_\\mathrm{eff}/A_{\\mathrm{eff,}1})^{1/2}";
      std::fprintf(
          fp, "Radius $R = %s\\Theta/(2\\varpi)$ (mode) & %s\\,$R_\\odot$ \\\\\n",
          rlab, tex_mode(mc->R, 3).c_str());
      std::fprintf(
          fp,
          "\\phantom{Radius $R = %s\\Theta/(2\\varpi)$ }(median) & "
          "%s\\,$R_\\odot$ \\\\\n",
          rlab, tex_median(mc->R, 3).c_str());
      std::fprintf(fp, "Mass $M = g R^2/G$ (mode) & %s\\,$M_\\odot$ \\\\\n",
                   tex_mode(mc->M, 3).c_str());
      std::fprintf(
          fp, "\\phantom{Mass $M = g R^2/G$ }(median) & %s\\,$M_\\odot$ \\\\\n",
          tex_median(mc->M, 3).c_str());
      std::fprintf(
          fp,
          "Luminosity $L = "
          "(R/R_\\odot)^2(T_\\mathrm{eff}/T_{\\mathrm{eff},\\odot})^4$ (mode) "
          "& %s\\,$L_\\odot$ \\\\\n",
          tex_mode(mc->L, 3).c_str());
      std::fprintf(
          fp,
          "\\phantom{Luminosity $L/L_\\odot = "
          "(R/R_\\odot)^2(T_\\mathrm{eff}/T_{\\mathrm{eff},\\odot})^4$ "
          "}(median) & %s \\\\\n",
          tex_median(mc->L, 3).c_str());
      std::fprintf(
          fp,
          "Gravitational redshift $\\varv_\\mathrm{grav} = GM/(Rc)$ & "
          "%s\\,km\\,s${}^{-1}$ \\\\\n",
          tex_mode(mc->vgrav, 3).c_str());
      std::fprintf(
          fp,
          "Escape velocity $\\varv_\\mathrm{esc} = \\sqrt{2gR}$ & "
          "%s\\,km\\,s${}^{-1}$ \\\\\n",
          tex_mode(mc->vesc, 3).c_str());
    }
  }
  // Blackbody blocks (photometry.sl:1240-1295): bb components are not
  // supported by this port, so the gates bb_teff!=0 && bb_sur_ratio!=0 are
  // always false and nothing is written.

  std::fputs("\\hline\n", fp);
  std::fprintf(
      fp,
      "Generic excess noise $\\delta_\\textnormal{excess}$ & $%.3f$\\,mag "
      "\\\\\n",
      in.norm_chi_red);
  std::fprintf(fp, "Reduced $\\chi^2$ at the best fit & $%.2f$ \\\\\n",
               in.chisqr_red);
  std::fputs("\\hline\n\\end{tabular}\n\\end{document}\n", fp);
  std::fclose(fp);
}

bool run_pdflatex(const std::string& outdir, const std::string& texname) {
  if (std::system("command -v pdflatex >/dev/null 2>&1") != 0) {
    std::fprintf(stderr,
                 "Warning: pdflatex not found in PATH; skipping "
                 "photometry_results.pdf\n");
    return false;
  }
  std::string stem = texname.substr(0, texname.rfind('.'));
  std::string dir = outdir.empty() ? "." : outdir;
  // photometry.sl:1303, run inside outdir so the .pdf lands there;
  // </dev/null guards against pdflatex prompting despite -halt-on-error
  std::string cmd = "cd '" + dir +
                    "' && pdflatex -halt-on-error -file-line-error '" +
                    texname + "' </dev/null | grep 'error' --color=always; " +
                    "rm -f '" + stem + ".aux' '" + stem + ".log' '" + stem +
                    ".out'";
  return std::system(cmd.c_str()) == 0;
}

}  // namespace sed
