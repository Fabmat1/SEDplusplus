#include "config.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace sed {

using nlohmann::json;

Config Config::load(const std::string& json_path) {
  std::ifstream in(json_path);
  if (!in) throw std::runtime_error("cannot open config " + json_path);
  json j;
  in >> j;

  Config c;
  c.star = j.at("star").get<std::string>();
  c.basename = j.value("basename", std::string());
  if (j.contains("coordinates") && !j["coordinates"].is_null()) {
    const auto& co = j["coordinates"];
    if (co.contains("ra") && !co["ra"].is_null())
      c.ra = co["ra"].get<double>();
    if (co.contains("dec") && !co["dec"].is_null())
      c.dec = co["dec"].get<double>();
  }
  if (j.contains("fix_distance") && !j["fix_distance"].is_null())
    c.fix_distance = j["fix_distance"].get<double>();
  if (j.contains("fix_distance_err") && !j["fix_distance_err"].is_null())
    c.fix_distance_err = j["fix_distance_err"].get<double>();

  auto load_parset = [&](const char* key, bool full) {
    ParSet ps;
    if (!j.contains(key)) return ps;
    const auto& p = j[key];
    ps.name = p.at("name").get<std::vector<std::string>>();
    ps.value = p.at("value").get<std::vector<double>>();
    ps.freeze = p.at("freeze").get<std::vector<int>>();
    if (full) {
      ps.min = p.at("min").get<std::vector<double>>();
      ps.max = p.at("max").get<std::vector<double>>();
    }
    size_t n = ps.name.size();
    if (ps.value.size() != n || ps.freeze.size() != n ||
        (full && (ps.min.size() != n || ps.max.size() != n)))
      throw std::runtime_error(std::string(key) + ": inconsistent lengths");
    return ps;
  };
  c.par = load_parset("par", false);
  c.par_full = load_parset("par_full", true);

  c.griddirectories = j.at("griddirectories").get<std::vector<std::string>>();
  c.bpaths = j.value("bpaths", std::vector<std::string>{"./"});

  c.conf_level = j.value("conf_level", 0);
  c.max_conf_restarts = j.value("max_conf_restarts", 1000);
  // Phase-3 toggles: default OFF (bulk fitting stays lean). Only write_model
  // (and plot, which implies it) has behaviour in Stage 1; the rest are
  // parse-only until later stages.
  c.write_model = j.value("write_model", 0) != 0;
  c.write_fits = j.value("write_fits", 0) != 0;
  c.write_tex = j.value("write_tex", 0) != 0;
  c.plot = j.value("plot", 0) != 0;
  c.save_MC = j.value("save_MC", 0) != 0;
  if (c.plot) {  // plot implies write_fits + write_model (photometry.sl)
    c.write_fits = true;
    c.write_model = true;
  }
  c.plot_script = j.value("plot_script", std::string());
  c.apply_ZPO_corr = j.value("apply_ZPO_corr", 1) != 0;
  c.remove_outliers = j.value("remove_outliers", 5.0);
  c.nMC = j.value("nMC", 2000000L);
  c.stilism_distance_simple = j.value("stilism_distance_simple", 1);
  c.stilism_ebmv_simple = j.value("stilism_ebmv_simple", 1);
  c.stilism_ebmv_rerun = j.value("stilism_ebmv_rerun", 1);
  c.mass_can = j.value("mass_can", 0.0);
  c.delta_mass_can = j.value("delta_mass_can", 0.05);
  c.derive_logg = j.value("derive_logg", 0);
  c.hb_distance = j.value("hb_distance", 0);
  if (c.hb_distance) c.derive_logg = 1;
  c.derive_logg_c2 = j.value("derive_logg_c2", 0);
  c.z_c2 = j.value("z_c2", -0.9);
  c.derive_sr = j.value("derive_sr", 0);
  c.sdOB_radius = j.value("sdOB_radius", 0.2);
  c.R1 = j.value("R1", 0.0);
  c.R1_err = j.value("R1_err", 0.01);

  c.refdata = j.at("refdata").get<std::string>();
  c.workdir = j.value("workdir", std::string("."));
  c.outdir = j.value("outdir", c.workdir);
  c.mc_seed = j.value("mc_seed", 42UL);

  // niche options not yet supported in the C++ port
  if (c.mass_can > 0 || c.derive_logg || c.hb_distance || c.derive_logg_c2 ||
      c.derive_sr || c.R1 > 0)
    throw std::runtime_error(
        "mass_can/derive_logg/hb_distance/derive_logg_c2/derive_sr/R1 are not "
        "implemented in this version");
  return c;
}

}  // namespace sed
