#include "results.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <ostream>
#include <stdexcept>

namespace sed {

std::string slang_repr(double v) {
  // Reproduce S-Lang's %S / string(Double_Type): shortest round-trip digits,
  // with %g-style format selection (fixed when the decimal exponent E is in
  // [-4, 5], scientific otherwise), a decimal point always present in fixed
  // notation, and lowercase 'e' with a signed >=2-digit exponent.
  if (std::isnan(v)) return "-nan";
  if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
  if (v == 0.0) return std::signbit(v) ? "-0.0" : "0.0";

  // shortest scientific form: "[-]d[.ddd]e[+-]dd"
  char buf[64];
  auto res = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::scientific);
  std::string sci(buf, res.ptr);

  std::string sign;
  size_t start = 0;
  if (sci[0] == '-') {
    sign = "-";
    start = 1;
  }
  size_t epos = sci.find('e', start);
  std::string mant = sci.substr(start, epos - start);  // "d" or "d.ddd"
  int E = std::atoi(sci.c_str() + epos + 1);

  if (E < -4 || E >= 6)
    return sign + mant + sci.substr(epos);  // scientific, exponent as-is

  // digits without the decimal point
  std::string digits;
  for (char c : mant)
    if (c != '.') digits += c;

  int point_pos = E + 1;  // decimal point sits after this many leading digits
  std::string out;
  if (point_pos <= 0) {
    out = "0." + std::string(-point_pos, '0') + digits;
  } else if (point_pos >= int(digits.size())) {
    out = digits + std::string(point_pos - int(digits.size()), '0') + ".0";
  } else {
    out = digits.substr(0, point_pos) + "." + digits.substr(point_pos);
  }
  return sign + out;
}

void write_table(std::ostream& os, const std::vector<std::string>& header,
                 const std::vector<std::vector<std::string>>& columns) {
  const size_t ncol = columns.size();
  const size_t nrow = columns.empty() ? 0 : columns[0].size();
  std::vector<size_t> width(ncol);
  for (size_t c = 0; c < ncol; ++c) {
    width[c] = header[c].size();
    for (const auto& cell : columns[c]) width[c] = std::max(width[c], cell.size());
  }
  auto emit_row = [&](const std::vector<std::string>& cells) {
    std::string line;
    for (size_t c = 0; c < ncol; ++c) {
      if (c) line += "   ";
      std::string cell = cells[c];
      line += std::string(width[c] - cell.size(), ' ') + cell;
    }
    os << line << "\n";
  };
  emit_row(header);
  for (size_t r = 0; r < nrow; ++r) {
    std::vector<std::string> cells(ncol);
    for (size_t c = 0; c < ncol; ++c) cells[c] = columns[c][r];
    emit_row(cells);
  }
}

void write_results_txt(const std::string& path, const std::string& header,
                       const FitResults& r) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("cannot write " + path);
  os << header;
  const size_t n = r.name.size();
  std::vector<std::vector<std::string>> cols(11);
  for (size_t i = 0; i < n; ++i) {
    cols[0].push_back(std::to_string(r.index[i]));
    cols[1].push_back(r.name[i]);
    cols[2].push_back(slang_repr(r.value[i]));
    cols[3].push_back(std::to_string(r.freeze[i]));
    cols[4].push_back(slang_repr(r.min[i]));
    cols[5].push_back(slang_repr(r.max[i]));
    cols[6].push_back(slang_repr(r.conf_min[i]));
    cols[7].push_back(slang_repr(r.conf_max[i]));
    cols[8].push_back(slang_repr(r.buf_below[i]));
    cols[9].push_back(slang_repr(r.buf_above[i]));
    cols[10].push_back(r.tex[i]);
  }
  write_table(os,
              {"index", "name", "value", "freeze", "min", "max", "conf_min",
               "conf_max", "buf_below", "buf_above", "tex"},
              cols);
}

void write_stellar_txt(const std::string& path,
                       const std::vector<StellarRow>& rows) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("cannot write " + path);
  std::vector<std::vector<std::string>> cols(7);
  for (const auto& row : rows) {
    cols[0].push_back(row.name);
    cols[1].push_back(slang_repr(row.value));
    cols[2].push_back(slang_repr(row.conf_min));
    cols[3].push_back(slang_repr(row.conf_max));
    cols[4].push_back(slang_repr(row.median_value));
    cols[5].push_back(slang_repr(row.median_conf_min));
    cols[6].push_back(slang_repr(row.median_conf_max));
  }
  write_table(os,
              {"name", "value", "conf_min", "conf_max", "median_value",
               "median_conf_min", "median_conf_max"},
              cols);
}

void write_mag_txt(const std::string& path,
                   const std::vector<MagFluxRow>& rows) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("cannot write " + path);
  std::vector<std::vector<std::string>> cols(12);
  for (const auto& r : rows) {
    cols[0].push_back(slang_repr(r.lambda_min));
    cols[1].push_back(slang_repr(r.lambda));
    cols[2].push_back(slang_repr(r.lambda_max));
    cols[3].push_back(slang_repr(r.flux_min));
    cols[4].push_back(slang_repr(r.flux));
    cols[5].push_back(slang_repr(r.flux_max));
    cols[6].push_back(slang_repr(r.diff));
    cols[7].push_back(slang_repr(r.diff_err));
    cols[8].push_back(r.passband);
    cols[9].push_back(r.system);
    cols[10].push_back(std::to_string(r.flag));
    cols[11].push_back(r.vizier_catalog);
  }
  write_table(os,
              {"lambda_min", "lambda", "lambda_max", "flux_min", "flux",
               "flux_max", "diff", "diff_err", "passband", "system", "flag",
               "VizieR_catalog"},
              cols);
}

void write_col_txt(const std::string& path,
                   const std::vector<ColFluxRow>& rows) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("cannot write " + path);
  std::vector<std::vector<std::string>> cols(6);
  for (const auto& r : rows) {
    cols[0].push_back(slang_repr(r.diff));
    cols[1].push_back(slang_repr(r.diff_err));
    cols[2].push_back(r.passband);
    cols[3].push_back(r.system);
    cols[4].push_back(std::to_string(r.flag));
    cols[5].push_back(r.vizier_catalog);
  }
  write_table(
      os, {"diff", "diff_err", "passband", "system", "flag", "VizieR_catalog"},
      cols);
}

void write_spectrum_txt(const std::string& path, const dvec& l, const dvec& f,
                        const std::vector<dvec>& f_comp) {
  std::ofstream os(path);
  if (!os) throw std::runtime_error("cannot write " + path);
  std::vector<std::string> header = {"l", "f"};
  std::vector<std::vector<std::string>> cols;
  cols.emplace_back();
  cols.emplace_back();
  for (size_t i = 0; i < l.size(); ++i) {
    cols[0].push_back(slang_repr(l[i]));
    cols[1].push_back(slang_repr(f[i]));
  }
  for (size_t c = 0; c < f_comp.size(); ++c) {
    header.push_back("f_c" + std::to_string(c + 1));
    cols.emplace_back();
    for (double v : f_comp[c]) cols.back().push_back(slang_repr(v));
  }
  write_table(os, header, cols);
}

}  // namespace sed
