#include "photometry_table.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "results.hpp"
#include "util.hpp"

namespace sed {

namespace {
double parse_double(const std::string& s) {
  // atof semantics: "-nan", "nan", "-999.0" etc.
  return std::atof(s.c_str());
}
}  // namespace

PhotometryTable PhotometryTable::read(const std::string& filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("cannot open " + filename);
  PhotometryTable t;
  std::string line;
  bool header_done = false;
  bool column_line_done = false;
  while (std::getline(in, line)) {
    if (!header_done && !line.empty() && line[0] == '#') {
      double a, b;
      if (2 == std::sscanf(line.c_str(), "# RA = %lf DEC = %lf", &a, &b)) {
        t.ra = a;
        t.dec = b;
      } else if (2 == std::sscanf(line.c_str(), "# meanSFD = %lf stdSFD = %lf",
                                  &a, &b)) {
        t.reddening.meanSFD = a;
        t.reddening.stdSFD = b;
      } else if (2 == std::sscanf(line.c_str(),
                                  "# meanSandF = %lf stdSandF = %lf", &a,
                                  &b)) {
        t.reddening.meanSandF = a;
        t.reddening.stdSandF = b;
      } else if (2 == std::sscanf(line.c_str(),
                                  "# meanStilism = %lf stdStilism = %lf", &a,
                                  &b)) {
        t.reddening.meanStilism = a;
        t.reddening.stdStilism = b;
      }
      continue;
    }
    header_done = true;
    std::string s = trim(line);
    if (s.empty()) continue;
    if (!column_line_done) {  // column-name header line
      column_line_done = true;
      continue;
    }
    auto tok = split_ws(s);
    if (tok.size() != 8)
      throw std::runtime_error("bad photometry.dat row: " + line);
    PhotEntry e;
    e.flag = std::atoi(tok[0].c_str());
    e.system = tok[1];
    e.passband = tok[2];
    e.magnitude = parse_double(tok[3]);
    e.uncertainty = parse_double(tok[4]);
    e.type = tok[5];
    e.angu_dist_arcsec = parse_double(tok[6]);
    e.vizier_catalog = tok[7];
    t.entries.push_back(std::move(e));
  }
  return t;
}

void PhotometryTable::write(const std::string& filename) const {
  std::ofstream os(filename);
  if (!os) throw std::runtime_error("cannot write " + filename);
  char buf[128];
  std::snprintf(buf, sizeof(buf), "# RA = %.10f DEC = %.10f\n",
                ra.value_or(0.0), dec.value_or(0.0));
  os << buf;
  if (reddening.meanSFD) {
    std::snprintf(buf, sizeof(buf), "# meanSFD = %.4f stdSFD = %.4f\n",
                  *reddening.meanSFD, reddening.stdSFD.value_or(0.0));
    os << buf;
  }
  if (reddening.meanSandF) {
    std::snprintf(buf, sizeof(buf), "# meanSandF = %.4f stdSandF = %.4f\n",
                  *reddening.meanSandF, reddening.stdSandF.value_or(0.0));
    os << buf;
  }
  if (reddening.meanStilism) {
    std::snprintf(buf, sizeof(buf), "# meanStilism = %.4f stdStilism = %.4f\n",
                  *reddening.meanStilism, reddening.stdStilism.value_or(0.0));
    os << buf;
  }
  std::vector<std::vector<std::string>> cols(8);
  for (const auto& e : entries) {
    cols[0].push_back(std::to_string(e.flag));
    cols[1].push_back(e.system);
    cols[2].push_back(e.passband);
    cols[3].push_back(slang_repr(e.magnitude));
    cols[4].push_back(slang_repr(e.uncertainty));
    cols[5].push_back(e.type);
    cols[6].push_back(slang_repr(e.angu_dist_arcsec));
    cols[7].push_back(e.vizier_catalog);
  }
  write_table(os, {"flag", "system", "passband", "magnitude", "uncertainty",
                   "type", "angu_dist_arcsec", "VizieR_catalog"},
              cols);
}

double ResultsHeader::num(const std::string& k) const {
  auto it = raw.find(k);
  if (it == raw.end()) throw std::runtime_error("missing header key " + k);
  return parse_double(it->second);
}

std::optional<ResultsHeader> ResultsHeader::read(const std::string& filename) {
  std::ifstream in(filename);
  if (!in) return std::nullopt;
  ResultsHeader h;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] != '#') break;
    // "# key = value; key = value;" -- mirror read_photometry_header: pad '='
    // with spaces, tokenize, take (key '=' value) triples, strip ';'.
    std::string body = line.substr(1);
    std::string padded;
    for (char c : body) {
      if (c == '=')
        padded += " = ";
      else
        padded += c;
    }
    auto tok = split_ws(padded);
    for (size_t i = 0; i + 2 < tok.size(); i += 3) {
      if (tok[i + 1] == "=") {
        std::string val = tok[i + 2];
        if (!val.empty() && val.back() == ';') val.pop_back();
        h.raw[tok[i]] = val;
      }
    }
  }
  if (h.raw.empty()) return std::nullopt;
  return h;
}

}  // namespace sed
