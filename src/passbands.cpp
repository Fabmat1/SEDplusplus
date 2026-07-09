#include "passbands.hpp"

#include <fitsio.h>

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "zp_table.inc"

namespace sed {

namespace {

// filter label aliases from filter_passbands.sl
std::string replace_all(std::string s, const std::string& a,
                        const std::string& b) {
  size_t pos = 0;
  while ((pos = s.find(a, pos)) != std::string::npos) {
    s.replace(pos, a.size(), b);
    pos += b.size();
  }
  return s;
}

bool file_exists(const std::string& p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0;
}

}  // namespace

PassbandDB::PassbandDB(const std::string& refdata_dir, bool apply_ZPO_corr,
                       const std::vector<std::string>& boxes)
    : refdata_(refdata_dir) {
  // zero-point table: refdata/zeropoint_offsets.txt when present (same row
  // format as the builtin table; '#' comments and the "system passband ..."
  // header line skipped), otherwise the builtin zp_table.inc
  std::vector<std::string> rows;
  const std::string zpo_file = refdata_ + "/zeropoint_offsets.txt";
  if (file_exists(zpo_file)) {
    std::ifstream in(zpo_file);
    std::string line;
    while (std::getline(in, line)) {
      size_t hash = line.find('#');
      if (hash != std::string::npos) line.erase(hash);
      size_t first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos) continue;
      if (line.compare(first, 6, "system") == 0) continue;  // header line
      rows.push_back(line);
    }
    if (rows.empty())
      throw std::runtime_error("no zero-point rows in " + zpo_file);
  } else {
    for (const char* row : ZP_TABLE) rows.push_back(row);
  }
  for (const std::string& row : rows) {
    BandEntry e;
    char sys[64], pb[64], type[32];
    double zp, zperr, corr, correrr;
    if (7 != std::sscanf(row.c_str(), "%63s %63s %lf %lf %lf %lf %31s", sys,
                         pb, &zp, &zperr, &corr, &correrr, type))
      throw std::runtime_error(std::string("bad ZP table row: ") + row);
    e.system = sys;
    e.passband = pb;
    e.ZP = zp;
    e.ZP_err = zperr;
    e.type = type;
    if (apply_ZPO_corr) {
      e.ZP = zp - corr;
      e.ZP_err = std::sqrt(zperr * zperr + correrr * correrr);
    }
    entries_.push_back(std::move(e));
  }
  for (const auto& b : boxes) {
    BandEntry e;
    e.system = "box";
    e.passband = b;
    e.ZP = 0.0;
    e.ZP_err = 0.015;  // boxes_ZP_err (same value below/above 3660 A)
    e.type = "magnitude";
    entries_.push_back(std::move(e));
  }
  n_mag_ = entries_.size();
  append_colors();
}

int PassbandDB::find(const std::string& system,
                     const std::string& passband) const {
  for (size_t i = 0; i < entries_.size(); ++i)
    if (entries_[i].system == system && entries_[i].passband == passband)
      return int(i);
  return -1;
}

const FilterCurve& PassbandDB::filter(size_t idx) const {
  const BandEntry& e = entries_.at(idx);
  std::string label = e.system + "_" + e.passband;
  auto it = curves_.find(label);
  if (it != curves_.end()) return it->second;

  FilterCurve c;
  bool is_td1_box = (e.system == "TD1" && (e.passband == "1565" ||
                                           e.passband == "1965" ||
                                           e.passband == "2365"));
  if (e.system == "box" || is_td1_box) {
    double l0, l1;
    if (e.system == "box") {
      if (2 != std::sscanf(e.passband.c_str(), "%lf_%lf", &l0, &l1))
        throw std::runtime_error("invalid box passband: " + e.passband);
    } else {
      double m = std::atof(e.passband.c_str());
      l0 = m - 165.0;
      l1 = m + 165.0;
    }
    double norm = 1.0 / (1.0 / l0 - 1.0 / l1);
    c.l = {l0, l1};
    c.f = {norm, norm};
  } else {
    std::string flabel = label;
    flabel = replace_all(flabel, "DECaPS", "DES");
    flabel = replace_all(flabel, "DES_DR2", "DES");
    flabel = replace_all(flabel, "DELVE_DR1", "DES");
    flabel = replace_all(flabel, "DELVE_DR2", "DES");
    flabel = replace_all(flabel, "SMASH", "SDSS");
    flabel = replace_all(flabel, "HST05", "HST");
    flabel = replace_all(flabel, "HST18", "HST");
    flabel = replace_all(flabel, "VegaHST", "HST");
    flabel = replace_all(flabel, "VegaINT", "INT");
    load_archive();
    auto ait = archive_.find(flabel);
    // DECam vs DES: newer archives name the filter-system extensions DECam_*,
    // older ones DES_*; accept either spelling for either system name
    if (ait == archive_.end() && flabel.compare(0, 4, "DES_") == 0)
      ait = archive_.find("DECam" + flabel.substr(3));
    if (ait == archive_.end() && flabel.compare(0, 6, "DECam_") == 0)
      ait = archive_.find("DES" + flabel.substr(5));
    if (ait == archive_.end())
      throw std::runtime_error("cannot open filter " + refdata_ +
                               "/filter_passbands.fits[.gz][" + flabel + "]");
    c = ait->second;
  }
  return curves_.emplace(label, std::move(c)).first->second;
}

void PassbandDB::load_archive() const {
  if (archive_loaded_) return;
  archive_loaded_ = true;
  // process-wide cache keyed by refdata path, so --multi decompresses the
  // archive once per process instead of once per star (single-threaded)
  static std::map<std::string, std::map<std::string, FilterCurve>> g_archives;
  auto git = g_archives.find(refdata_);
  if (git != g_archives.end()) {
    archive_ = git->second;
    return;
  }
  std::string path = refdata_ + "/filter_passbands.fits.gz";
  if (!file_exists(path)) {
    const std::string plain = refdata_ + "/filter_passbands.fits";
    if (file_exists(plain)) path = plain;
  }
  fitsfile* fp = nullptr;
  int status = 0;
  if (fits_open_file(&fp, path.c_str(), READONLY, &status))
    throw std::runtime_error("cannot open filter archive " + path);
  int nhdus = 0;
  fits_get_num_hdus(fp, &nhdus, &status);
  for (int hdu = 2; hdu <= nhdus && !status; ++hdu) {
    int hdutype = 0;
    fits_movabs_hdu(fp, hdu, &hdutype, &status);
    char extname[FLEN_VALUE] = "";
    int kstat = 0;
    fits_read_key(fp, TSTRING, "EXTNAME", extname, nullptr, &kstat);
    if (kstat || !extname[0]) continue;
    int lcol = 0, fcol = 0, cstat = 0;
    fits_get_colnum(fp, CASEINSEN, const_cast<char*>("l"), &lcol, &cstat);
    fits_get_colnum(fp, CASEINSEN, const_cast<char*>("f"), &fcol, &cstat);
    if (cstat) continue;  // not a filter-curve table
    long nrows = 0;
    fits_get_num_rows(fp, &nrows, &status);
    FilterCurve c;
    c.l.resize(nrows);
    c.f.resize(nrows);
    double nulval = 0.0;
    int anynul = 0;
    fits_read_col(fp, TDOUBLE, lcol, 1, 1, nrows, &nulval, c.l.data(), &anynul,
                  &status);
    fits_read_col(fp, TDOUBLE, fcol, 1, 1, nrows, &nulval, c.f.data(), &anynul,
                  &status);
    archive_.emplace(extname, std::move(c));
  }
  fits_close_file(fp, &status);
  if (status) throw std::runtime_error("error reading " + path);
  g_archives.emplace(refdata_, archive_);
}

void PassbandDB::append_colors() {
  struct ColorDef {
    const char* system;
    const char* passband;
  };
  static const ColorDef defs[21] = {
      {"Johnson", "UmB"},    {"Johnson", "BmV"},
      {"Johnson", "VmR"},    {"Johnson", "VmI"},
      {"Johnson", "RmI"},    {"Stroemgren", "bmy"},
      {"Stroemgren", "umb"}, {"Stroemgren", "c1"},
      {"Stroemgren", "m1"},  {"Stroemgren", "Hbeta_B"},
      {"Stroemgren", "Hbeta_AF"}, {"Geneva", "UmB"},
      {"Geneva", "B1mB"},    {"Geneva", "B2mB"},
      {"Geneva", "V1mB"},    {"Geneva", "VmB"},
      {"Geneva", "GmB"},     {"Geneva", "X"},
      {"Geneva", "Y"},       {"Geneva", "Z"},
      {"AstroSat", "F148WmF169M"},
  };
  auto midx = [&](const char* sys, const char* pb) {
    int i = find(sys, pb);
    if (i < 0) throw std::runtime_error("missing magnitude row");
    return i;
  };
  const int U = midx("Johnson", "U"), B = midx("Johnson", "B"),
            V = midx("Johnson", "V"), R = midx("Johnson", "R"),
            I = midx("Johnson", "I");
  const int sb = midx("Stroemgren", "b"), sy = midx("Stroemgren", "y"),
            sv = midx("Stroemgren", "v"), su = midx("Stroemgren", "u"),
            hn = midx("Stroemgren", "Hbeta_narrow"),
            hw = midx("Stroemgren", "Hbeta_wide");
  const int GU = midx("Geneva", "U"), GB = midx("Geneva", "B"),
            GV = midx("Geneva", "V"), GG = midx("Geneva", "G"),
            GB1 = midx("Geneva", "B1"), GB2 = midx("Geneva", "B2"),
            GV1 = midx("Geneva", "V1");
  const int AF148 = midx("AstroSat", "F148W"), AF169 = midx("AstroSat", "F169M");
  cidx_ = ColorIdx{U,  B,  V,  R,  I,   sb,  sy,  sv,    su,   hn,
                   hw, GU, GB, GV, GG,  GB1, GB2, GV1,   AF148, AF169};

  // dependency lists (mag_depend), in the same row order as defs[]
  deps_ = {
      {U, B},        {B, V},          {V, R},
      {V, I},        {R, I},          {sb, sy},
      {su, sb},      {su, sv, sb},    {sv, sb, sy},
      {hn, hw},      {hn, hw},        {GU, GB},
      {GB1, GB},     {GB2, GB},       {GV1, GB},
      {GV, GB},      {GG, GB},        {GU, GB, GB1, GB2, GV1, GG},
      {GU, GB, GB1, GB2, GV1, GG},    {GU, GB, GB1, GB2, GV1, GG},
      {AF148, AF169},
  };

  auto zp = [&](int i) { return entries_[i].ZP; };
  auto zperr = [&](int i) { return entries_[i].ZP_err; };
  auto q = [](std::initializer_list<double> v) {
    double s = 0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
  };

  for (const auto& d : defs) {
    BandEntry e;
    e.system = d.system;
    e.passband = d.passband;
    e.type = "color";
    entries_.push_back(e);
  }
  // color zero points (linear combinations of magnitude ZPs)
  auto& E = entries_;
  const size_t c0 = n_mag_;
  auto setzp = [&](size_t k, double z, double ze) {
    E[c0 + k].ZP = z;
    E[c0 + k].ZP_err = ze;
  };
  setzp(0, zp(U) - zp(B), q({zperr(U), zperr(B)}));
  setzp(1, zp(B) - zp(V), q({zperr(B), zperr(V)}));
  setzp(2, zp(V) - zp(R), q({zperr(V), zperr(R)}));
  setzp(3, zp(V) - zp(I), q({zperr(V), zperr(I)}));
  setzp(4, zp(R) - zp(I), q({zperr(R), zperr(I)}));
  setzp(5, zp(sb) - zp(sy), q({zperr(sb), zperr(sy)}));
  setzp(6, zp(su) - zp(sb), q({zperr(su), zperr(sb)}));
  setzp(7, zp(su) - 2 * zp(sv) + zp(sb), q({zperr(su), 2 * zperr(sv), zperr(sb)}));
  setzp(8, zp(sv) - 2 * zp(sb) + zp(sy), q({zperr(sv), 2 * zperr(sb), zperr(sy)}));
  double hzp = zp(hn) - zp(hw);
  double hzperr = q({zperr(hn), zperr(hw)});
  setzp(9, 1.305 * hzp, 1.305 * hzperr);
  setzp(10, 1.368 * hzp, 1.368 * hzperr);
  double gUmB = zp(GU) - zp(GB), gB1mB = zp(GB1) - zp(GB),
         gB2mB = zp(GB2) - zp(GB), gV1mB = zp(GV1) - zp(GB),
         gVmB = zp(GV) - zp(GB), gGmB = zp(GG) - zp(GB);
  double eUmB = q({zperr(GU), zperr(GB)}), eB1mB = q({zperr(GB1), zperr(GB)}),
         eB2mB = q({zperr(GB2), zperr(GB)}), eV1mB = q({zperr(GV1), zperr(GB)}),
         eVmB = q({zperr(GV), zperr(GB)}), eGmB = q({zperr(GG), zperr(GB)});
  setzp(11, gUmB, eUmB);
  setzp(12, gB1mB, eB1mB);
  setzp(13, gB2mB, eB2mB);
  setzp(14, gV1mB, eV1mB);
  setzp(15, gVmB, eVmB);
  setzp(16, gGmB, eGmB);
  setzp(17,
        1.3764 * gUmB - 1.2162 * gB1mB - 0.8498 * gB2mB - 0.1554 * gV1mB +
            0.8450 * gGmB,
        q({1.3764 * eUmB, 1.2162 * eB1mB, 0.8498 * eB2mB, 0.1554 * eV1mB,
           0.8450 * eGmB}));
  setzp(18,
        0.3235 * gUmB - 2.3228 * gB1mB + 2.3363 * gB2mB + 0.7495 * gV1mB -
            1.0865 * gGmB,
        q({0.3235 * eUmB, 2.3228 * eB1mB, 2.3363 * eB2mB, 0.7495 * eV1mB,
           1.0865 * eGmB}));
  setzp(19,
        0.0255 * gUmB - 0.1740 * gB1mB + 0.4696 * gB2mB - 1.1205 * gV1mB +
            0.7994 * gGmB,
        q({0.0255 * eUmB, 0.1740 * eB1mB, 0.4696 * eB2mB, 1.1205 * eV1mB,
           0.7994 * eGmB}));
  setzp(20, zp(AF148) - zp(AF169), q({zperr(AF148), zperr(AF169)}));
}

void PassbandDB::compute_colors(dvec& mags) const {
  auto m = [&](int i) { return mags[i]; };
  const int U = cidx_.U, B = cidx_.B, V = cidx_.V, R = cidx_.R, I = cidx_.I;
  const int sb = cidx_.sb, sy = cidx_.sy, sv = cidx_.sv, su = cidx_.su,
            hn = cidx_.hn, hw = cidx_.hw;
  const int GU = cidx_.GU, GB = cidx_.GB, GV = cidx_.GV, GG = cidx_.GG,
            GB1 = cidx_.GB1, GB2 = cidx_.GB2, GV1 = cidx_.GV1;
  const int AF148 = cidx_.AF148, AF169 = cidx_.AF169;
  const size_t c0 = n_mag_;
  mags[c0 + 0] = m(U) - m(B);
  mags[c0 + 1] = m(B) - m(V);
  mags[c0 + 2] = m(V) - m(R);
  mags[c0 + 3] = m(V) - m(I);
  mags[c0 + 4] = m(R) - m(I);
  mags[c0 + 5] = m(sb) - m(sy);
  mags[c0 + 6] = m(su) - m(sb);
  mags[c0 + 7] = m(su) - 2 * m(sv) + m(sb);
  mags[c0 + 8] = m(sv) - 2 * m(sb) + m(sy);
  double Hbeta = m(hn) - m(hw);
  mags[c0 + 9] = 2.525 + 1.305 * Hbeta;
  mags[c0 + 10] = 2.506 + 1.368 * Hbeta;
  double gUmB = m(GU) - m(GB), gB1mB = m(GB1) - m(GB), gB2mB = m(GB2) - m(GB),
         gV1mB = m(GV1) - m(GB), gVmB = m(GV) - m(GB), gGmB = m(GG) - m(GB);
  mags[c0 + 11] = gUmB;
  mags[c0 + 12] = gB1mB;
  mags[c0 + 13] = gB2mB;
  mags[c0 + 14] = gV1mB;
  mags[c0 + 15] = gVmB;
  mags[c0 + 16] = gGmB;
  mags[c0 + 17] = 1.3764 * gUmB - 1.2162 * gB1mB - 0.8498 * gB2mB -
                  0.1554 * gV1mB + 0.8450 * gGmB + 0.3788;
  mags[c0 + 18] = 0.3235 * gUmB - 2.3228 * gB1mB + 2.3363 * gB2mB +
                  0.7495 * gV1mB - 1.0865 * gGmB - 0.8288;
  mags[c0 + 19] = 0.0255 * gUmB - 0.1740 * gB1mB + 0.4696 * gB2mB -
                  1.1205 * gV1mB + 0.7994 * gGmB - 0.4572;
  mags[c0 + 20] = m(AF148) - m(AF169);
}

}  // namespace sed
