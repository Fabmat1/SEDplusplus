#include "maggrid.hpp"

#include <fitsio.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

#include "extinction.hpp"
#include "fitsio_util.hpp"
#include "passbands.hpp"
#include "synthmag.hpp"

namespace sed {

namespace {

namespace fs = std::filesystem;

std::uint64_t fnv1a_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::uint64_t h = 1469598103934665603ULL;
  char buf[1 << 16];
  while (in) {
    in.read(buf, sizeof buf);
    for (std::streamsize i = 0; i < in.gcount(); ++i)
      h = (h ^ std::uint8_t(buf[i])) * 1099511628211ULL;
  }
  return h;
}

std::string mag_dir(const std::string& griddir) {
  std::string g = griddir;
  if (g.empty() || g.back() != '/') g += '/';
  return g + "mag/";
}

}  // namespace

std::uint64_t mag_filters_hash(const std::string& refdata) {
  static std::map<std::string, std::uint64_t> cache;
  auto it = cache.find(refdata);
  if (it != cache.end()) return it->second;
  std::string path = refdata + "/filter_passbands.fits.gz";
  if (!fs::exists(path)) path = refdata + "/filter_passbands.fits";
  std::uint64_t h = fnv1a_file(path);
  h = (h ^ std::uint64_t(MAG_GRID_VERSION)) * 1099511628211ULL;
  cache.emplace(refdata, h);
  return h;
}

std::vector<MagBandInfo> read_mag_bands(const std::string& bands_fits) {
  fitsfile* fp = nullptr;
  int status = 0;
  if (fits_open_file(&fp, (bands_fits + "[1]").c_str(), READONLY, &status))
    throw std::runtime_error("cannot open " + bands_fits);
  long nrows = 0;
  fits_get_num_rows(fp, &nrows, &status);
  auto colnum = [&](const char* name) {
    int c = 0;
    fits_get_colnum(fp, CASEINSEN, const_cast<char*>(name), &c, &status);
    return c;
  };
  const int c_sys = colnum("SYSTEM"), c_pb = colnum("PASSBAND"),
            c_hb = colnum("HBETA"), c_le = colnum("LEFF");
  std::vector<MagBandInfo> bands(nrows);
  std::vector<int> hb(nrows);
  dvec le(nrows);
  {
    char val[128];
    char* valp = val;
    for (long r = 1; r <= nrows && !status; ++r) {
      int anynul = 0;
      fits_read_col(fp, TSTRING, c_sys, r, 1, 1, nullptr, &valp, &anynul,
                    &status);
      bands[r - 1].system = val;
      fits_read_col(fp, TSTRING, c_pb, r, 1, 1, nullptr, &valp, &anynul,
                    &status);
      bands[r - 1].passband = val;
    }
    int anynul = 0;
    if (nrows) {
      fits_read_col(fp, TINT, c_hb, 1, 1, nrows, nullptr, hb.data(), &anynul,
                    &status);
      fits_read_col(fp, TDOUBLE, c_le, 1, 1, nrows, nullptr, le.data(),
                    &anynul, &status);
    }
  }
  fits_close_file(fp, &status);
  if (status) throw std::runtime_error("error reading " + bands_fits);
  for (long r = 0; r < nrows; ++r) {
    bands[r].hbeta = hb[r] != 0;
    bands[r].leff = le[r];
  }
  return bands;
}

bool read_mag_manifest(const std::string& manifest_txt, int& version,
                       std::uint64_t& filters_hash) {
  std::ifstream in(manifest_txt);
  if (!in) return false;
  version = 0;
  filters_hash = 0;
  std::string key;
  while (in >> key) {
    if (key == "version")
      in >> version;
    else if (key == "filters")
      in >> std::hex >> filters_hash >> std::dec;
    else
      in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return version != 0;
}

// ---------------------------------------------------------------------------
// --premag: walk the ATLAS12/ node tree and write the mag/ mirror
// ---------------------------------------------------------------------------

int premag_main(int argc, char** argv) {
  std::vector<std::string> boxes;
  bool force = false;
  std::string griddir, refdata;
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--force")
      force = true;
    else if (griddir.empty())
      griddir = a;
    else if (refdata.empty())
      refdata = a;
    else
      boxes.push_back(a);
  }
  if (griddir.empty() || refdata.empty()) {
    std::fprintf(stderr,
                 "usage: sedfit --premag <griddir> <refdata> [--force] "
                 "[box ...]\n");
    return 1;
  }
  if (griddir.back() != '/') griddir += '/';
  const std::string mdir = mag_dir(griddir);
  const std::uint64_t fhash = mag_filters_hash(refdata);

  // up-to-date check: manifest version + filter hash match and the existing
  // band table already covers every requested box
  std::set<std::string> box_set(boxes.begin(), boxes.end());
  {
    int ver = 0;
    std::uint64_t h = 0;
    if (!force && read_mag_manifest(mdir + "manifest.txt", ver, h) &&
        ver == MAG_GRID_VERSION && h == fhash &&
        fs::exists(mdir + "bands.fits")) {
      std::set<std::string> have;
      for (const auto& b : read_mag_bands(mdir + "bands.fits"))
        if (b.system == "box") have.insert(b.passband);
      if (std::includes(have.begin(), have.end(), box_set.begin(),
                        box_set.end())) {
        std::printf("PREMAG UPTODATE\n");
        return 0;
      }
      // rebuild with the union so alternating campaigns do not thrash
      box_set.insert(have.begin(), have.end());
    }
  }

  boxes.assign(box_set.begin(), box_set.end());
  PassbandDB db(refdata, /*apply_ZPO_corr=*/false, boxes);
  dvec lambda =
      fits_read_double_column(griddir + "ATLAS12/lambda.fits[1]", "l");
  dvec k0, s;
  extinction_curve_components(lambda, k0, s);
  std::vector<int> rows(db.n_magnitudes());
  for (size_t i = 0; i < rows.size(); ++i) rows[i] = int(i);
  SynthMag synth(db, lambda, rows);
  const auto& filters = synth.filters();
  const size_t nb = filters.size();
  if (nb == 0) {
    std::fprintf(stderr, "premag: no filter overlaps the grid wavelengths\n");
    return 1;
  }

  // per-band premultiplied moment weights: coeff, coeff*k0, coeff*s, ...
  struct BandW {
    std::vector<int> idx;
    dvec c, ck, cs, ckk, cks, css;
  };
  std::vector<BandW> bw(nb);
  std::vector<MagBandInfo> bands(nb);
  for (size_t j = 0; j < nb; ++j) {
    const auto& fw = filters[j];
    const BandEntry& be = db.entries()[fw.row];
    const size_t n = fw.idx.size();
    BandW& w = bw[j];
    w.idx = fw.idx;
    w.c.resize(n);
    w.ck.resize(n);
    w.cs.resize(n);
    w.ckk.resize(n);
    w.cks.resize(n);
    w.css.resize(n);
    double csum = 0.0, clsum = 0.0;
    for (size_t k = 0; k < n; ++k) {
      const int i = fw.idx[k];
      const double c = fw.coeff[k];
      w.c[k] = c;
      w.ck[k] = c * k0[i];
      w.cs[k] = c * s[i];
      w.ckk[k] = c * k0[i] * k0[i];
      w.cks[k] = c * k0[i] * s[i];
      w.css[k] = c * s[i] * s[i];
      csum += c;
      clsum += c * lambda[i];
    }
    bands[j] = {be.system, be.passband, fw.hbeta,
                csum > 0 ? clsum / csum : 0.0};
  }

  // node list from the actual files (grids can be ragged)
  std::vector<std::string> nodes;  // relpath below ATLAS12/, without .fits
  {
    const fs::path root = fs::path(griddir) / "ATLAS12";
    for (const auto& e : fs::recursive_directory_iterator(root)) {
      if (!e.is_regular_file() || e.path().extension() != ".fits") continue;
      if (e.path().filename() == "lambda.fits") continue;
      std::string rel = fs::relative(e.path(), root).string();
      nodes.push_back(rel.substr(0, rel.size() - 5));
    }
    std::sort(nodes.begin(), nodes.end());
  }
  const int total = int(nodes.size());
  if (total == 0) {
    std::fprintf(stderr, "premag: no node files under %sATLAS12/\n",
                 griddir.c_str());
    return 1;
  }

  // bands.fits
  fs::create_directories(mdir);
  {
    std::vector<FitsCol> cols(4);
    cols[0] = {"SYSTEM", 'A', {}, {}, {}};
    cols[1] = {"PASSBAND", 'A', {}, {}, {}};
    cols[2] = {"HBETA", 'J', {}, {}, {}};
    cols[3] = {"LEFF", 'D', {}, {}, {}};
    for (const auto& b : bands) {
      cols[0].s.push_back(b.system);
      cols[1].s.push_back(b.passband);
      cols[2].i.push_back(b.hbeta ? 1 : 0);
      cols[3].d.push_back(b.leff);
    }
    FitsWriter w(mdir + "bands.fits");
    w.primary_hdu();
    w.binary_table("BANDS", cols);
    w.close();
  }

  // nodes
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  dvec I0(nb), K(nb), S(nb), V00(nb), V01(nb), V11(nb), MAG(nb), LEFF(nb);
  for (size_t j = 0; j < nb; ++j) LEFF[j] = bands[j].leff;
  int done = 0;
  for (const std::string& rel : nodes) {
    dvec f = fits_read_double_column(griddir + "ATLAS12/" + rel + ".fits[1]",
                                     "c");
    if (f.size() != lambda.size())
      throw std::runtime_error("flux/lambda size mismatch in " + rel);
    for (size_t j = 0; j < nb; ++j) {
      const BandW& w = bw[j];
      double i0 = 0, mk = 0, ms = 0, mkk = 0, mks = 0, mss = 0;
      for (size_t k = 0; k < w.idx.size(); ++k) {
        const double fv = f[w.idx[k]];
        i0 += w.c[k] * fv;
        mk += w.ck[k] * fv;
        ms += w.cs[k] * fv;
        mkk += w.ckk[k] * fv;
        mks += w.cks[k] * fv;
        mss += w.css[k] * fv;
      }
      if (!(i0 > 0) || !std::isfinite(i0)) {
        I0[j] = K[j] = S[j] = V00[j] = V01[j] = V11[j] = MAG[j] = NaN;
        continue;
      }
      const double kb = mk / i0, sb = ms / i0;
      I0[j] = i0;
      K[j] = kb;
      S[j] = sb;
      V00[j] = std::max(0.0, mkk / i0 - kb * kb);
      V01[j] = mks / i0 - kb * sb;
      V11[j] = std::max(0.0, mss / i0 - sb * sb);
      MAG[j] = -2.5 * std::log10(i0) + (bands[j].hbeta ? 0.0 : MAG_AB_REF);
    }
    const fs::path out = fs::path(mdir) / (rel + ".fits");
    fs::create_directories(out.parent_path());
    std::vector<FitsCol> cols = {
        {"I0", 'D', I0, {}, {}},   {"K", 'D', K, {}, {}},
        {"S", 'D', S, {}, {}},     {"V00", 'D', V00, {}, {}},
        {"V01", 'D', V01, {}, {}}, {"V11", 'D', V11, {}, {}},
        {"MAG", 'D', MAG, {}, {}}, {"LEFF", 'D', LEFF, {}, {}}};
    FitsWriter w(out.string());
    w.primary_hdu();
    w.binary_table("MAGS", cols);
    w.close();
    if (++done % 200 == 0 || done == total) {
      std::printf("PREMAG %d %d\n", done, total);
      std::fflush(stdout);
    }
  }

  // manifest last: an interrupted run never validates
  {
    std::ofstream m(mdir + "manifest.txt");
    m << "version " << MAG_GRID_VERSION << "\n"
      << "filters " << std::hex << fhash << std::dec << "\n"
      << "nbands " << nb << "\n"
      << "nnodes " << total << "\n";
  }
  std::printf("PREMAG DONE %d\n", total);
  return 0;
}

}  // namespace sed
