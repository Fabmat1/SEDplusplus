#include "spectra.hpp"

#include <zlib.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "fitsio_util.hpp"
#include "http.hpp"
#include "query.hpp"
#include "texval.hpp"
#include "util.hpp"
#include "votable.hpp"

namespace sed {

namespace {

// AB zero-point constant: 2.5*log10(Const_c*1e8*10^(48.6/-2.5)); the box
// magnitude is -2.5*log10(int*norm) + kIntRef (photometric_magnitudes.sl:813).
constexpr double kIntRef = -2.407948242680184;

// ---------------------------------------------------------------------------
// Box magnitude core (photometric_magnitudes box path / magnitude_from_box).
// ---------------------------------------------------------------------------

// The strict-interior grid: u = union(wmin, wave[wmin<w<wmax], wmax) plus the
// indices of the interior points (aligned with u[1..n-2]).
struct BoxGrid {
  dvec u;
  std::vector<size_t> ind;  // indices into wave of the interior points
};

BoxGrid box_grid(const dvec& wave, double wmin, double wmax) {
  BoxGrid g;
  g.u.push_back(wmin);
  for (size_t i = 0; i < wave.size(); ++i)
    if (wmin < wave[i] && wave[i] < wmax) {
      g.ind.push_back(i);
      g.u.push_back(wave[i]);
    }
  g.u.push_back(wmax);
  return g;
}

// AB magnitude in a top-hat box [wmin, wmax]. NaN when the box is not fully
// inside the spectrum's coverage (photometric_magnitudes.sl:756 requires
// l_filter[0] >= l[0] and l_filter[-1] <= l[-1]) or no point falls in the
// box -- that box yields no row.
//
// Faithful to photometric_magnitudes.sl:768-813 including its padding quirk:
// ind = where(wmin <= l <= wmax) is padded by one point on the LEFT only --
// the right-pad condition `ind[-1] < length(ind)-1` compares the last index
// against the length of ind itself and thus can never be true (indices are
// strictly increasing from >= 0). Consequently the flux resample
// `interpol(u, l[ind], f[ind])` linearly EXTRAPOLATES at u == wmax from the
// last two in-box points instead of interpolating across the edge. This is
// observable in the fixture (e.g. LWP07276LL 2000_2500 = 10.319, not the
// 10.326 an edge-interpolating integral gives).
double box_magnitude(const dvec& wave, const dvec& flux, double wmin,
                     double wmax) {
  if (wave.size() < 2) return std::nan("");
  if (!(wmin >= wave.front() && wmax <= wave.back())) return std::nan("");
  size_t lo = wave.size(), hi = 0;
  dvec inbox;
  for (size_t i = 0; i < wave.size(); ++i)
    if (wmin <= wave[i] && wave[i] <= wmax) {
      if (i < lo) lo = i;
      hi = i;
      inbox.push_back(wave[i]);
    }
  if (inbox.empty()) return std::nan("");
  dvec u = sorted_union({wmin, wmax}, inbox);
  if (lo > 0) --lo;  // left pad (photometric_magnitudes.sl:772-773)
  dvec xp(wave.begin() + lo, wave.begin() + hi + 1);
  dvec yp(flux.begin() + lo, flux.begin() + hi + 1);
  if (xp.size() < 2) return std::nan("");  // interpol needs 2 points
  dvec fu = interp_linear(u, xp, yp);
  // fu_filter = interpol(u, [wmin,wmax], [norm,norm]) == norm exactly;
  // integrand is fu*fu_filter per point, integrated by trapezoid.
  const double norm = 1.0 / (1.0 / wmin - 1.0 / wmax);
  for (double& v : fu) v *= norm;
  double integ = integrate_trapez(u, fu);
  return -2.5 * std::log10(integ) + kIntRef;
}

// Gaussian error propagation of the trapezoid integral:
// int_error = sqrt(sum((0.5*(u[k+2]-u[k]) * sigma[ind])^2)) over the interior
// points; returned as |(-2.5/ln10)/int * int_error| (unrounded).
double box_mag_error(const dvec& wave, const dvec& flux, const dvec& sigma,
                     double wmin, double wmax) {
  BoxGrid g = box_grid(wave, wmin, wmax);
  dvec fu = interp_linear(g.u, wave, flux);
  double integ = integrate_trapez(g.u, fu);
  double s2 = 0.0;
  for (size_t k = 0; k < g.ind.size(); ++k) {
    double t = 0.5 * (g.u[k + 2] - g.u[k]) * sigma[g.ind[k]];
    s2 += t * t;
  }
  return std::abs(-2.5 / std::log(10.0) / integ * std::sqrt(s2));
}

bool parse_box(const std::string& box, double& wmin, double& wmax) {
  return std::sscanf(box.c_str(), "%lf_%lf", &wmin, &wmax) == 2;
}

// ---------------------------------------------------------------------------
// average_boxes (query_photometry.sl:354-386)
// ---------------------------------------------------------------------------

// ISIS `moment` (isis/src/math.c, mean_stddev_doubles): Welford's online
// mean/variance. Replicated exactly -- the FP value of `ave` feeds round2 and
// can sit on a rounding boundary (e.g. the fixture's 10.0565 -> 10.057).
void welford(const dvec& x, double& ave, double& sdev) {
  double mean_i = 0.0, variance_i = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    double diff = x[i] - mean_i;
    mean_i += diff / double(i + 1);
    variance_i += diff * (x[i] - mean_i);
  }
  ave = mean_i;
  sdev = x.size() > 1 ? std::sqrt(variance_i / double(x.size() - 1)) : 0.0;
}

// Collapse multiple flag-0 rows per box into one "Average" row; the
// contributing rows are kept at flag 2.
void average_boxes(std::vector<PhotEntry>& entries,
                   const std::vector<std::string>& boxes, double thres,
                   double mag_err_sys) {
  for (const auto& box : boxes) {
    std::vector<size_t> ind;
    for (size_t i = 0; i < entries.size(); ++i)
      if (entries[i].passband == box && entries[i].flag == 0) ind.push_back(i);
    const size_t ndata = ind.size();
    if (ndata <= 1) continue;
    for (size_t i : ind) entries[i].flag = 2;
    // array_sort is stable; sort the row indices by magnitude ascending.
    std::vector<size_t> order = ind;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      return entries[a].magnitude < entries[b].magnitude;
    });
    // idx_cut = [0 : max(nint(ndata*thres), 1)] (inclusive); nint is
    // round-half-away-from-zero (verified against the ISIS binary). For
    // thres < 0.5 this never exceeds ndata; clamp defensively anyway.
    size_t ncut =
        size_t(std::max<long>(std::lround(double(ndata) * thres), 1)) + 1;
    if (ncut > ndata) ncut = ndata;
    dvec mags(ncut);
    for (size_t k = 0; k < ncut; ++k) mags[k] = entries[order[k]].magnitude;
    double ave, sdev;
    welford(mags, ave, sdev);
    double mag_error = std::sqrt(sdev * sdev + mag_err_sys * mag_err_sys);
    RoundedErr re = round_err(mag_error);
    double dist = 0.0;
    for (size_t k = 0; k < ncut; ++k)
      dist += entries[order[k]].angu_dist_arcsec;
    dist /= double(ncut);
    // ISIS table.add_entry has no duplicate resolution -- plain append.
    entries.push_back({0, "box", box, round2(ave, int(re.digit) - 1), re.value,
                       "magnitude", dist, "Average"});
  }
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

// Write `data` as a gzip file (stands in for the S-Lang's wget + `gzip -f`).
bool write_gzip(const std::string& path, const std::string& data) {
  gzFile g = gzopen(path.c_str(), "wb");
  if (!g) return false;
  bool ok = data.empty() ||
            gzwrite(g, data.data(), unsigned(data.size())) == int(data.size());
  return gzclose(g) == Z_OK && ok;
}

// Download `url` into <dir>/<filename>.gz unless that file already exists
// (the S-Lang caches the same way: only wget when stat_file(...)==NULL). A
// failed transfer still leaves a (bogus) file, like `wget -q -O` does; the
// FITS read then fails and the spectrum is skipped.
void download_cached(const std::string& dir, const std::string& filename,
                     const std::string& url) {
  std::filesystem::create_directories(dir);
  const std::string gz = dir + "/" + filename + ".gz";
  if (std::filesystem::exists(gz)) return;
  std::string body;
  try {
    http::Options opt;
    opt.timeout_s = 120;
    auto resp = http::get(url, opt);
    body = resp.body;  // non-2xx bodies are written too (wget -O semantics)
  } catch (const std::exception&) {
    body.clear();
  }
  write_gzip(gz, body);
}

std::string strip_suffix(std::string s, const std::string& suf) {
  if (s.size() >= suf.size() &&
      s.compare(s.size() - suf.size(), suf.size(), suf) == 0)
    s.erase(s.size() - suf.size());
  return s;
}

// One CSV line -> cells (RFC-4180-ish: quotes, doubled quotes, no newlines).
std::vector<std::string> csv_split(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < line.size(); ++i) {
    char ch = line[i];
    if (inq) {
      if (ch == '"') {
        if (i + 1 < line.size() && line[i + 1] == '"') {
          cur += '"';
          ++i;
        } else {
          inq = false;
        }
      } else {
        cur += ch;
      }
    } else if (ch == '"') {
      inq = true;
    } else if (ch == ',') {
      out.push_back(cur);
      cur.clear();
    } else if (ch != '\r') {
      cur += ch;
    }
  }
  out.push_back(cur);
  return out;
}

std::string lower(std::string s) {
  for (char& c : s) c = char(std::tolower((unsigned char)c));
  return s;
}

// ---------------------------------------------------------------------------
// IUE (query_photometry.sl:1862-1941)
// ---------------------------------------------------------------------------

void iue_boxes(double ra, double dec, const SpectraOptions& opt,
               std::vector<PhotEntry>& entries,
               std::vector<std::string>& boxes) {
  double radius_iue = std::max(opt.search_radius_deg, 70.0 / 3600.0);
  if (opt.force_search_radius) radius_iue = opt.search_radius_deg;
  const double radius_arcsec = radius_iue * 3600.0;

  votable::Document doc;
  try {
    doc = vizier_fetch("VI/110/inescat", ra, dec, radius_arcsec);
  } catch (const std::exception&) {
  }
  const votable::Table* t = doc.find("VI/110/inescat");
  auto col = [&](const char* name) {
    int i = t ? t->field_index(name) : -1;
    if (i < 0 && t) i = t->field_index_ci(name);
    return i;
  };
  std::vector<size_t> keep;
  int c_cam = -1, c_img = -1, c_disp = -1, c_aper = -1, c_abn = -1, c_ra = -1,
      c_de = -1;
  if (t) {
    c_cam = col("Camera");
    c_img = col("Image");
    c_disp = col("Disp");
    c_aper = col("Aper");
    c_abn = col("AbnCode");
    c_ra = col("RAJ2000");
    c_de = col("DEJ2000");
    if (c_cam >= 0 && c_img >= 0 && c_disp >= 0 && c_aper >= 0 && c_abn >= 0 &&
        c_ra >= 0 && c_de >= 0) {
      // the VizieR parser keeps rows within the flat-distance radius
      const double d2r = 3.14159265358979323846 / 180.0;
      for (size_t r = 0; r < t->rows.size(); ++r) {
        double rr = std::atof(t->rows[r][c_ra].c_str());
        double dd = std::atof(t->rows[r][c_de].c_str());
        double dra = (rr - ra) * std::cos(dec * d2r);
        double dde = dd - dec;
        if (3600.0 * std::sqrt(dra * dra + dde * dde) <= radius_arcsec)
          keep.push_back(r);
      }
    }
  }
  if (keep.empty()) {
    std::fprintf(stderr,
                 "Information from 'query_photometry': Queried catalogue "
                 "%s: Found nothing.\n",
                 "VI/110/inescat");
    return;
  }
  std::fprintf(stderr,
               "Information from 'query_photometry': Queried catalogue %s: "
               "Found %zu spectra.\n",
               "VI/110/inescat", keep.size());

  const std::string dir = opt.cache_dir + "/IUE";
  for (size_t r : keep) {
    const auto& row = t->rows[r];
    int flag = 0;
    const std::string abn = trim(row[c_abn]);
    const std::string aper = trim(row[c_aper]);
    if (abn.find_first_of("ABCRTZ") != std::string::npos || aper == "S")
      flag = -1;  // http://cdsarc.u-strasbg.fr/viz-bin/ReadMe/VI/110
    std::string disp = trim(row[c_disp]);
    for (char& ch : disp)
      if (ch == 'H') ch = 'R';
    char fn[128];
    std::snprintf(fn, sizeof fn, "%s%05.0f%s%s.FITS", trim(row[c_cam]).c_str(),
                  std::atof(row[c_img].c_str()), disp.c_str(), aper.c_str());
    const std::string filename = fn;
    std::fprintf(stderr,
                 "Information from 'query_photometry': Downloading IUE "
                 "spectrum %s\n",
                 filename.c_str());
    // https does not work anymore (05/2022) -- same endpoint as the S-Lang
    download_cached(dir, filename,
                    "http://sdc.cab.inta-csic.es/ines/jsp/"
                    "SingleDownload.jsp?filename=" +
                        filename);
    const std::string gzpath = dir + "/" + filename + ".gz";
    dvec wave, flux, sigma, quality;
    try {
      wave = fits_read_double_column(gzpath + "[1]", "WAVELENGTH");
      flux = fits_read_double_column(gzpath + "[1]", "FLUX");
      sigma = fits_read_double_column(gzpath + "[1]", "SIGMA");
      quality = fits_read_double_column(gzpath + "[1]", "QUALITY");
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "Information from 'query_photometry': Cannot read IUE "
                   "spectrum %s.gz\n",
                   filename.c_str());
      continue;  // some IUE files are corrupt -> skip them
    }
    double angu_dist_arcsec =
        3600.0 * angular_separation_deg(ra, dec, std::atof(row[c_ra].c_str()),
                                        std::atof(row[c_de].c_str()));
    // remove bad pixels (flux>0 and quality==0),
    // https://archive.stsci.edu/iue/manual/dacguide/node60.html
    dvec w2, f2, s2;
    for (size_t i = 0; i < wave.size(); ++i)
      if (flux[i] > 0 && quality[i] == 0) {
        w2.push_back(wave[i]);
        f2.push_back(flux[i]);
        s2.push_back(sigma[i]);
      }
    boxes = {"1300_1800", "2000_2500", "2500_3000"};
    if (w2.empty() || f2.empty()) {
      boxes.clear();  // S-Lang: boxes = String_Type[0] (persists!)
      continue;
    }
    for (const auto& box : boxes) {
      double wmin, wmax;
      if (!parse_box(box, wmin, wmax)) continue;
      double mag = box_magnitude(w2, f2, wmin, wmax);
      if (std::isnan(mag)) continue;
      double mag_error = box_mag_error(w2, f2, s2, wmin, wmax);
      // add generic uncertainty to IUE
      mag_error = std::sqrt(mag_error * mag_error + 0.04 * 0.04);
      RoundedErr re = round_err(mag_error);
      entries.push_back({flag, "box", box,
                         round2(mag, int(re.digit) - 1), re.value, "magnitude",
                         angu_dist_arcsec,
                         strip_suffix(filename + ".gz", ".FITS.gz") +
                             "_VI/110/inescat"});
    }
  }
  // Replace magnitudes from individual spectra by a single measurement
  const double thres_iue = 0.32;  // take brightest "thres_iue" fraction
  average_boxes(entries, boxes, thres_iue, /*mag_err_sys=*/0.03);
}

// ---------------------------------------------------------------------------
// MAST FUSE/FOS (query_photometry.sl:1942-2052) via STILTS SSAP cone search
// ---------------------------------------------------------------------------

// Run the exact STILTS command ssapquery.sl uses and return the CSV lines.
std::vector<std::string> stilts_ssap(double ra, double dec, double radius) {
  char cmd[512];
  std::snprintf(cmd, sizeof cmd,
                "stilts cone servicetype=ssa serviceurl='%s' lon=%.5f "
                "lat=%.5f radius=%.5f ofmt=csv",
                "https://archive.stsci.edu/ssap/search2.php", ra, dec, radius);
  std::vector<std::string> lines;
  FILE* p = ::popen(cmd, "r");
  if (!p) return lines;
  std::string cur;
  char buf[4096];
  while (std::fgets(buf, sizeof buf, p)) {
    cur += buf;
    size_t pos;
    while ((pos = cur.find('\n')) != std::string::npos) {
      lines.push_back(cur.substr(0, pos));
      cur.erase(0, pos + 1);
    }
  }
  if (!cur.empty()) lines.push_back(cur);
  ::pclose(p);
  return lines;
}

void mast_boxes(double ra, double dec, const SpectraOptions& opt,
                std::vector<PhotEntry>& entries,
                std::vector<std::string>& boxes) {
  double radius_mast = std::max(opt.search_radius_deg, 60.0 / 3600.0);
  if (opt.force_search_radius) radius_mast = opt.search_radius_deg;

  auto lines = stilts_ssap(ra, dec, radius_mast);
  auto nothing = [] {
    std::fprintf(stderr,
                 "Information from 'query_photometry': Queried catalogue "
                 "MAST: Found nothing.\n");
  };
  if (lines.size() <= 1) {
    nothing();
    return;
  }
  auto header = csv_split(lines[0]);
  auto idx = [&](const char* name) {
    for (size_t i = 0; i < header.size(); ++i)
      if (header[i] == name) return int(i);
    return -1;
  };
  const int c_url = idx("url"), c_tel = idx("telescop"),
            c_ins = idx("instrume"), c_aper = idx("aperture"),
            c_cri = idx("cr_ident"), c_cal = idx("fluxcal"),
            c_snr = idx("der_snr"), c_sep = idx("ang_sep"),
            c_sax = idx("spectralaxisname"), c_fax = idx("fluxaxisname");
  if (c_url < 0 || c_tel < 0 || c_ins < 0 || c_aper < 0 || c_cri < 0 ||
      c_cal < 0 || c_snr < 0 || c_sep < 0 || c_sax < 0 || c_fax < 0) {
    nothing();
    return;
  }
  std::vector<std::vector<std::string>> rows;
  for (size_t i = 1; i < lines.size(); ++i) {
    if (trim(lines[i]).empty()) continue;
    auto cells = csv_split(lines[i]);
    if (cells.size() < header.size()) continue;
    // keep FUSE or FOS ...
    if (cells[c_tel] != "FUSE" && cells[c_ins] != "FOS") continue;
    // ... with der_snr >= 2. NOTE: the S-Lang also *intends* to require
    // object != "WAVE" and aperture != "RFPT", but writes
    // `string(query.object) != "WAVE"` -- string() of an *array* yields
    // "String_Type[n]", so both conditions are scalar-true no-ops. Replicated
    // faithfully by NOT filtering on object/aperture here (this is why FOS
    // wavelength-calibration exposures get downloaded at all).
    if (!(std::atof(cells[c_snr].c_str()) >= 2.0)) continue;
    rows.push_back(std::move(cells));
  }
  if (rows.empty()) {
    nothing();
    return;
  }
  std::fprintf(stderr,
               "Information from 'query_photometry': Queried catalogue MAST: "
               "Found %zu spectra.\n",
               rows.size());

  const std::string dir = opt.cache_dir + "/MAST";
  for (const auto& row : rows) {
    int flag = 0;
    const std::string& fluxcal = row[c_cal];
    const std::string& aperture = row[c_aper];
    // string_match(fluxcal,"ABSOLUTE")!=1 == "does not match at position 1"
    if (fluxcal.rfind("ABSOLUTE", 0) != 0 ||
        aperture == "1.25x20" ||               // FUSE-HRS
        aperture == "4x20" ||                  // FUSE-MRS
        aperture.rfind("unknown", 0) == 0) {
      flag = -1;
    }
    std::string filename = row[c_tel] + "-" + row[c_ins] + "-" + aperture +
                           "_" + row[c_cri] + ".FITS";
    filename.erase(std::remove(filename.begin(), filename.end(), ' '),
                   filename.end());
    std::fprintf(stderr,
                 "Information from 'query_photometry': Downloading %s "
                 "spectrum %s\n",
                 row[c_ins].c_str(), filename.c_str());
    download_cached(dir, filename, row[c_url]);
    const std::string gzpath = dir + "/" + filename + ".gz";

    std::string filetype;
    dvec wave, flux, sigma;
    try {
      filetype = fits_read_key_string(gzpath, "FILETYPE");
      wave = fits_read_double_column(gzpath + "[1]", lower(row[c_sax]));
      flux = fits_read_double_column(gzpath + "[1]", lower(row[c_fax]));
      sigma = fits_read_double_column(gzpath + "[1]", "sigma");
    } catch (const std::exception&) {
      std::fprintf(stderr,
                   "Information from 'query_photometry': Cannot read %s "
                   "spectrum %s.gz\n",
                   row[c_ins].c_str(), filename.c_str());
      continue;  // some files are corrupt -> skip them
    }
    if (filetype == "WAV") continue;  // exclude calibration exposures
    double angu_dist_arcsec = 3600.0 * std::atof(row[c_sep].c_str());
    // remove bad pixels and EUV
    dvec w2, f2, s2;
    for (size_t i = 0; i < wave.size(); ++i)
      if (flux[i] > 0 && sigma[i] > 0 && wave[i] > 920.0) {
        w2.push_back(wave[i]);
        f2.push_back(flux[i]);
        s2.push_back(sigma[i]);
      }
    if (row[c_tel] == "FUSE") {
      boxes = {"950_1000", "1000_1080", "1091_1120"};
      flag = -1;  // always flag FUSE since the pointing is not reliable
    }
    if (row[c_ins] == "FOS") {
      // dynamic 4-point grid [wmin+15 : wmax-15 : #4] -> 3 boxes
      if (w2.empty()) continue;  // (unreachable in ISIS; would throw there)
      double wmin = *std::min_element(w2.begin(), w2.end());
      double wmax = *std::max_element(w2.begin(), w2.end());
      double lo = wmin + 15.0, hi = wmax - 15.0;
      boxes.clear();
      for (int k = 0; k < 3; ++k) {
        char b[64];
        std::snprintf(b, sizeof b, "%.0f_%.0f", lo + (hi - lo) * k / 3.0,
                      lo + (hi - lo) * (k + 1) / 3.0);
        boxes.push_back(b);
      }
    }
    if (w2.size() < 2) continue;  // guard (ISIS would throw on empty input)
    for (const auto& box : boxes) {
      double wmin, wmax;
      if (!parse_box(box, wmin, wmax)) continue;
      double mag = box_magnitude(w2, f2, wmin, wmax);
      if (std::isnan(mag)) continue;
      RoundedErr re = round_err(box_mag_error(w2, f2, s2, wmin, wmax));
      entries.push_back({flag, "box", box, round2(mag, int(re.digit) - 1),
                         re.value, "magnitude", angu_dist_arcsec,
                         strip_suffix(filename + ".gz", ".FITS.gz") + "_MAST"});
    }
  }
  // Replace magnitudes from individual spectra by a single measurement.
  // NOTE the latent S-Lang quirk: `boxes` here is whatever the *last processed
  // spectrum* left behind (FUSE list, FOS dynamic list, or even the IUE list
  // if no MAST spectrum was processed) -- replicated via the shared `boxes`.
  const double thres_mast = 0.16;  // take brightest "thres_mast" fraction
  average_boxes(entries, boxes, thres_mast, /*mag_err_sys=*/0.005);
}

// ---------------------------------------------------------------------------
// Gaia DR3 XP-sampled (query_gaia_spectrum + box_from_gaia, lines 1-275)
// ---------------------------------------------------------------------------
// Faithful to *current ISIS behaviour*, which in practice always yields zero
// rows (PHASE2_SPECTRA.md §0): the mirror probe prefers gaia.aip.de and then
// bails, and even the esac branch uses the broken `<tap>data?...` datalink URL
// (404). Deliberately NOT fixed to the working data-server endpoint -- parity
// with the ISIS fixtures rules; see PHASE2_SPECTRA.md for the opt-in idea.

struct GaiaSpectrum {
  dvec wavelength, flux, flux_error;
};

std::optional<GaiaSpectrum> query_gaia_spectrum(double ra, double dec,
                                                const std::string& cache_dir) {
  static const char* kUrls[] = {"https://gaia.aip.de/tap",
                                "https://gea.esac.esa.int/tap-server/tap"};
  std::string url_tap;
  for (const char* u : kUrls) {
    // S-Lang probe: `curl -s --connect-timeout 4 -m 1 <url>` exit code == 0
    if (http::probe(u, 4, 1)) {
      url_tap = u;
      break;
    }
    std::fprintf(stderr,
                 "Information from 'query_gaia_spectrum': URL offline: %s\n",
                 u);
  }
  if (url_tap.empty()) return std::nullopt;

  // Resolve the DR3 source_id from (ra, dec) on the chosen mirror.
  char adql[512];
  std::snprintf(adql, sizeof adql,
                "SELECT designation,ra,dec FROM gaiadr3.gaia_source WHERE "
                "1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE ('ICRS', %f, %f, "
                "%f ))",
                ra, dec, 0.01);
  // FORMAT=votable/td: both mirrors are/emulate DaCHS-style services that
  // default to BINARY serialization (same deviation as astrometry.cpp).
  std::string body =
      "REQUEST=doQuery&LANG=ADQL&FORMAT=votable%2Ftd&QUERY=" +
      http::url_encode(adql);
  std::string id;
  try {
    http::Options hopt;
    hopt.timeout_s = 120;
    auto resp = http::post(url_tap + "/sync", body,
                           "application/x-www-form-urlencoded", hopt);
    if (!resp.ok()) return std::nullopt;
    auto doc = votable::parse(resp.body);
    const votable::Table* t =
        doc.tables.empty() ? nullptr : &doc.tables.front();
    if (!t || t->rows.empty()) return std::nullopt;
    int c_des = t->field_index_ci("designation");
    int c_ra = t->field_index_ci("ra");
    int c_de = t->field_index_ci("dec");
    if (c_des < 0 || c_ra < 0 || c_de < 0) return std::nullopt;
    double best = 1e300;
    int besti = -1;
    for (size_t r = 0; r < t->rows.size(); ++r) {
      double a = angular_separation_deg(ra, dec,
                                        std::atof(t->rows[r][c_ra].c_str()),
                                        std::atof(t->rows[r][c_de].c_str()));
      if (a < best) {
        best = a;
        besti = int(r);
      }
    }
    if (besti < 0) return std::nullopt;
    std::string designation = t->rows[besti][c_des];  // "Gaia DR3 <id>"
    auto tok = split_ws(designation);
    if (tok.empty()) return std::nullopt;
    id = tok.back();
  } catch (const std::exception&) {
    std::fprintf(stderr,
                 "Warning in 'query_gaia_spectrum': Failed to obtain XP "
                 "spectrum with TAP.\n");
    return std::nullopt;
  }

  // ISIS guard (query_photometry.sl:121): only proceed when the answering
  // mirror is esac -- with AIP normally reachable this returns NULL here.
  if (url_tap.find("gea.esac.esa.int") == std::string::npos)
    return std::nullopt;

  // ISIS datalink URL (line 125): url_tap + "data?..." =
  // ".../tap-server/tapdata?..." which 404s; kept verbatim for parity.
  const std::string query_str = url_tap + "data?ID=Gaia+DR3+" + id +
                                "&RETRIEVAL_TYPE=XP_SAMPLED&format=FITS";
  const std::string fits_name = cache_dir + "/XP_SAMPLED.fits";
  std::string dl;
  try {
    http::Options hopt;
    hopt.timeout_s = 120;
    dl = http::get(query_str, hopt).body;  // curl -s ... > XP_SAMPLED.fits
  } catch (const std::exception&) {
  }
  {
    std::FILE* f = std::fopen(fits_name.c_str(), "wb");
    if (f) {
      if (!dl.empty()) std::fwrite(dl.data(), 1, dl.size(), f);
      std::fclose(f);
    }
  }
  GaiaSpectrum s;
  try {
    s.wavelength = fits_read_double_column(fits_name + "[1]", "wavelength");
    s.flux = fits_read_double_column(fits_name + "[1]", "flux");
    s.flux_error = fits_read_double_column(fits_name + "[1]", "flux_error");
  } catch (const std::exception&) {
    std::remove(fits_name.c_str());  // FitsError branch: s = NULL
    return std::nullopt;
  }
  std::remove(fits_name.c_str());
  // nm -> Angstrom; W m^-2 nm^-1 -> erg cm^-2 s^-1 A^-1; +30% blanket error
  for (double& w : s.wavelength) w *= 10.0;
  for (double& f : s.flux) f *= 100.0;
  for (double& e : s.flux_error) e *= 100.0 * 1.30;
  // u-band error bumps: +0.5/0.6/0.6 of the median positive error in each
  // sub-range, applied in sequence (<3800, <3660, <3500).
  const double cuts[3] = {3800.0, 3660.0, 3500.0};
  const double fracs[3] = {0.5, 0.6, 0.6};
  for (int k = 0; k < 3; ++k) {
    dvec good;
    std::vector<size_t> idx_u;
    for (size_t i = 0; i < s.wavelength.size(); ++i)
      if (s.wavelength[i] < cuts[k]) {
        idx_u.push_back(i);
        if (s.flux_error[i] > 0) good.push_back(s.flux_error[i]);
      }
    if (idx_u.empty() || good.empty()) continue;
    double med = median(good);
    for (size_t i : idx_u) s.flux_error[i] += fracs[k] * med;
  }
  return s;
}

void gaia_xp_boxes(double ra, double dec, const SpectraOptions& opt,
                   std::vector<PhotEntry>& entries) {
  auto sg = query_gaia_spectrum(ra, dec, opt.cache_dir);
  if (!sg) return;  // box_from_gaia(table, NULL) is a no-op
  dvec w2, f2, e2;
  for (size_t i = 0; i < sg->wavelength.size(); ++i)
    if (sg->flux[i] > 0 && sg->flux_error[i] > 0) {
      w2.push_back(sg->wavelength[i]);
      f2.push_back(sg->flux[i]);
      e2.push_back(sg->flux_error[i]);
    }
  if (w2.empty()) return;
  static const char* kBoxes[] = {
      "3400_3660", "3660_4000", "4000_4300", "4300_4600", "4600_4740",
      "4740_5200", "5200_5600", "5600_6000", "6000_6400", "6400_6800",
      "6800_8400", "8400_8800", "8800_9200", "9200_10200"};
  for (const char* box : kBoxes) {
    double wmin, wmax;
    if (!parse_box(box, wmin, wmax)) continue;
    double mag = box_magnitude(w2, f2, wmin, wmax);
    if (std::isnan(mag)) continue;
    RoundedErr re = round_err(box_mag_error(w2, f2, e2, wmin, wmax));
    int flag = wmax < 3700.0 ? -1 : 0;
    entries.push_back({flag, "box", box, round2(mag, int(re.digit) - 1),
                       re.value, "magnitude", 0.0, "Gaia_DR3_XP"});
  }
}

}  // namespace

void add_spectra_boxes(double ra, double dec, const SpectraOptions& opt,
                       std::vector<PhotEntry>& entries) {
  // Gaia XP is attempted unconditionally (query_photometry.sl:1854-1859).
  gaia_xp_boxes(ra, dec, opt, entries);
  // `boxes` is one shared function-local in the S-Lang; the MAST averaging
  // step sees whatever the last processed spectrum (IUE or MAST) set.
  std::vector<std::string> boxes;
  if (opt.iue) iue_boxes(ra, dec, opt, entries, boxes);
  if (opt.mast) mast_boxes(ra, dec, opt, entries, boxes);
}

}  // namespace sed
