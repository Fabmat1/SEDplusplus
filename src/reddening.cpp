#include "reddening.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "http.hpp"

namespace sed {

namespace {

// Extract the text between <tag> and </tag>, strip whitespace and the "(mag)"
// unit suffix, and parse as a double. Returns false if the tag is absent.
// Mirrors query_reddening.sl, which takes the line after the tag and strips
// "(mag)" plus whitespace.
bool extract_mag(const std::string& xml, const std::string& tag, double& out) {
  const std::string open = "<" + tag + ">";
  const std::string close = "</" + tag + ">";
  auto a = xml.find(open);
  if (a == std::string::npos) return false;
  a += open.size();
  auto b = xml.find(close, a);
  if (b == std::string::npos) return false;
  std::string inner = xml.substr(a, b - a);
  // Drop the "(mag)" unit and any surrounding whitespace, then atof.
  std::string cleaned;
  cleaned.reserve(inner.size());
  for (size_t i = 0; i < inner.size();) {
    if (inner.compare(i, 5, "(mag)") == 0) {
      i += 5;
      continue;
    }
    char c = inner[i];
    if (c != ' ' && c != '\t' && c != '\n' && c != '\r') cleaned += c;
    ++i;
  }
  if (cleaned.empty()) return false;
  out = std::atof(cleaned.c_str());
  return true;
}

}  // namespace

Reddening query_reddening(double ra, double dec) {
  char url[256];
  std::snprintf(url, sizeof(url),
                "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=%f+"
                "%f+equ+j2000",
                ra, dec);
  auto resp = sed::http::get(url);
  if (!resp.ok())
    throw std::runtime_error("query_reddening: IRSA DUST HTTP " +
                             std::to_string(resp.status));
  // The service reports success via <results status="ok"> (or an error status).
  auto st = resp.body.find("status=");
  if (st == std::string::npos ||
      resp.body.find("\"ok\"", st) == std::string::npos)
    throw std::runtime_error(
        "query_reddening: IRSA DUST did not report status ok");

  Reddening r;
  double v;
  if (extract_mag(resp.body, "meanValueSFD", v)) r.meanSFD = v;
  if (extract_mag(resp.body, "stdSFD", v)) r.stdSFD = v;
  if (extract_mag(resp.body, "meanValueSandF", v)) r.meanSandF = v;
  if (extract_mag(resp.body, "stdSandF", v)) r.stdSandF = v;
  return r;
}

}  // namespace sed
