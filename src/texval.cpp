#include "texval.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace sed {

namespace {
double round2(double x, int digits) {
  // S-Lang round2(x, -14): round at decimal place 'digits'
  double f = std::pow(10.0, -digits);
  return std::round(x * f) / f;
}

double atof_fmt(double x, double digits) {
  char fmt[16], buf[64];
  if (digits < 0)
    std::snprintf(fmt, sizeof fmt, "%%.%df", int(std::lround(-digits)));
  else
    std::snprintf(fmt, sizeof fmt, "%%.0f");
  std::snprintf(buf, sizeof buf, fmt, x);
  return std::atof(buf);
}
}  // namespace

RoundedErr round_err(double x) {
  RoundedErr r{x, 0.0};
  if (x == 0.0) return r;
  double digits = std::floor(std::log10(x));
  if (round2(x / std::pow(10.0, digits), -14) < 3) digits -= 1;
  double v = std::round(x * std::pow(10.0, -digits + 2)) /
             std::pow(10.0, -digits + 2);
  v = std::ceil(v * std::pow(10.0, -digits)) / std::pow(10.0, -digits);
  r.value = atof_fmt(v, digits);
  r.digit = digits;
  return r;
}

RoundedErr round_err_digits(double x, double digits) {
  RoundedErr r{x, digits};
  if (x == 0.0) return r;
  double v = std::round(x * std::pow(10.0, -digits + 2)) /
             std::pow(10.0, -digits + 2);
  v = std::ceil(v * std::pow(10.0, -digits)) / std::pow(10.0, -digits);
  r.value = atof_fmt(v, digits);
  r.digit = digits;
  return r;
}

RoundedConf round_conf(double lo, double val, double hi) {
  RoundedErr elo = round_err(val - lo);
  RoundedErr ehi = round_err(hi - val);
  double min_digit = std::min(elo.digit, ehi.digit);
  if (elo.digit != ehi.digit) {
    elo = round_err_digits(val - lo, min_digit);
    ehi = round_err_digits(hi - val, min_digit);
  }
  double v = std::round(val * std::pow(10.0, -min_digit)) *
             std::pow(10.0, min_digit);
  v = atof_fmt(v, min_digit);
  return RoundedConf{elo.value, v, ehi.value, min_digit};
}

namespace {
std::string tex_impl(double val, double mn, double mx, bool has_limits,
                     double min_limit, double max_limit, int sci) {
  bool min_is_at_limit = false, max_is_at_limit = false;
  if (has_limits) {
    if (mn == min_limit) min_is_at_limit = true;
    if (mx == max_limit) max_is_at_limit = true;
  }
  RoundedConf s = round_conf(mn, val, mx);
  val = s.value;
  mn = val - s.err_lo;
  mx = val + s.err_hi;
  double dig = s.digit;
  if (has_limits) {
    mn = std::max(mn, min_limit);
    mx = std::min(mx, max_limit);
  }
  double max_abs = std::max({std::fabs(mn), std::fabs(mx), std::fabs(val)});
  if (max_abs == 0.0) return "$\\equiv0$";

  int scismall = -sci, scilarge = sci;
  std::string power_str;
  double power = 0.0;
  double val4power = (val != 0.0 ? std::fabs(val) : max_abs);
  if (!(std::pow(10.0, scismall) <= val4power &&
        val4power < std::pow(10.0, scilarge))) {
    power = double(int(std::log10(val4power)));
    if (power < 0) power -= 1;
    double pow10 = std::pow(10.0, -power);
    val *= pow10;
    mn *= pow10;
    mx *= pow10;
    char buf[64];
    std::snprintf(buf, sizeof buf, "\\times10^{%d}", int(power));
    power_str = buf;
  }
  int n = int(std::lround(power - dig));
  if (n < 0) n = 0;

  char fmt[16];
  std::snprintf(fmt, sizeof fmt, "%%.%df", n);
  char b1[64], b2[64], b3[64];
  std::snprintf(b1, sizeof b1, fmt, val - mn);
  std::snprintf(b2, sizeof b2, fmt, mx - val);
  std::string m_err = b1, p_err = b2;
  // phantom padding for alignment
  int diff = int(p_err.size()) - int(m_err.size());
  if (diff != 0) {
    std::string zeros(std::abs(diff), '0');
    if (diff > 0) m_err = "\\phantom{" + zeros + "}" + m_err;
    if (diff < 0) p_err = "\\phantom{" + zeros + "}" + p_err;
  }
  if (min_is_at_limit) m_err = "{\\color{red}" + m_err + "}";
  if (max_is_at_limit) p_err = "{\\color{red}" + p_err + "}";

  std::string err = (p_err == m_err) ? (" \\pm " + p_err)
                                     : ("^{+" + p_err + "}_{-" + m_err + "}");
  std::snprintf(b3, sizeof b3, fmt, val);
  std::string vstr = b3;
  char b4[64];
  std::snprintf(b4, sizeof b4, fmt, 0.0);
  std::string null = b4;
  if (mn >= 0 && vstr == null && m_err == null) {
    char b5[96];
    std::string f = std::string("$\\le") + fmt + "%s$";
    std::snprintf(b5, sizeof b5, f.c_str(), mx, power_str.c_str());
    return b5;
  }
  bool paren = !power_str.empty();
  return "$" + std::string(paren ? "\\left(" : "") + vstr + err +
         (paren ? "\\right)" : "") + power_str + "$";
}
}  // namespace

std::string tex_value_pm_error(double value, double mn, double mx, int sci) {
  return tex_impl(value, mn, mx, false, 0, 0, sci);
}

std::string tex_value_pm_error(double value, double mn, double mx,
                               double min_limit, double max_limit, int sci) {
  return tex_impl(value, mn, mx, true, min_limit, max_limit, sci);
}

}  // namespace sed
