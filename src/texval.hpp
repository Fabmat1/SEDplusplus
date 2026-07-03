// Ports of round_err, round_conf and TeX_value_pm_error (DIN 1333 rounding,
// TeX formatting of value +- error), used for the 'tex' results column.
#pragma once

#include <string>

namespace sed {

struct RoundedErr {
  double value;
  double digit;
};

RoundedErr round_err(double x);
RoundedErr round_err_digits(double x, double digits);

struct RoundedConf {
  double err_lo, value, err_hi, digit;
};

RoundedConf round_conf(double lo, double val, double hi);

// TeX_value_pm_error(value, min, max [, min_limit, max_limit]; sci=6)
std::string tex_value_pm_error(double value, double mn, double mx, int sci);
std::string tex_value_pm_error(double value, double mn, double mx,
                               double min_limit, double max_limit, int sci);

}  // namespace sed
