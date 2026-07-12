// Deterministic fast 10^x for the extinction factor (config fast_ext=1).
// Accuracy ~1e-15 relative (vs ~0.5 ulp for libm pow(10,x)): utterly
// negligible against photometric uncertainties, but NOT bit-identical to
// pow -- which is why it is opt-in. Branch-light and self-contained, so the
// result does not depend on the libm version and the loop auto-vectorizes.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>

namespace sed {

// Vectorized batch kernel (fastext.cpp, compiled with -mavx2 -mfma):
// out[k] = 10^((-0.4*E) * (k0[k] + sl[k]*(R-3.02) + R))
void ext_factor_fast(const double* k0, const double* sl, double E, double R,
                     double* out, std::size_t n);

inline double fast_exp10(double x) {
  // 10^x = 2^n * e^s  with  t = x*log2(10), n = round(t),
  // s = (t - n)*ln(2)  in [-0.347, 0.347]
  constexpr double LOG2_10 = 3.3219280948873623478703194294894;
  constexpr double LN2 = 0.69314718055994530941723212145818;
  const double t = x * LOG2_10;
  const double n = std::nearbyint(t);
  const double s = (t - n) * LN2;
  // e^s Taylor series to s^13 (remainder < 1.6e-16 for |s| <= 0.347)
  double p = 1.0 / 6227020800.0;  // 1/13!
  p = p * s + 1.0 / 479001600.0;
  p = p * s + 1.0 / 39916800.0;
  p = p * s + 1.0 / 3628800.0;
  p = p * s + 1.0 / 362880.0;
  p = p * s + 1.0 / 40320.0;
  p = p * s + 1.0 / 5040.0;
  p = p * s + 1.0 / 720.0;
  p = p * s + 1.0 / 120.0;
  p = p * s + 1.0 / 24.0;
  p = p * s + 1.0 / 6.0;
  p = p * s + 0.5;
  p = p * s + 1.0;
  p = p * s + 1.0;
  // scale by 2^n via exponent bits (n is clamped to the finite range; the
  // extinction exponent never comes near it in practice)
  double ni = n;
  if (ni > 1023.0) ni = 1023.0;
  if (ni < -1022.0) ni = -1022.0;  // no denormals: 10^-307 is flush-to-0 land
  const std::uint64_t bits = (std::uint64_t)(std::int64_t)(ni + 1023.0) << 52;
  double scale;
  std::memcpy(&scale, &bits, 8);
  return p * scale;
}

}  // namespace sed
