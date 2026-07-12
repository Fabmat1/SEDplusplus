// Vectorized extinction-factor kernel for config fast_ext=1 (see
// fastmath.hpp for the accuracy statement: ~3e-16 relative vs libm pow,
// deterministic, NOT byte-exact with the S-Lang reference -- opt-in only).
// This file is compiled with -O3 -mavx2 -mfma (see CMakeLists.txt); nothing
// on the byte-parity default path lives here.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace sed {

// out[k] = 10^((-0.4*E) * (k0[k] + sl[k]*(R-3.02) + R))
void ext_factor_fast(const double* k0, const double* sl, double E, double R,
                     double* out, std::size_t n) {
  constexpr double LOG2_10 = 3.3219280948873623478703194294894;
  constexpr double LN2 = 0.69314718055994530941723212145818;
  constexpr double SHIFT = 6755399441055744.0;  // 1.5*2^52: |t|<2^51 rounder
  const double a = -0.4 * E;
  const double dR = R - 3.02;
  for (std::size_t k = 0; k < n; ++k) {
    const double x = a * (k0[k] + sl[k] * dR + R);
    const double t =
        std::min(1023.0, std::max(-1022.0, x * LOG2_10));  // branchless clamp
    const double nf = t + SHIFT;  // round-to-nearest-even into the mantissa
    std::uint64_t nbits;
    std::memcpy(&nbits, &nf, 8);
    const double nd = nf - SHIFT;
    const double s = (t - nd) * LN2;
    // e^s Taylor to s^13 (|s| <= 0.347 -> remainder < 1.6e-16)
    double p = 1.0 / 6227020800.0;
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
    // 2^nd from the exponent bits: mantissa low 52 bits of nf hold 2^51 + n
    const std::int64_t ni =
        (std::int64_t)(nbits & 0xFFFFFFFFFFFFFULL) - (1LL << 51);
    const std::uint64_t ebits = (std::uint64_t)(ni + 1023) << 52;
    double scale;
    std::memcpy(&scale, &ebits, 8);
    out[k] = p * scale;
  }
}

}  // namespace sed
