// Port of R's Mersenne-Twister unif_rand and rbinom; see r_rng.hpp.
// Source: R 4.5.2, src/main/RNG.c and src/nmath/rbinom.c
// (Copyright R Core Team and others, GPL-2 or later).
#include "r_rng.hpp"

#include <cmath>
#include <climits>
#include <cstdint>
#include <stdexcept>

namespace saige {
namespace rrng {

namespace {

typedef std::uint32_t Int32;  // R: typedef unsigned int Int32

// RNG.c: #define i2_32m1 2.328306437080797e-10  /* = 1/(2^32 - 1) */
const double i2_32m1 = 2.328306437080797e-10;

// MT period parameters (RNG.c)
const int N = 624;
const int M = 397;
const Int32 MATRIX_A   = 0x9908b0dfU;
const Int32 UPPER_MASK = 0x80000000U;
const Int32 LOWER_MASK = 0x7fffffffU;
const Int32 TEMPERING_MASK_B = 0x9d2c5680U;
const Int32 TEMPERING_MASK_C = 0xefc60000U;

// RNG.c: static Int32 dummy[625]; dummy[0] is mti, mt = dummy + 1.
// R starts with mti = N+1 (not seeded); an unseeded draw then falls back to
// MT_sgenrand(4357). R itself would have seeded from the clock first
// (Randomize), so draws without set_seed() are well-defined here but not R's.
Int32 dummy[N + 1] = { (Int32)(N + 1) };
Int32* const mt = dummy + 1;

// RNG.c: fixup() -- ensure 0 and 1 are never returned
double fixup(double x) {
  if (x <= 0.0) return 0.5 * i2_32m1;
  if ((1.0 - x) <= 0.0) return 1.0 - 0.5 * i2_32m1;
  return x;
}

// RNG.c: MT_sgenrand()
void MT_sgenrand(Int32 seed) {
  for (int i = 0; i < N; i++) {
    mt[i] = seed & 0xffff0000U;
    seed = 69069 * seed + 1;
    mt[i] |= (seed & 0xffff0000U) >> 16;
    seed = 69069 * seed + 1;
  }
  dummy[0] = N;
}

// RNG.c: MT_genrand()
double MT_genrand() {
  static const Int32 mag01[2] = { 0x0U, MATRIX_A };
  Int32 y;
  int mti = (int)dummy[0];

  if (mti >= N) {
    int kk;
    if (mti == N + 1) MT_sgenrand(4357);
    for (kk = 0; kk < N - M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
    mti = 0;
  }

  y = mt[mti++];
  y ^= (y >> 11);
  y ^= (y << 7) & TEMPERING_MASK_B;
  y ^= (y << 15) & TEMPERING_MASK_C;
  y ^= (y >> 18);
  dummy[0] = (Int32)mti;

  return ((double)y * 2.3283064365386963e-10);  // [0,1)
}

}  // namespace

// RNG.c: do_setseed -> RNG_Init(MERSENNE_TWISTER, (Int32) seed) -> FixupSeeds
void set_seed(int seed_in) {
  Int32 seed = (Int32)seed_in;
  for (int j = 0; j < 50; j++)          // initial scrambling
    seed = (69069 * seed + 1);
  for (int j = 0; j < N + 1; j++) {     // n_seed = 625; i_seed[0] is mti
    seed = (69069 * seed + 1);
    dummy[j] = seed;
  }
  // FixupSeeds(MERSENNE_TWISTER, initial = 1)
  dummy[0] = N;
  bool notallzero = false;
  for (int j = 1; j <= N; j++)
    if (dummy[j] != 0) { notallzero = true; break; }
  // R would call Randomize() (a clock seed) here; no 32-bit seed reaches this.
  if (!notallzero) throw std::logic_error("rrng::set_seed: all-zero MT state");
}

// RNG.c: unif_rand(), case MERSENNE_TWISTER
double unif_rand() { return fixup(MT_genrand()); }

// rbinom.c, restricted to the inversion branch. R caches qn in statics keyed on
// (p, n); recomputing it every call gives the same value.
double rbinom(double nin, double pp) {
  if (!std::isfinite(nin)) return NAN;
  double r = std::nearbyint(nin);      // R_forceint
  if (r != nin) return NAN;
  if (!std::isfinite(pp) || r < 0 || pp < 0. || pp > 1.) return NAN;

  if (r == 0 || pp == 0.) return 0;
  if (pp == 1.) return r;

  if (r >= INT_MAX)
    throw std::domain_error("rrng::rbinom: n >= INT_MAX (qbinom branch) not ported");
  const int n = (int)r;

  const double p = std::fmin(pp, 1. - pp);
  const double q = 1. - p;
  const double np = n * p;
  r = p / q;
  const double g = r * (n + 1);

  if (!(np < 30.0))
    throw std::domain_error("rbinom: n*p >= 30 (BTPE branch) not ported");

  // R_pow_di(q, n)
  double qn = 1.;
  {
    double x = q;
    int k = n;
    for (;;) {
      if (k & 01) qn *= x;
      if (k >>= 1) x *= x; else break;
    }
  }

  int ix;
  for (;;) {
    ix = 0;
    double f = qn;
    double u = unif_rand();
    for (;;) {
      if (u < f) goto finis;
      if (ix > 110) break;
      u -= f;
      ix++;
      f *= (g / ix - r);
    }
  }
finis:
  if (pp > 0.5) ix = n - ix;           // psave > 0.5
  return (double)ix;
}

}  // namespace rrng
}  // namespace saige
