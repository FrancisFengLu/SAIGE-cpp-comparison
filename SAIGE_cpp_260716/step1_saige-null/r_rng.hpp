// Standalone port of the parts of R's random number generation that step 1
// uses, so that saige-null no longer needs an embedded R runtime to draw the
// Hutchinson trace probes. Draws are bit-identical to R (>= 3.6) with the
// default RNGkind("Mersenne-Twister", "Inversion", "Rejection").
//
// Ported from R 4.5.2 sources (GPL-2 or later, same as this repository):
//   src/main/RNG.c     RNG_Init (initial LCG scrambling), FixupSeeds,
//                      MT_sgenrand, MT_genrand, fixup, unif_rand
//   src/nmath/rbinom.c inversion branch (n*p < 30)
//
// The state is one process-wide stream, like R's .Random.seed: set_seed()
// reseeds it and every later draw continues it. R's GetRNGstate()/PutRNGstate()
// only copy that state to and from .Random.seed, so they have no counterpart
// here. Not thread-safe (neither is R's generator); call from one thread.
#pragma once

namespace saige {
namespace rrng {

// set.seed(seed) with the default RNG kind (Mersenne-Twister).
void set_seed(int seed);

// unif_rand(): uniform on (0,1), never exactly 0 or 1.
double unif_rand();

// rbinom(n, p) for n*p < 30 (R's inversion branch). n must be a
// non-negative integer value < INT_MAX and 0 <= p <= 1; the BTPE branch
// (n*min(p,1-p) >= 30) is not ported and throws std::domain_error.
double rbinom(double n, double p);

}  // namespace rrng
}  // namespace saige
