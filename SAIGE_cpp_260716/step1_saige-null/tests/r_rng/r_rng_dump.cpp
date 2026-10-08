// Dump draws from the r_rng port as raw doubles, to compare byte-for-byte with
// R (see run_r_rng_check.sh). Usage: r_rng_dump <seed> <n> <unif|rbinom> > out.bin
#include "../../r_rng.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
int main(int argc, char** argv) {
  if (argc != 4) { std::fprintf(stderr, "usage: %s seed n unif|rbinom\n", argv[0]); return 2; }
  const int seed = std::atoi(argv[1]);
  const long n = std::atol(argv[2]);
  const bool unif = std::strcmp(argv[3], "unif") == 0;
  std::vector<double> v(n);
  saige::rrng::set_seed(seed);
  for (long i = 0; i < n; ++i) v[i] = unif ? saige::rrng::unif_rand() : saige::rrng::rbinom(1.0, 0.5);
  std::fwrite(v.data(), sizeof(double), n, stdout);
  return 0;
}
