// spa_gpu_stub.cpp -- the no-CUDA build of spa_gpu.hpp. create() refuses, so
// every caller keeps SPA on the CPU; nothing else is reachable.
#include "spa_gpu.hpp"

namespace saige {
namespace spa_gpu {

bool supports(saige::gpu2::Prec t_p) { return t_p == saige::gpu2::Prec::FP64; }
Spa* create(const CreateArgs&) { return nullptr; }
void destroy(Spa*) {}
PairIn*  in(Spa*)  { return nullptr; }
PairOut* out(Spa*) { return nullptr; }
double*  pairLut(Spa*) { return nullptr; }
bool run(Spa*, const Geno&, int) { return false; }
bool uploadPacked(Spa*, const unsigned char*, std::size_t, const double*, int, Geno*) { return false; }
bool uploadDense(Spa*, const double*, int, Geno*) { return false; }
void timings(const Spa*, double* k, double* h, double* d, long long* n)
{
    if (k) *k = 0; if (h) *h = 0; if (d) *d = 0; if (n) *n = 0;
}
std::size_t deviceBytes(const Spa*) { return 0; }
const char* lastError() { return "built without CUDA (spa_gpu_stub.cpp)"; }
bool debugErfc(int, const double*, int, double*, double*, double*) { return false; }

}  // namespace spa_gpu
}  // namespace saige
