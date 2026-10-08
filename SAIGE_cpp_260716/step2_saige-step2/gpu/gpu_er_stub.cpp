// gpu_er_stub.cpp — the no-CUDA implementation of gpu_er.hpp. Linked when the
// Makefile is invoked without USE_CUDA=1; erCreate() says no, so the gpuER
// switch can never take effect and ER stays on the CPU scalar path.
#include "gpu_er.hpp"

namespace saige {
namespace gpu2 {

bool        erSupports(Prec t_p)                             { return t_p == Prec::FP64; }
Er*         erCreate(const ErCreateArgs&)                    { return nullptr; }
void        erDestroy(Er*)                                   {}
bool        erRun(Er*, const ErPairIn*, int, const uint32_t*, const double*,
                  const unsigned char*, long long, ErPairOut*) { return false; }
void        erTimings(const Er*, double* t, long long* p)    { if (t) *t = 0; if (p) *p = 0; }
const char* erLastError()                                    { return "no CUDA build"; }
bool        erMathCheck(int, const double*, long long, double*, double*) { return false; }

}  // namespace gpu2
}  // namespace saige
