// gpu_firth_stub.cpp — the no-CUDA implementation of gpu_firth.hpp. Linked when
// the Makefile is invoked without USE_CUDA=1; firthCreate() says no, so the
// gpuFirth switch can never take effect and Firth stays on the CPU scalar path.
#include "gpu_firth.hpp"

namespace saige {
namespace gpu2 {

bool          firthSupports(Prec t_p)                          { return t_p == Prec::FP64; }
Firth*        firthCreate(const FirthCreateArgs&)              { return nullptr; }
void          firthDestroy(Firth*)                             {}
FirthPairIn*  firthIn(Firth*)                                  { return nullptr; }
FirthPairOut* firthOut(Firth*)                                 { return nullptr; }
double*       firthPairLut(Firth*)                             { return nullptr; }
bool          firthRun(Firth*, const Reducer*, int, int)            { return false; }
void          firthTimings(const Firth*, double*, long long*)  {}
std::size_t   firthDeviceBytes(const Firth*)                   { return 0; }
const char*   firthLastError()                                 { return "no CUDA build"; }

}  // namespace gpu2
}  // namespace saige
