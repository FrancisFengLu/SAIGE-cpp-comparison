// gpu_spa_stub.cpp — the no-CUDA implementation of gpu_spa.hpp. Linked when
// the Makefile is invoked without USE_CUDA=1; spaCreate() says no, so the
// gpuSpa switch can never take effect and SPA stays on the CPU scalar path.
#include "gpu_spa.hpp"

namespace saige {
namespace gpu2 {

bool        spaSupports(Prec t_p)                        { return t_p == Prec::FP64; }
Spa*        spaCreate(const SpaCreateArgs&)              { return nullptr; }
void        spaDestroy(Spa*)                             {}
SpaPairIn*  spaIn(Spa*)                                  { return nullptr; }
SpaPairOut* spaOut(Spa*)                                 { return nullptr; }
bool        spaRun(Spa*, const Reducer*, int, int)            { return false; }
void        spaTimings(const Spa*, double*, long long*)  {}
std::size_t spaDeviceBytes(const Spa*)                   { return 0; }
const char* spaLastError()                               { return "built without CUDA"; }

}  // namespace gpu2
}  // namespace saige
