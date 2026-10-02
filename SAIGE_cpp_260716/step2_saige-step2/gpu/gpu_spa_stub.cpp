// gpu_spa_stub.cpp — the no-CUDA implementation of gpu_spa.hpp. Linked when
// the Makefile is invoked without USE_CUDA=1; spaCreate() says no, so the
// gpuSpa switch can never take effect and SPA stays on the CPU scalar path.
#include "gpu_spa.hpp"

namespace saige {
namespace gpu2 {

Spa*        spaCreate(const SpaCreateArgs&)              { return nullptr; }
void        spaDestroy(Spa*)                             {}
SpaPairIn*  spaIn(Spa*)                                  { return nullptr; }
SpaPairOut* spaOut(Spa*)                                 { return nullptr; }
bool        spaRun(Spa*, const Reducer*, int)            { return false; }
void        spaTimings(const Spa*, double*, long long*)  {}
std::size_t spaDeviceBytes(const Spa*)                   { return 0; }

}  // namespace gpu2
}  // namespace saige
