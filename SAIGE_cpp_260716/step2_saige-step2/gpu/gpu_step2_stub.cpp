// gpu_step2_stub.cpp — the no-CUDA implementation of the gpu_step2.hpp facade.
//
// Built instead of gpu_step2.cu when the Makefile is invoked without USE_CUDA=1
// (the default). Everything reports "unavailable", so main()'s GPU gate is
// simply never satisfied and the CPU path runs exactly as it did before the
// GPU work existed. This is what keeps the default build byte-identical
// without a single #ifdef in main.cpp.

#include "gpu_step2.hpp"

namespace saige {
namespace gpu2 {

bool available(int, std::string* t_why)
{
    if (t_why) *t_why = "built without CUDA (rebuild with: make USE_CUDA=1)";
    return false;
}

std::string describe(int) { return std::string(); }

Reducer* create(int, int, int, const double*, int, bool) { return nullptr; }
void     destroy(Reducer*) {}

unsigned char* packed(Reducer*)            { return nullptr; }
double*        lut(Reducer*)               { return nullptr; }
std::size_t    bytesPerSlot(const Reducer*) { return 0; }
bool           reduce(Reducer*, int)       { return false; }
bool           isFp64(const Reducer*)       { return false; }
const float*   outCf(const Reducer*)       { return nullptr; }
const double*  outCd(const Reducer*)       { return nullptr; }
std::size_t    ldC(const Reducer*)         { return 0; }
std::size_t    deviceBytes(const Reducer*) { return 0; }
void           timings(const Reducer*, double*, double*, double*, double*) {}

}  // namespace gpu2
}  // namespace saige
