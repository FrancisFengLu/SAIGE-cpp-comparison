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

int maskWords(int t_N) { return (t_N + 31) / 32; }

bool     scanSupports(Prec t_p) { return t_p == Prec::FP64; }
Reducer* create(const CreateArgs&) { return nullptr; }
void     destroy(Reducer*) {}

unsigned char*  packed(Reducer*, int)         { return nullptr; }
double*         lut(Reducer*, int)            { return nullptr; }
std::size_t     bytesPerSlot(const Reducer*)  { return 0; }
int             stagingSets(const Reducer*)   { return 0; }
bool            reduce(Reducer*, int, int, int) { return false; }
bool            isFp64(const Reducer*)        { return false; }
const float*    outCf(const Reducer*, int)    { return nullptr; }
const double*   outCd(const Reducer*, int)    { return nullptr; }
const float*    outC2f(const Reducer*, int)   { return nullptr; }
const double*   outC2d(const Reducer*, int)   { return nullptr; }
std::size_t     ldC(const Reducer*)           { return 0; }
const uint32_t* outCounts(const Reducer*, int) { return nullptr; }
std::size_t     deviceBytes(const Reducer*)   { return 0; }
void            timings(const Reducer*, double*, double*, double*, double*, double*) {}
const void*     devicePacked(const Reducer*, int) { return nullptr; }
const void*     deviceLut(const Reducer*, int) { return nullptr; }
void*           deviceStream(const Reducer*)  { return nullptr; }
int             deviceSets(const Reducer*)    { return 0; }
bool            bindDevice(int)               { return false; }
std::string     setBlockingSync(int)          { return "built without CUDA"; }

bool                 statsSetup(Reducer*, const StatsArgs&)  { return false; }
void                 statsDisable(Reducer*)                  {}
std::string          statsLastError()                        { return "built without CUDA"; }
double*              statsVr(Reducer*, int)                  { return nullptr; }
double*              statsAf(Reducer*, int)                  { return nullptr; }
const double*        statsS(const Reducer*, int)             { return nullptr; }
const double*        statsVar2(const Reducer*, int)          { return nullptr; }
const double*        statsP(const Reducer*, int)             { return nullptr; }
const unsigned char* statsFlags(const Reducer*, int)         { return nullptr; }
double               statsSeconds(const Reducer*)            { return 0.0; }

}  // namespace gpu2
}  // namespace saige
