// gpu_sparse_stub.cpp — the no-CUDA implementation of gpu_sparse.hpp. Linked
// when the Makefile is invoked without USE_CUDA=1; spqCreate() says no, so the
// gpuSparse switch can never take effect.
#include "gpu_sparse.hpp"

namespace saige {
namespace gpu2 {

bool          spqSupports(Prec t_p)                                { return t_p == Prec::FP64; }
SpQuad*       spqCreate(const SpQuadCreateArgs&)                   { return nullptr; }
void          spqDestroy(SpQuad*)                                  {}
bool          spqRun(SpQuad*, const Reducer*, int, int)            { return false; }
bool          spqRunOwn(SpQuad*, const Reducer*, int, const double*, int) { return false; }
bool          spqRunOwnDev(SpQuad*, const Reducer*, int, const void*, int) { return false; }
const double* spqOut(const SpQuad*, int)                           { return nullptr; }
const void*   spqDeviceOut(const SpQuad*)                          { return nullptr; }
void          spqTimings(const SpQuad*, double*, long long*)       {}
std::size_t   spqDeviceBytes(const SpQuad*)                        { return 0; }
const char*   spqLastError()                                       { return "no CUDA build"; }

}  // namespace gpu2
}  // namespace saige
