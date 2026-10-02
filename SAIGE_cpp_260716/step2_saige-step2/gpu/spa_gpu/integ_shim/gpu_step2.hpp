// integ_shim/gpu_step2.hpp -- stand-in for the step-2 reducer header, used
// ONLY to compile the integrator's gpu/gpu_spa.cu against spa_gpu_test
// without the reducer. Its kernel reads the genotype rows through these three
// accessors; here a Reducer is just the device pointers spa_gpu_test
// uploaded. Selected by putting this directory first on the include path
// (Makefile target spa_gpu_test_integ).
#pragma once
#include <cstddef>

namespace saige {
namespace gpu2 {

struct Reducer {
    const void*  pk  = nullptr;   // device: nSlots x bpv packed rows
    const void*  lut = nullptr;   // device: nSlots x 4 doubles
    std::size_t  bpv = 0;
};

inline const void* devicePacked(const Reducer* r) { return r ? r->pk : nullptr; }
inline const void* deviceLut(const Reducer* r)    { return r ? r->lut : nullptr; }
inline std::size_t bytesPerSlot(const Reducer* r) { return r ? r->bpv : 0; }

}  // namespace gpu2
}  // namespace saige
