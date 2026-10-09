// gpu_er_fp32.cuh — internal interface between gpu_er.cu (context, uploads,
// dispatch) and gpu_er_fp32.cu (the fp32 kernels of gpuPrecisionER: fp32).
// Not part of the public gpu_er.hpp contract.
#pragma once

#include <cstddef>
#include <cstdint>

#include <cuda_runtime.h>

#include "gpu_er.hpp"

namespace saige {
namespace gpu2 {
namespace erfp32 {

// Bytes of one large-k (k > 12) hand-off record between the two kernels.
std::size_t prepBytes();

// Launches the fp32 pair kernel on all pairs and the fp32 block kernel on the
// pairs with k > 12 (bigPair / bigSlot as in the fp64 path), on stream st.
// mu / muOff / nOf / ncaseOf are the fp64 path's resident arrays; muSum holds
// per trait the sum of its mu as a float pair (hi, lo), filled at erCreate.
// carG is the double dosage array of erRun. ErPairOut::pval is written as double.
// Returns the first launch error (cudaSuccess otherwise).
cudaError_t launch(cudaStream_t st,
                   const double* mu, const long long* muOff, const int* nOf, const int* ncaseOf,
                   const float* muSum,
                   const ErPairIn* in, int nPairs,
                   const uint32_t* carIdx, const double* carG, const unsigned char* carCase,
                   ErPairOut* out, const int* bigSlot, const int* bigPair, int nBig, void* prep);

}  // namespace erfp32
}  // namespace gpu2
}  // namespace saige
