// gpu_precision.hpp — per-stage arithmetic precision of the step-2 GPU path.
//
// Four stages, each chosen independently by its own config key (main.cpp);
// any combination is allowed:
//
//   stage   config key          modes               module(s)
//   scan    gpuPrecisionScan    fp64 | fp32 | int8  gpu_step2.cu (decode + GEMMs),
//           (alias gpuPrecision)                    gpu_sparse.cu (sparse-GRM cross terms)
//   SPA     gpuPrecisionSPA     fp64 | fp32         spa_gpu/spa_gpu.cu (gpuSpaImpl: lib),
//                                                   gpu_spa.cu (gpuSpaImpl: own)
//   ER      gpuPrecisionER      fp64 | fp32         gpu_er.cu
//   Firth   gpuPrecisionFirth   fp64 | fp32         gpu_firth.cu
//
// int8 (scan only) is the Ozaki-style split: genotypes are exact in int8, the
// trait-side matrix is split into gpuInt8Slices (1..8, default 7) int8 slices
// with per-column scaling, int8 x int8 -> int32 GEMMs via cublasGemmEx, the
// partial products recombined in fp64. The reducer hands back double results
// in that mode.
//
// fp64 everywhere is the default and is today's code, bit for bit.
//
// Every module exports `bool <module>Supports(Prec)`: true when its create()
// accepts that mode. main.cpp asks before it creates anything and stops the
// run with "useGPU: <stage> precision <mode> is not implemented yet" when the
// answer is no -- a requested precision is never silently replaced by another.
// A kernel agent that implements a mode flips its module's Supports() and
// plugs the variant in at the module's dispatch point (marked
// "TODO(precision:<stage>)"). Nothing else needs to change in main.cpp.
//
// Header-only and CUDA-free, so main.cpp, the stubs and the .cu files can all
// include it.
#pragma once

namespace saige {
namespace gpu2 {

enum class Prec : int { FP64 = 0, FP32 = 1, INT8 = 2 };

inline const char* precName(Prec p)
{
    switch (p) {
        case Prec::FP64: return "fp64";
        case Prec::FP32: return "fp32";
        case Prec::INT8: return "int8";
    }
    return "?";
}

constexpr int kInt8SlicesDefault = 7;
constexpr int kInt8SlicesMin = 1;
constexpr int kInt8SlicesMax = 8;

}  // namespace gpu2
}  // namespace saige
