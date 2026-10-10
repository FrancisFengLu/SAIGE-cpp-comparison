// s2_gpu_sparse.hpp -- host half of the step-2 GPU path's sparse-GRM variance
// (config key gpuSparse; gpu/gpu_sparse.hpp has the math).
//
// Built once per run, before the reducer is created: for every trait whose
// model carries a sparse GRM, Sigma (as SAIGEClass holds it, tau1*K with
// 1/mu2 -- binary -- or tau0 -- quantitative -- on the diagonal) is partitioned into its connected components and
// every block is inverted densely with the same class the CPU's
// blockSparseSigma switch uses (block_sigma.hpp, fp64), so B = Sigma^-1 is
// held explicitly.
//
// The CPU's sparse score test (SAIGEClass::scoreTest) forms
//     g~ = g - Y (XV g),  Y = XXVX_inv (N x p), XV (p x N)
//     S  = g~' res / tau0,   var2 = g~' B g~
// -- with XV / XXVX_inv, NOT the XVX_inv_XV / S_a the dense fast test uses;
// step 1 stores those separately and they agree only to fp32-level, and S_a is
// the fit-time X'res while scoreTest dots with res itself. So with z = XV g:
//     S    = (g'res - z'(Y'res)) / tau0
//     var2 = g'Bg - 2 z'(BY)'g + z'(Y'BY)z
// The run needs, per trait:
//   XV'    N x p    one block of columns in the reducer's first GEMM -> z
//   BY     N x p    another block                                    -> (BY)'g
//   Y'res  p, Y'BY p x p    host contractions
//   diag B N        one more column of the second GEMM, (G % G)' diag(B)
//   2 B_ij          for i < j inside a block: the cross terms, on the device
//                   (gpu_sparse.cu); the pair list is the union over traits
//                   (a trait with tau1 = 0 has no off-diagonal entries left in
//                   its Sigma and contributes zero weights)
// so that every marker's exact sparse variance comes without a per-marker solve.
#ifndef S2_GPU_SPARSE_HPP
#define S2_GPU_SPARSE_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "saige_mt.hpp"
#include "saige_test.hpp"

namespace s2gs {

struct Plan {
    int N = 0, P = 0, nSp = 0, sumPsp = 0;
    // Sparse traits, binary or quantitative, in internal trait order: spIdx[t]
    // is the trait's sparse index s (-1 without a sparse GRM), spTrait[s] the
    // trait, spOff[s] the first column of its p-block inside XVt / BY.
    std::vector<int> spIdx, spTrait, spOff;
    arma::mat XVt;                 // N x sumPsp, trait block at spOff: XV'
    arma::mat BY;                  // N x sumPsp, trait block at spOff: B XXVX_inv
    arma::mat Bdiag;               // N x nSp
    std::vector<arma::vec> Yres;   // [P] p: XXVX_inv' res; empty when not sparse
    std::vector<arma::mat> YBY;    // [P] p x p: XXVX_inv' B XXVX_inv; empty when not sparse
    std::vector<int> pi, pj;       // within-block pairs, i < j, union over traits
    std::vector<double> w;         // nSp x nPairs, sparse-trait-major, 2 B_ij (0 where absent)
    long long nPairs = 0;
    int nBlocksMax = 0, maxBlock = 0, nOn = 0;
    // Different sample lists (gpuOwnSampleSets) only. Every per-sample array
    // above is union-length with exact zeros outside the trait's samples, and
    // pi / pj are union indices. A trait with its own list is scored on its
    // own genotype vector g_t = a g + b 1_t + d m (MTBlockAdj), so its tail
    // also needs the sums of its columns over its samples:
    std::vector<arma::vec> sumXV;  // [P] p: XV 1_t
    std::vector<arma::vec> sumBY;  // [P] p: (B XXVX_inv)' 1_t
    std::vector<double>    trB;    // [P] sum_i B_ii
    double secs = 0.0;
};

// Builds the plan for every binary or quantitative trait with flagSparseGRM
// (the quantitative Sigma carries tau0 on its diagonal, the binary one 1/mu2;
// both come from the loader as SAIGEClass holds them). Returns "" on
// success, otherwise the reason the device sparse variance is not available
// (the caller then refuses or keeps the CPU path). refreshBudget_s is the
// block inverse's one-refresh cost gate (blockSparseSigmaRefreshBudget_s).
// maxPairs caps the within-block pair count (the device kernel's cost is
// pairs x markers). With different sample lists (ctx.sampleSetsDiffer) each
// trait's Sigma is over its own samples; its blocks are mapped onto union
// indices through ctx.samp[t].pos.
std::string build(const SAIGE::MTContext& ctx,
                  const std::vector<SAIGE::SAIGEClass*>& objs,
                  double refreshBudget_s, long long maxPairs, Plan& out);

}  // namespace s2gs

#endif  // S2_GPU_SPARSE_HPP
