// spa_gpu.hpp -- SAIGE's binary-trait saddlepoint approximation (SPA) on the
// device, for a table of flagged (marker, trait) pairs. Self-contained: no
// CUDA type crosses this header, and the genotype source is handed over as
// device pointers (the step-2 reducer's resident packed rows, or a dense
// buffer), so the library neither owns nor needs the reducer.
//
// Contract: the per-pair inputs and outputs are the ones gpu/gpu_spa.hpp
// (the integrator's draft) uses -- slot, trait, fast, logp, Tstat, var1,
// var2, pno in; pval, conv, s1, s2, niter, roots, m1 out -- plus a status
// word that says which variant ran and which fallback, if any, the pair
// took. A caller written against gpu_spa.hpp swaps in with the adapter in
// SPA_GPU_LIB.md (section "Swapping in").
//
// ---------------------------------------------------------------------------
// What one pair computes (SAIGEClass::getMarkerPval, saige_test.cpp, binary)
// ---------------------------------------------------------------------------
//   g~      = g - XXVX_inv (XV g)             getadjGFast: XV g summed over
//                                             carriers (g != 0) only
//   m1      = mu . g~
//   q       = Tstat / sqrt(var1 / var2) + m1,  qinv mirrored about m1
//   variant fast  iff  |{g == 0}| / N >= 0.5 (and no sparse GRM; host's call):
//             sums over carriers, non-carriers folded into
//             NAmu = m1 - sum_c mu g~,  NAsigma = var2 - sum_c mu (1-mu) g~^2
//   roots   two Newton solves of K1(t) = q resp. qinv from t = 0, with the
//           prevJump halving on a sign change, |dt| < tol, maxiter; root =
//           +inf, niter 0 when q lies outside (sum g~<0, sum g~>0)
//   tails   two saddlepoint probabilities; a tail that is not a saddle is
//           replaced by pno/2 (pno - log 2 in the log domain) and withdraws
//           convergence (spa.cpp SPA / SPA_fast); |p1| + |p2|, or add_logp
//   pno     when either root fails: pval = pno, conv = false
//
// spa_binary.cpp and spa.cpp are the reference. The kernel carries their
// control flow over statement by statement, one CUDA block per pair, the
// block's threads sharing every length-N (or length-carriers) sum, so each
// pair runs exactly as many Newton steps as SAIGE's loop. The two variants
// keep their own quirks: only the full-N root-find has the K2 ~ 0 guard and
// the !isfinite(tnew) test (the fast one tests isnan), and they detect a sign
// change differently (sign() != sign() vs product < 0).
//
// What differs from the CPU: the order of the additions inside each sum
// (tree instead of Armadillo's two interleaved accumulators / OpenBLAS), and
// exp / log from CUDA's libm instead of glibc's. Everything else is matched:
//   * per-term arithmetic is associated exactly as Armadillo's expression
//     templates evaluate it, and the build uses --fmad=false (the CPU is
//     -std=c++17, so it does not contract either);
//   * the normal tail: the CPU's boost::math::cdf(normal) evaluates erfc in
//     long double and narrows, i.e. it is the correctly rounded double down
//     into the subnormals, 0 from z ~ 27.226 (p ~ 2.5e-324). The device uses
//     Boost's 53-bit rational approximations with exp(-z^2) scaled so the
//     only rounding into the subnormal range is the final one
//     (erfc_boost53.cuh), which reproduces the subnormal p and the p == 0
//     cutoff; erfcMode selects CUDA's erfc or the literal port for measuring
//     the difference (SPA_GPU_LIB.md section 4);
//   * K2's sum skips non-finite terms like sum_arma1 does; Korg and K1 keep
//     them like arma::sum does.
// fp64 throughout (FP64): in plain fp32 the squared denominator of K2
// underflows; the FP32 variant below avoids it by rewriting the terms.
//
// What stays on the host, after run(): the quantile step that turns the SPA
// p-value into seBeta (and can itself withdraw convergence), the p == 0 rule,
// the "%.6E" / "%.1fE%d" string, the Firth decision. spa_gpu_test.cpp applies
// those rules on both sides so the printed strings can be compared.
//
// Thread safety: one Spa per host thread; run() is synchronous.
//
// Precision (config key gpuPrecisionSPA, ../gpu_precision.hpp): CreateArgs::
// precision, FP64 (default, everything above) or FP32. Only modes supports()
// accepts may be passed; create() refuses any other. CONTRACT for any non-fp64
// variant: the pair table in and out keeps its types -- PairIn's Tstat, var1,
// var2, pno arrive as double, and every double field of PairOut (pval, roots,
// p1, p2, m1, q, qinv) is handed back as double (the fp32 value converted),
// with conv / status / reason / niter / nnz / s1 / s2 meaning exactly what
// they mean in fp64. The host post-rules (quantile step, p == 0, the Firth
// decision, the formatting) read those doubles unchanged and do not know which
// mode ran.
//
// FP32 (spa_gpu.cu, "fp32 variant"): same control flow, same block per pair,
// same fused / dynamic / ownSamples / occupancy options. Meant for GPUs whose
// fp64 rate is 1/64 of fp32, so every per-sample and per-reduction operation
// is float: passes A and B (b = XV g, g~ = g - XXVX_inv b, gpos / gneg, m1,
// the carrier sums) read float copies of XV, XXVX_inv and mu made once at
// create(); the Newton and tail passes (exp, expm1, log1p, the K', K'', K
// sums); and the block reductions (float-float). Double: only the per-pair
// scalars -- the sums read out, q, NAmu / NAsigma, the Newton step (t, K1,
// K2, prevJump) and the whole saddlepoint tail (w, v, Ztest, the erfc port,
// the log domain), so a p-value cannot underflow because of float.
// Safeguards:
//   * centred sums: the passes return K' - m1 and K - t m1, and K' - q,
//     zeta q - K are formed in double from them and q - m1 -- no difference
//     of two large float numbers;
//   * only exp(-|g~ t|) <= 1 is formed per term (p, 1 - p, p - mu and
//     log(1 - mu + mu e^x) - mu x via expm1 / log1p), so float never
//     overflows; mu is stored as min(mu, 1 - mu) with g~ negated when
//     mu > 1/2 (each term is invariant under that), so 1 - mu is exact enough
//     for mu near 1;
//   * per-thread sums are Neumaier-compensated floats; products that feed
//     m1 and b are added error-free (TwoProd via fma); block reductions merge
//     the (sum, compensation) pairs with TwoSum, in a fixed order; a sum is
//     read out in double once per pair / pass;
//   * m1 is summed from exactly the stored float (g~, mu), so it is also the
//     centring the passes need;
//   * NAsigma = var2 - sum_c mu (1-mu) g~^2 cancels by up to ~1e3 for very
//     rare markers; that carrier sum is float-float: mu, XV and XXVX_inv
//     are also kept as lo floats (x - float(x)), read for carriers only, so
//     b and the carriers' g~ come out as hi + lo pairs. Without it the root
//     moves by ~1e-4 relative and can cross fp64's Korg overflow threshold.
//     Device memory for FP32: 8 N (mu hi, lo) + 16 traitStride (XV, XXVX_inv
//     hi, lo) bytes per trait, on top of the fp64 tables;
//   * g~ t is formed as g~ th + g~ tl with t = th + tl split into two floats,
//     so a large t is not rounded to float before the product;
//   * Newton tolerance: max(tol, 1e-5, 2^-17 max(|t|, |tnew|)) (kTolFp32,
//     kRelTolFp32), used for both the |dt| test and the prevJump test. The
//     default tol, eps^(1/4) = 1.2e-4, is what decides almost every pair; the
//     relative term takes over only for |t| > ~16 and exists because the
//     float sums resolve a Newton step to about 1e-7 |t|: roots of 1e3 .. 1e5
//     occur for very rare markers (full-N variant, sparse GRM) and would
//     otherwise run to maxiter. maxiter still bounds the loop. The root error
//     this allows does not reach p at first order (zeta q - K is stationary at
//     the root).
//   * fp64's overflow of exp(g~ t) in Korg (g~ t > 709.78), which makes a tail
//     "not a saddle" in SAIGE, is reproduced as Korg = +inf, so the route
//     (Is.SPA) is fp64's.
// The fp64 path's arithmetic is unchanged by the variant (bit for bit).
#pragma once

#include <cstddef>
#include <cstdint>

#include "../gpu_precision.hpp"

namespace saige {
namespace spa_gpu {

struct Spa;   // opaque; defined in spa_gpu.cu

constexpr int PMAX = 8;   // covariate columns per trait, incl. intercept

// Per-trait null-model pieces, host pointers, copied at create().
struct TraitArgs {
    const double* mu       = nullptr;   // N
    const double* XV       = nullptr;   // p x N column-major: sample i's p values at XV + i*p
    const double* XXVX_inv = nullptr;   // N x p column-major: XXVX_inv + j*N + i
    int p = 0;                          // 1..PMAX
    // CreateArgs::ownSamples only: the trait's samples among the N, 1 bit per
    // sample (bit i of word i >> 6, (N + 63) / 64 words); nullptr = all N.
    // A sample outside it has dosage 0, and the caller embeds mu / XV /
    // XXVX_inv with exact zeros there, so it adds exact zeros to every sum.
    const uint64_t* mask   = nullptr;
};

struct CreateArgs {
    int device = 0;
    int N = 0;
    int nTraits = 0;
    const TraitArgs* traits = nullptr;  // nTraits entries, indexed by PairIn::trait
    int maxPairs = 0;                   // per run()
    double tol = 0.0;                   // |dt| convergence; getMarkerPval uses pow(eps, 0.25)
    int maxiter = 1000;
    // Resident blocks. Each needs 16*N bytes of device scratch (g~ and mu,
    // interleaved; the fast variant stores carriers only). 256 blocks keep a
    // V100 busy; at N = 435k that is 1.8 GB, so lower it if memory is tight.
    int blocks = 256;
    // 1 (default): Boost's 53-bit erfc with the underflow deferred to one
    // final rounding, so subnormal p and the p == 0 cutoff are the CPU's
    // (erfc_boost53.cuh). 2: the literal Boost 53-bit port. 0: CUDA libm.
    int erfcMode = 1;
    // Pass fusion (default off). The two Newton solves of a pair (q and qinv)
    // and the two saddlepoint tails read the same stored (g~, mu) pairs; with
    // fusedRoots != 0 one pass over them serves both solves, so a pair makes
    // about half as many passes. Each root still visits exactly the sequence
    // of t the unfused solve visits and each per-thread partial sum is
    // accumulated in the same order, so every output field is bit-identical
    // to fusedRoots = 0; only the memory traffic changes
    // (S2_KERNEL_ROOFLINE.md).
    int fusedRoots = 0;
    // Dynamic pair assignment (default off): blocks take the next pair from a
    // device counter instead of a fixed grid stride, so the blocks in flight
    // always work on a window of consecutive pairs. With the caller's pair
    // table sorted by trait, the trait's XXVX_inv and mu then stay in L2
    // across the blocks that need them. Per-pair arithmetic is unchanged.
    int dynamicPairs = 0;
    // Occupancy hint (default 0 = none): the kernel's __launch_bounds__
    // minimum blocks per SM, 1..4. With no hint the unfused kernel compiles to
    // 128 registers on a V100 (2 blocks of 8 warps per SM) and the fused one to
    // 194 (1 block); an explicit 1 lets the unfused one grow to 162 (1 block);
    // 2 / 3 / 4 cap the registers at 128 / 80 / 64 and spill the rest to local
    // memory for 16 / 24 / 32 resident warps. Arithmetic unchanged.
    int minBlocksPerSM = 0;
    // Per-trait sample sets (default off). The N samples are the union of the
    // traits' sample lists; each trait carries a mask (TraitArgs::mask) and
    // each pair its own 4-entry dosage table, filled by the caller in
    // pairLut() next to in() (a trait's imputed value and flip can differ from
    // the union column's). Off: the slot's table, every sample counts.
    int ownSamples = 0;
    // Arithmetic of the kernel (gpuPrecisionSPA); see the contract above.
    saige::gpu2::Prec precision = saige::gpu2::Prec::FP64;
};

// True when create() accepts CreateArgs::precision = t_p. main.cpp asks first
// and stops the run on false ("useGPU: SPA precision <mode> is not implemented
// yet").
bool supports(saige::gpu2::Prec t_p);

// Device-resident genotype source for the slots PairIn::slot names. Exactly
// one of the two forms is set.
struct Geno {
    // PLINK 2-bit codes, 4 samples per byte, sample i in bits 2*(i&3) of byte
    // i>>2; slot s at packed + s*bpv, its code->dosage table at lut + 4*s
    // (the reducer's layout, gpu_step2.hpp). A carrier is a sample whose
    // table entry is != 0, as fillOneMarkerFusedDense_ts defines it.
    const unsigned char* packed = nullptr;
    std::size_t          bpv    = 0;
    const double*        lut    = nullptr;
    // Dense dosages, N doubles per slot, slot s at dense + s*ld.
    const double*        dense  = nullptr;
    std::size_t          ld     = 0;
};

struct PairIn {
    int    slot;     // index into Geno
    int    trait;    // index into CreateArgs::traits
    int    fast;     // 1: SPA_fast variant; 0: full-N SPA; -1: decide on the
                     //    device by the CPU's rule |{g==0}|/N >= 0.5
    int    logp;     // 1: pno is log(p) and pval comes back as log(p)
    double Tstat, var1, var2, pno;
};

// Status bits (PairOut::status).
enum : unsigned {
    ST_FAST         = 1u << 0,   // the fast (carriers) variant ran
    ST_LOGP         = 1u << 1,   // log domain
    ST_ROOT1_INF    = 1u << 2,   // q outside (gneg, gpos): root1 = +inf, niter 0, "converged"
    ST_ROOT2_INF    = 1u << 3,   // same for qinv
    ST_ROOT1_FAIL   = 1u << 4,   // Newton did not converge (reason1 says why)
    ST_ROOT2_FAIL   = 1u << 5,
    ST_SADDLE1_FAIL = 1u << 6,   // tail 1 not a saddle: p1 = pno/2, conv withdrawn
    ST_SADDLE2_FAIL = 1u << 7,
    ST_PNO          = 1u << 8    // a root failed: pval = pno, conv = false
};
// Root-find failure reasons (PairOut::reason1/2).
enum : int {
    RS_OK        = 0,
    RS_K2_GUARD  = 1,   // full variant only: K2 non-finite or |K2| < 1e-15
    RS_TNEW      = 2,   // tnew non-finite (full) / NaN (fast)
    RS_MAXITER   = 3
};

struct PairOut {
    double   pval;       // SPA p-value (log when logp), or pno when not converged
    int      conv;       // the wrapper's Isconverge (spa.cpp SPA / SPA_fast)
    unsigned status;     // ST_* bits
    int      reason1, reason2;
    int      niter1, niter2;
    double   root1, root2;
    double   p1, p2;     // the two tail terms as combined (after the pno/2 substitution)
    double   m1, q, qinv;
    int      nnz;        // carriers (g != 0)
    int      s1, s2;     // isSaddle per tail; -1 when the roots did not both converge
};

// nullptr on any failure (no device, allocation refused, bad args); lastError()
// says why. The caller then keeps SPA on the CPU.
Spa* create(const CreateArgs& t_args);
void destroy(Spa* t_s);

// Pinned tables of maxPairs entries, owned by the Spa. Fill in(), call run(),
// read out().
PairIn*  in(Spa* t_s);
PairOut* out(Spa* t_s);
// ownSamples: pinned maxPairs x 4 doubles, pair k's dosage table at 4k.
// nullptr without ownSamples.
double*  pairLut(Spa* t_s);

// Solve pairs [0, t_nPairs) of in() against t_geno; results in out().
// Synchronous. false on a CUDA error (results undefined; lastError()).
bool run(Spa* t_s, const Geno& t_geno, int t_nPairs);

// Convenience for callers that hold genotypes on the host (tests, a standalone
// step 2.2): copy nSlots packed rows (+ 4-entry tables) or dense columns into a
// buffer the Spa owns and return the Geno view for them. Reallocates as needed;
// the previous view is invalid after the call.
bool uploadPacked(Spa* t_s, const unsigned char* t_packed, std::size_t t_bpv,
                  const double* t_lut, int t_nSlots, Geno* t_view);
bool uploadDense(Spa* t_s, const double* t_g, int t_nSlots, Geno* t_view);

// Cumulative device seconds since create() -- kernel, H2D of the pair table,
// D2H of the results -- and pairs solved.
void timings(const Spa* t_s, double* t_kernel, double* t_h2d, double* t_d2h,
             long long* t_pairs);
std::size_t deviceBytes(const Spa* t_s);
const char* lastError();

// Instrumentation for the tail study: evaluate, on the device, the scaled
// port (erfcMode 1), the literal port (erfcMode 2) and CUDA's erfc (erfcMode
// 0) at n points. Needs no Spa.
bool debugErfc(int t_device, const double* t_z, int t_n, double* t_scaledPort,
               double* t_literalPort, double* t_cuda);

}  // namespace spa_gpu
}  // namespace saige
