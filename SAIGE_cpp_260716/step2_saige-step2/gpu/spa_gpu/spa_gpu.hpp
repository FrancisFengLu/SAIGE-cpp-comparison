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
//   * the normal tail is Boost's own 53-bit erfc (erf.hpp, Boost 1.85),
//     ported operation for operation including the error-compensated
//     exp(-z^2) and the "z >= 28 is zero" cutoff, so the smallest p the CPU
//     can print and the point where it becomes 0 (or -inf in the log domain)
//     are the same on both sides; erfcMode = 0 switches to CUDA's erfc for
//     measuring the difference (SPA_GPU_LIB.md section "The tail");
//   * K2's sum skips non-finite terms like sum_arma1 does; Korg and K1 keep
//     them like arma::sum does.
// fp64 throughout: in fp32 the squared denominator of K2 underflows.
//
// What stays on the host, after run(): the quantile step that turns the SPA
// p-value into seBeta (and can itself withdraw convergence), the p == 0 rule,
// the "%.6E" / "%.1fE%d" string, the Firth decision. spa_gpu_test.cpp applies
// those rules on both sides so the printed strings can be compared.
//
// Thread safety: one Spa per host thread; run() is synchronous.
#pragma once

#include <cstddef>
#include <cstdint>

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
    // 1: Boost's 53-bit erfc, ported (default). 0: CUDA libm erfc.
    int erfcMode = 1;
};

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

// Instrumentation for the tail study: evaluate, on the device, the ported
// Boost erfc and CUDA's erfc at n points. Needs no Spa.
bool debugErfc(int t_device, const double* t_z, int t_n, double* t_boostPort,
               double* t_cuda);

}  // namespace spa_gpu
}  // namespace saige
