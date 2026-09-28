#ifndef SAIGE_PHASE_TIMING_HPP
#define SAIGE_PHASE_TIMING_HPP
// ---------------------------------------------------------------------------
// Compile-time-gated stage timers for the multi-trait single-variant loop.
//
// Purpose: split step 2's multi-trait work into the dense "fast-path scan"
// (BED read + QC + impute, the batched score GEMMs, the gate decisions, the
// output writer) and the irregular "slow-path correction" (the scalar
// getMarkerPval fallback: score recompute, SPA, ER, Firth), so the share of
// each in wall clock can be measured rather than inferred by differencing two
// runs with different configurations.
//
// NOTHING here is compiled unless -DMT_PHASE_TIMING is passed.  The shipped
// binary (`make`) expands every macro below to `((void)0)` and links no extra
// symbol; `make saige-step2.phase` builds the instrumented copy under a
// different name.  Same idiom as the existing MTFOLD_PROF / MTVEC_PROF blocks.
//
// Accumulation is per OpenMP thread into a padded slot array, so the timers
// add one omp_get_wtime() pair per region and no synchronisation.  The
// multi-trait loop's parallel dimension is the block and one thread owns a
// block end to end, so thread-seconds divided by the thread count is the
// wall-clock share of that stage.
// ---------------------------------------------------------------------------

#ifdef MT_PHASE_TIMING

#include <omp.h>
#include <cstdio>
#include <string>

namespace PT {

// Time slots.  S_FINAL nests S_GATE and S_FALL; S_FALL nests S_SCORE, S_SPA,
// S_ER and S_FIRTH.  Exclusive times are differences taken in the report.
enum Slot {
    S_READ = 0,   // BED read + QC + impute + variance ratio, per block
    S_GEMM,       // scoreTestBatchMT calls (the batched score reduction)
    S_FINAL,      // the whole per-(marker,trait) finalize loop
    S_GATE,       // gate decision: fastRecomputeSameCtx probe + needSPA/Firth/Fast
    S_FALL,       // the scalar fallback branch (phase 2)
    S_SCORE,      // inside getMarkerPval: the score test recompute
    S_SPA,        // inside getMarkerPval: the SPA block
    S_ER,         // inside getMarkerPval: the ER block
    S_FIRTH,      // inside getMarkerPval: the Firth fit
    S_NSLOT
};

// Count slots.
enum Cnt {
    C_GATED = 0,  // (marker,trait) pairs that reached the batch gate
    C_NEEDSPA,    // ... of those, gated to the fallback because needSPA
    C_NEEDFIRTH,  // ... because needFirth
    C_NEEDFAST,   // ... because needFast
    C_PAIR_BIN,   // binary (marker,trait) pairs finalized
    C_ER_LOWMAC,  // binary pairs with MAC <= MACCutoffforER (never batched)
    C_NSLOT
};

struct alignas(128) Acc {
    double    t[S_NSLOT];
    long long n[S_NSLOT];
    long long c[C_NSLOT];
    char      pad[128];
};

inline Acc g_acc[512];

inline void add(int s, double t0) {
    Acc& a = g_acc[omp_get_thread_num()];
    a.t[s] += omp_get_wtime() - t0;
    a.n[s] += 1;
}
inline void cnt(int c) { g_acc[omp_get_thread_num()].c[c] += 1; }
inline void cntif(int c, bool p) { if (p) g_acc[omp_get_thread_num()].c[c] += 1; }

struct Scope {
    int s; double t0;
    explicit Scope(int slot) : s(slot), t0(omp_get_wtime()) {}
    ~Scope() { add(s, t0); }
};

// Sums over threads and prints one "[phase]" line per slot.  `tWrite` is the
// already-measured wall time of the output writer, passed in so the report can
// state phase 1 and phase 2 in the same units.
void report(double tWrite, int nThreads);

}  // namespace PT

#define PT_T0(v)        const double v = omp_get_wtime()
#define PT_ADD(s, v)    PT::add(PT::s, (v))
#define PT_SCOPE(s)     PT::Scope pt_scope_guard_(PT::s)
#define PT_CNT(c)       PT::cnt(PT::c)
#define PT_CNTIF(c, p)  PT::cntif(PT::c, (p))
#define PT_REPORT(w, n) PT::report((w), (n))

#else   // !MT_PHASE_TIMING -- the shipped binary

#define PT_T0(v)        ((void)0)
#define PT_ADD(s, v)    ((void)0)
#define PT_SCOPE(s)     ((void)0)
#define PT_CNT(c)       ((void)0)
#define PT_CNTIF(c, p)  ((void)0)
#define PT_REPORT(w, n) ((void)0)

#endif  // MT_PHASE_TIMING
#endif  // SAIGE_PHASE_TIMING_HPP
