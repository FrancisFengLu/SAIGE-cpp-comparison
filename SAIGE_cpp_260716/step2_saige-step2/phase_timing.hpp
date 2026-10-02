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
    // ---- finalize sub-slots (SPA_CLAIM_GAPS gap 1) ----
    S_FIN_STAT,   // per-(marker,trait) stat copy: MAC/altFreq/flip + O.altFreq.. writes
    S_FIN_BATCH,  // useBatch branch: copy Beta/se/Tstat/var + two std::string p-value copies
    S_FIN_AF,     // AF_case / AF_ctrl: the two index gathers (or pcSumPair)
    S_FIN_SLOT,   // the rest of the slot writes (O.Beta..varT, AF/N_case/N_ctrl, isSPAConverge); nests S_FIN_AF
    // ---- fallback sub-slots (S2_REMAINDER): what the scalar path does besides the root-find ----
    S_FB_IDX,     // main.cpp: {i : g[i]==0} / complement index build (once per flagged marker)
    S_FB_ALLOC,   // getMarkerPval: constructing gNB/gNA/muNB/muNA (N-sized; mmap under the 64 KB mallopt)
    S_SPA_PREP,   // getMarkerPval: prepare_spa_inputs (gtilde, m1, the four subsets, NAmu/NAsigma)
    S_SPA_ROOT,   // getMarkerPval: the SPA_fast / SPA call (two root-finds + two saddle probabilities)
    S_SPA_GPOS,   // spa_binary.cpp: the gpos/gneg accu(g.elem(find(..))) at the top of each root-find
    S_PARWALL,    // main.cpp, master thread only: wall of the block-parallel region per chunk
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
