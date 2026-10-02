// See phase_timing.hpp.  Whole file is empty unless -DMT_PHASE_TIMING.
#include "phase_timing.hpp"

#ifdef MT_PHASE_TIMING
#include <iostream>
#include <iomanip>

namespace PT {

static const char* kSlotName[S_NSLOT] = {
    "read+QC+impute", "batch GEMM", "finalize(total)", "gate", "fallback(total)",
    "  fb: score recompute", "  fb: SPA", "  fb: ER", "  fb: Firth fit",
    "  fin: stat copy", "  fin: batch copy", "  fin: AF gather", "  fin: slot writes(+AF)",
    "  fb: index build", "  fb: 4 vec allocs", "  spa: prep inputs", "  spa: root+saddle",
    "  spa: gpos/gneg", "par-region wall(master)"
};
static const char* kCntName[C_NSLOT] = {
    "pairs reaching gate", "gate needSPA", "gate needFirth", "gate needFast",
    "binary pairs", "binary pairs MAC<=ERcut"
};

void report(double tWrite, int nThreads) {
    double T[S_NSLOT] = {0};
    long long N[S_NSLOT] = {0};
    long long C[C_NSLOT] = {0};
    for (int i = 0; i < 512; i++) {
        for (int s = 0; s < S_NSLOT; s++) { T[s] += g_acc[i].t[s]; N[s] += g_acc[i].n[s]; }
        for (int c = 0; c < C_NSLOT; c++) C[c] += g_acc[i].c[c];
    }
    const double nt = (nThreads > 0 ? (double)nThreads : 1.0);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  [phase] threads " << nThreads
              << "   (cpu-s = summed over threads; wall-s = cpu-s / threads)\n";
    for (int s = 0; s < S_NSLOT; s++) {
        std::cout << "  [phase] " << std::left << std::setw(22) << kSlotName[s]
                  << std::right << " cpu-s " << std::setw(12) << T[s]
                  << "  wall-s " << std::setw(10) << (T[s] / nt)
                  << "  calls " << N[s] << "\n";
    }
    // Exclusive decomposition of the block-parallel region.
    const double finalExcl = T[S_FINAL] - T[S_GATE] - T[S_FALL];
    const double fbExcl    = T[S_FALL] - T[S_SCORE] - T[S_SPA] - T[S_ER] - T[S_FIRTH];
    const double phase1cpu = T[S_READ] + T[S_GEMM] + T[S_GATE] + finalExcl;
    const double phase2cpu = T[S_FALL];
    std::cout << "  [phase] finalize excl(gate,fallback) cpu-s " << finalExcl << "\n";
    {
        const double slotExcl = T[S_FIN_SLOT] - T[S_FIN_AF];
        const double unattr = finalExcl - T[S_FIN_STAT] - T[S_FIN_BATCH] - T[S_FIN_SLOT];
        std::cout << "  [phase]   finalize excl = stat copy " << T[S_FIN_STAT]
                  << " + batch copy " << T[S_FIN_BATCH]
                  << " + AF gather " << T[S_FIN_AF]
                  << " + slot writes excl AF " << slotExcl
                  << " + unattributed(loop/timer overhead) " << unattr << "\n";
    }
    std::cout << "  [phase] fallback excl(score,SPA,ER,Firth) cpu-s " << fbExcl << "\n";
    {
        // S_FB_IDX sits inside S_FALL but outside getMarkerPval; S_FB_ALLOC is
        // inside getMarkerPval but outside S_SCORE / S_SPA. The remainder is the
        // call plumbing, StdStat, the q/qinv branch, the post-SPA qnorm +
        // sprintf, and the munmap of the four vectors at function exit.
        const double fbRest = fbExcl - T[S_FB_IDX] - T[S_FB_ALLOC];
        std::cout << "  [phase]   fallback excl = index build " << T[S_FB_IDX]
                  << " + 4 vec allocs " << T[S_FB_ALLOC]
                  << " + rest(call/StdStat/qnorm/sprintf/munmap) " << fbRest << "\n";
        const double spaRest = T[S_SPA] - T[S_SPA_PREP] - T[S_SPA_ROOT];
        const double rootExcl = T[S_SPA_ROOT] - T[S_SPA_GPOS];
        std::cout << "  [phase]   SPA = prep inputs " << T[S_SPA_PREP]
                  << " + root+saddle " << T[S_SPA_ROOT]
                  << " (of which gpos/gneg " << T[S_SPA_GPOS]
                  << ", Newton+saddle " << rootExcl << ")"
                  << " + rest(q/qinv, post qnorm+sprintf) " << spaRest << "\n";
        const double busy = T[S_READ] + T[S_GEMM] + T[S_FINAL];
        const double idle = nt * T[S_PARWALL] - busy;
        std::cout << "  [phase] par-region wall " << T[S_PARWALL] << " s x " << nThreads
                  << " threads = " << (nt * T[S_PARWALL]) << " thread-s; busy(read+GEMM+finalize) "
                  << busy << "; idle(barrier imbalance) " << idle
                  << " thread-s = " << (idle / nt) << " wall-s\n";
    }
    std::cout << "  [phase] PHASE1 cpu-s " << phase1cpu
              << "  wall-s " << (phase1cpu / nt) << " (+ writer " << tWrite << " s wall)\n";
    std::cout << "  [phase] PHASE2 cpu-s " << phase2cpu
              << "  wall-s " << (phase2cpu / nt) << "\n";
    for (int c = 0; c < C_NSLOT; c++)
        std::cout << "  [phase] count " << std::left << std::setw(26) << kCntName[c]
                  << std::right << " " << C[c] << "\n";
    std::cout.unsetf(std::ios::floatfield);
}

}  // namespace PT
#endif
