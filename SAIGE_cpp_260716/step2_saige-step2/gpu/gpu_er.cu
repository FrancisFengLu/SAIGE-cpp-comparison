// gpu_er.cu — the exact test (ER) for binary traits, one thread per (marker,
// trait) pair. Contract and provenance: gpu_er.hpp. Reference: er_binary.cpp
// (SKATExactBin_Work -> SKATExactBin_ComputProb_New -> HyperGeo; SKAT_Exact ->
// ComputeExact::Init / Run). Every floating-point statement below is the CPU's,
// in the CPU's order. Compiled with --fmad=false, so nothing is fused unless
// written as fma(): the CPU build (g++ -O3 -march=native, whose default for C++
// is -ffp-contract=fast even with -std=c++17) fused exactly five expressions in
// er_binary.o -- lch + lweight*j (both HyperGeo tables), stat += one*one
// (CalTestStat / _INV) and pval - pval_same/2 -- and those five are fma() here,
// read off the disassembly (S2_ER_GPU.md §2). The other fma() calls are the
// glibc exp / log ports' (gpu_er_glibc.cuh), the fused operations of the CPU's libm.

#include "gpu_er.hpp"
#include "gpu_er_glibc.cuh"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace saige {
namespace gpu2 {

namespace {

thread_local std::string lastErrEr;
#define CKE(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    lastErrEr = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)

constexpr int KM = kErMaxCarriers;
constexpr int NBIN = 10;            // SKATExactBin_ComputeProb_Group's ngroup1

// The 2^k case/control assignments of the carriers in ComputeExact::Run's
// order, calling f(j, raw, stat) for each: j = the number of cases, raw = the
// un-normalised Fisher weight (CalFisherProb / CalFisherProb_INV), stat = the
// statistic (CalTestStat / CalTestStat_INV, m = 1), computed only when
// WANT_STAT. Strata j <= k/2 + 1 enumerate the j-subsets of the cases, the
// others the (k-j)-subsets of the controls, both lexicographically -- the
// order SKAT_Exact_Recurse / SKAT_Exact_Recurse_INV visit them in.
// The CPU folds each subset left to right from 1 (pprod) and from the
// all-control (all-case) statistic; consecutive subsets in lexicographic order
// share a prefix, and the folded value of a prefix does not depend on what
// follows it, so pr[q] / ts[q] (the folds over a[0..q-1]) are kept and only the
// changed suffix is refolded -- the same operations on the same operands.
template <bool WANT_STAT, typename F>
__device__ __forceinline__ void forEachConfig(int k, const double* Z0, const double* Z1, const double* odds,
                                              double tsZ0, double tsZ1, double pprod, F f)
{
    int a[KM];
    double pr[KM + 1], ts[KM + 1];
    for (int j = 0; j <= k; j++) {
        const bool inv = !(j <= k / 2 + 1);
        const int r = inv ? (k - j) : j;
        for (int q = 0; q < r; q++) a[q] = q;
        pr[0] = inv ? pprod : 1.0;
        ts[0] = inv ? tsZ1 : tsZ0;
        int from = 0;   // pr / ts are valid up to index `from`
        while (true) {
            for (int q = from; q < r; q++) {
                const int l = a[q];
                if (!inv) {
                    pr[q + 1] = pr[q] * odds[l];
                    if (WANT_STAT) ts[q + 1] = ts[q] + (Z1[l] - Z0[l]);
                } else {
                    pr[q + 1] = pr[q] / odds[l];
                    if (WANT_STAT) ts[q + 1] = ts[q] + (Z0[l] - Z1[l]);
                }
            }
            const double stat = WANT_STAT ? fma(ts[r], ts[r], 0.0) : 0.0;   // stat += one*one, contracted
            f(j, pr[r], stat);
            // next r-subset of [0, k) in lexicographic order
            int q = r - 1;
            while (q >= 0 && a[q] == k - r + q) q--;
            if (q < 0) break;
            a[q]++;
            for (int u = q + 1; u < r; u++) a[u] = a[u - 1] + 1;
            from = q;
        }
    }
}

__global__ void __launch_bounds__(32)
er_pairs(const double* __restrict__ mu, const long long* __restrict__ muOff,
         const int* __restrict__ nOf, const int* __restrict__ ncaseOf,
         const double* __restrict__ LT,
         const ErPairIn* __restrict__ in, int nPairs,
         const uint32_t* __restrict__ carIdx, const double* __restrict__ carG,
         const unsigned char* __restrict__ carCase,
         ErPairOut* __restrict__ out)
{
    for (int pq = blockIdx.x * blockDim.x + threadIdx.x; pq < nPairs; pq += gridDim.x * blockDim.x) {
        const ErPairIn P = in[pq];
        const int t = P.trait, k = P.k;
        const int n = nOf[t], ncase = ncaseOf[t];
        const double* m = mu + muOff[t];
        int idx[KM];
        double g[KM], p1[KM];
        bool cs[KM];
        for (int i = 0; i < k; i++) {
            idx[i] = (int)carIdx[P.off + i];
            g[i] = carG[P.off + i];
            cs[i] = carCase[P.off + i] != 0;
            p1[i] = m[idx[i]];
        }

        // ---- SKATExactBin_ComputeProb_Group ----
        // p2temp = mean(pi1(idxCompVec)): armadillo's accumulate, two
        // interleaved accumulators over the non-carriers in index order.
        double p2temp;
        {
            double acc1 = 0.0, acc2 = 0.0;
            long long pos = 0;
            int u = 0;
            for (int c = 0; c <= k; c++) {
                const int end = (c < k) ? idx[c] : n;
                for (; u < end; u++, pos++) {
                    if (pos & 1) acc2 = acc2 + m[u];
                    else         acc1 = acc1 + m[u];
                }
                u = end + 1;
            }
            p2temp = (acc1 + acc2) / (double)(n - k);
        }
        double weight[NBIN + 1];
        int group[NBIN + 1];
        int ng = 0;
        for (int b = 0; b < NBIN; b++) {
            const double a1 = double(b) / NBIN;
            const double a2 = double(b + 1) / NBIN;
            double acc1 = 0.0, acc2 = 0.0;
            int cnt = 0;
            for (int i = 0; i < k; i++) {
                const double pc = (p1[i] >= 1) ? 0.999 : p1[i];
                const bool inBin = (pc >= a1) && ((b + 1) < NBIN ? (pc < a2) : (pc <= a2));
                if (!inBin) continue;
                if (cnt & 1) acc2 = acc2 + pc;
                else         acc1 = acc1 + pc;
                cnt++;
            }
            if (cnt > 0) {
                const double p1temp = (acc1 + acc2) / (double)cnt;
                weight[ng] = p1temp / (1 - p1temp);
                group[ng] = cnt;
                ng++;
            }
        }
        const double p2oddtemp = p2temp / (1 - p2temp);
        weight[ng] = p2oddtemp;
        for (int i = 0; i <= ng; i++) weight[i] = weight[i] / p2oddtemp;
        group[ng] = n - k;
        const int ngroup = ng + 1;

        // ---- HyperGeo::Run ----
        double lw[NBIN + 1];
        for (int i = 0; i < ngroup; i++) lw[i] = glx::log_(weight[i]);
        // carrier groups: temp[j] = lCombinations(group, j) + lweight * j
        double tabC[KM + NBIN];
        int tabOff[NBIN];
        {
            int o = 0;
            for (int i = 0; i < ngroup - 1; i++) {
                tabOff[i] = o;
                for (int j = 0; j <= group[i]; j++) {
                    double lch = 0.0;
                    if (!(j > group[i])) {
                        int nn = group[i];
                        for (int d = 1; d <= j; ++d) { lch = lch + LT[nn--]; lch = lch - LT[d]; }
                    }
                    tabC[o + j] = fma((double)j, lw[i], lch);   // gcc contracted lch + lweight*j
                }
                o += group[i] + 1;
            }
        }
        // last group: temp[j] = lCombinations(n - k, ncase - j) + lweight * (ncase - j).
        // lCombinations(nn, kk) is the running sum after kk steps of one loop
        // that does not depend on kk, so all k + 1 values come from one pass.
        double tabL[KM + 1];
        double mref = 0.0;
        {
            double lch[KM + 1];
            for (int j = 0; j <= k; j++) lch[j] = 0.0;
            const int nn0 = group[ngroup - 1];
            const int dmax = (ncase < nn0) ? ncase : nn0;
            double r = 0.0;
            int nn = nn0;
            for (int d = 1; d <= dmax; ++d) {
                r = r + LT[nn--];
                r = r - LT[d];
                const int j = ncase - d;
                if (j >= 0 && j <= k) lch[j] = r;
            }
            for (int j = 0; j <= k; j++) {
                tabL[j] = fma((double)ncase - j, lw[ngroup - 1], lch[j]);   // contracted, as above
                mref = (tabL[j] > mref) ? tabL[j] : mref;
            }
        }
        // Recursive(0, 0, 0), depth first, as an explicit stack
        double kprob[KM + 1];
        for (int j = 0; j <= k; j++) kprob[j] = 0.0;
        {
            const int G = ngroup - 1;
            int sel[NBIN + 1];
            double ps[NBIN + 1];
            int nc[NBIN + 1];
            int lev = 0;
            ps[0] = 0.0; nc[0] = 0; sel[0] = -1;
            while (lev >= 0) {
                if (lev == G) {
                    const double lp = ps[G] + tabL[nc[G]];
                    kprob[nc[G]] = kprob[nc[G]] + glx::exp_(lp - mref);
                    lev--;
                    continue;
                }
                sel[lev]++;
                if (sel[lev] > group[lev]) { lev--; continue; }
                const int i = sel[lev];
                if (nc[lev] + i <= ncase) {
                    ps[lev + 1] = ps[lev] + tabC[tabOff[lev] + i];
                    nc[lev + 1] = nc[lev] + i;
                    lev++;
                    sel[lev] = -1;
                }
            }
        }
        double prob[KM + 1];
        {
            double sum1 = 0.0;
            for (int i = 0; i <= k; i++) sum1 = sum1 + kprob[i];
            for (int i = 0; i <= k; i++) prob[i] = kprob[i] / sum1;
        }

        // ---- SKATExactBin_Work + ComputeExact ----
        double Z0[KM], Z1[KM], odds[KM];
        for (int i = 0; i < k; i++) {
            Z0[i] = g[i] * (-p1[i]);
            Z1[i] = g[i] * (1 - p1[i]);
            odds[i] = p1[i] / (1 - p1[i]);
        }
        double pprod = 1.0, tsZ0 = 0.0, tsZ1 = 0.0;
        for (int i = 0; i < k; i++) pprod = pprod * odds[i];
        for (int i = 0; i < k; i++) { tsZ0 = tsZ0 + Z0[i]; tsZ1 = tsZ1 + Z1[i]; }
        double Q;
        {
            double tq = tsZ0;
            for (int i = 0; i < k; i++) if (cs[i]) tq = tq + (Z1[i] - Z0[i]);
            Q = fma(tq, tq, 0.0);
        }
        double denomi[KM + 1];
        for (int j = 0; j <= k; j++) denomi[j] = 0.0;
        forEachConfig<false>(k, Z0, Z1, odds, tsZ0, tsZ1, pprod,
                             [&](int j, double raw, double) { denomi[j] = denomi[j] + raw; });
        double total = 0.0;
        forEachConfig<false>(k, Z0, Z1, odds, tsZ0, tsZ1, pprod,
                             [&](int j, double raw, double) { total = total + raw / denomi[j] * prob[j]; });
        double n_num = 0.0, n_same = 0.0;
        const double eps = 1e-6;
        forEachConfig<true>(k, Z0, Z1, odds, tsZ0, tsZ1, pprod,
                            [&](int j, double raw, double stat) {
                                const double f = raw / denomi[j] * prob[j] / total;
                                double temp1 = Q - stat;
                                if (fabs(temp1) <= eps) temp1 = 0;
                                if (temp1 <= 0) {
                                    n_num = n_num + f;
                                    if (temp1 == 0) n_same = n_same + f;
                                }
                            });
        out[pq].pval = fma(-n_same, 0.5, n_num);   // pval - pval_same / 2: gcc emits vfnmadd with 0.5
    }
}

__global__ void er_math(const double* __restrict__ x, long long n, double* __restrict__ e, double* __restrict__ l)
{
    for (long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x; i < n; i += (long long)gridDim.x * blockDim.x) {
        e[i] = glx::exp_(x[i]);
        l[i] = glx::log_(x[i]);
    }
}

template <typename T>
bool growDev(T** p, std::size_t* cap, std::size_t need)
{
    if (need <= *cap) return true;
    if (*p) cudaFree(*p);
    *p = nullptr; *cap = 0;
    const std::size_t c = need + need / 2;
    if (cudaMalloc(p, c * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
    *cap = c;
    return true;
}

}  // namespace

struct Er {
    int device = 0;
    int nTraits = 0;
    double* dMu = nullptr;
    long long* dMuOff = nullptr;
    int* dN = nullptr;
    int* dNcase = nullptr;
    double* dLT = nullptr;
    ErPairIn* dIn = nullptr;      std::size_t capIn = 0;
    ErPairOut* dOut = nullptr;    std::size_t capOut = 0;
    uint32_t* dIdx = nullptr;     std::size_t capIdx = 0;
    double* dG = nullptr;         std::size_t capG = 0;
    unsigned char* dCase = nullptr; std::size_t capCase = 0;
    cudaStream_t st = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    double tKernel = 0.0;
    long long nPairsDone = 0;
};

Er* erCreate(const ErCreateArgs& a)
{
    if (a.nTraits <= 0 || a.traits == nullptr || a.logTable == nullptr || a.maxN <= 0) {
        lastErrEr = "erCreate: bad arguments";
        return nullptr;
    }
    if (cudaSetDevice(a.device) != cudaSuccess) { lastErrEr = "cudaSetDevice failed"; return nullptr; }
    Er* s = new Er();
    s->device = a.device;
    s->nTraits = a.nTraits;
    auto fail = [&](const char* w) -> Er* { lastErrEr = w; erDestroy(s); return nullptr; };
    std::vector<long long> off(a.nTraits);
    std::vector<int> nv(a.nTraits), cv(a.nTraits);
    long long tot = 0;
    for (int t = 0; t < a.nTraits; t++) {
        if (a.traits[t].mu == nullptr || a.traits[t].n <= 0 || a.traits[t].n > a.maxN)
            return fail("erCreate: a trait's mu / n is unusable");
        off[t] = tot; tot += a.traits[t].n;
        nv[t] = a.traits[t].n; cv[t] = a.traits[t].ncase;
    }
    if (cudaMalloc(&s->dMu, (std::size_t)tot * sizeof(double)) != cudaSuccess) return fail("cudaMalloc mu");
    for (int t = 0; t < a.nTraits; t++)
        if (cudaMemcpy(s->dMu + off[t], a.traits[t].mu, (std::size_t)nv[t] * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            return fail("cudaMemcpy mu");
    if (cudaMalloc(&s->dMuOff, a.nTraits * sizeof(long long)) != cudaSuccess ||
        cudaMalloc(&s->dN, a.nTraits * sizeof(int)) != cudaSuccess ||
        cudaMalloc(&s->dNcase, a.nTraits * sizeof(int)) != cudaSuccess ||
        cudaMalloc(&s->dLT, (std::size_t)(a.maxN + 1) * sizeof(double)) != cudaSuccess)
        return fail("cudaMalloc trait tables");
    if (cudaMemcpy(s->dMuOff, off.data(), a.nTraits * sizeof(long long), cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(s->dN, nv.data(), a.nTraits * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(s->dNcase, cv.data(), a.nTraits * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(s->dLT, a.logTable, (std::size_t)(a.maxN + 1) * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return fail("cudaMemcpy trait tables");
    if (cudaStreamCreateWithFlags(&s->st, cudaStreamNonBlocking) != cudaSuccess) return fail("cudaStreamCreate");
    if (cudaEventCreate(&s->e0) != cudaSuccess || cudaEventCreate(&s->e1) != cudaSuccess) return fail("cudaEventCreate");
    return s;
}

void erDestroy(Er* s)
{
    if (!s) return;
    cudaSetDevice(s->device);
    if (s->e0) cudaEventDestroy(s->e0);
    if (s->e1) cudaEventDestroy(s->e1);
    if (s->st) cudaStreamDestroy(s->st);
    cudaFree(s->dMu); cudaFree(s->dMuOff); cudaFree(s->dN); cudaFree(s->dNcase); cudaFree(s->dLT);
    cudaFree(s->dIn); cudaFree(s->dOut); cudaFree(s->dIdx); cudaFree(s->dG); cudaFree(s->dCase);
    delete s;
}

bool erRun(Er* s, const ErPairIn* in, int nPairs, const uint32_t* carIdx, const double* carG,
           const unsigned char* carCase, long long nCar, ErPairOut* out)
{
    if (!s) { lastErrEr = "erRun: no context"; return false; }
    if (nPairs <= 0) return true;
    for (int q = 0; q < nPairs; q++)
        if (in[q].k < 1 || in[q].k > kErMaxCarriers || in[q].trait < 0 || in[q].trait >= s->nTraits ||
            in[q].off < 0 || in[q].off + in[q].k > nCar) { lastErrEr = "erRun: bad pair"; return false; }
    CKE(cudaSetDevice(s->device));
    const std::size_t nc = (std::size_t)(nCar > 0 ? nCar : 1);
    if (!growDev(&s->dIn, &s->capIn, (std::size_t)nPairs) || !growDev(&s->dOut, &s->capOut, (std::size_t)nPairs) ||
        !growDev(&s->dIdx, &s->capIdx, nc) || !growDev(&s->dG, &s->capG, nc) || !growDev(&s->dCase, &s->capCase, nc)) {
        lastErrEr = "erRun: cudaMalloc";
        return false;
    }
    CKE(cudaMemcpyAsync(s->dIn, in, (std::size_t)nPairs * sizeof(ErPairIn), cudaMemcpyHostToDevice, s->st));
    CKE(cudaMemcpyAsync(s->dIdx, carIdx, (std::size_t)nCar * sizeof(uint32_t), cudaMemcpyHostToDevice, s->st));
    CKE(cudaMemcpyAsync(s->dG, carG, (std::size_t)nCar * sizeof(double), cudaMemcpyHostToDevice, s->st));
    CKE(cudaMemcpyAsync(s->dCase, carCase, (std::size_t)nCar, cudaMemcpyHostToDevice, s->st));
    CKE(cudaEventRecord(s->e0, s->st));
    const int nt = 32;   // few pairs per call: spread them over the SMs
    const int nb = (nPairs + nt - 1) / nt;
    er_pairs<<<nb, nt, 0, s->st>>>(s->dMu, s->dMuOff, s->dN, s->dNcase, s->dLT, s->dIn, nPairs,
                                   s->dIdx, s->dG, s->dCase, s->dOut);
    CKE(cudaGetLastError());
    CKE(cudaEventRecord(s->e1, s->st));
    CKE(cudaMemcpyAsync(out, s->dOut, (std::size_t)nPairs * sizeof(ErPairOut), cudaMemcpyDeviceToHost, s->st));
    CKE(cudaStreamSynchronize(s->st));
    float ms = 0.f;
    if (cudaEventElapsedTime(&ms, s->e0, s->e1) == cudaSuccess) s->tKernel += ms * 1e-3;
    s->nPairsDone += nPairs;
    return true;
}

void erTimings(const Er* s, double* tk, long long* np)
{
    if (tk) *tk = s ? s->tKernel : 0.0;
    if (np) *np = s ? s->nPairsDone : 0;
}

const char* erLastError() { return lastErrEr.c_str(); }

bool erMathCheck(int device, const double* x, long long n, double* e, double* l)
{
    CKE(cudaSetDevice(device));
    double *dx = nullptr, *de = nullptr, *dl = nullptr;
    const std::size_t b = (std::size_t)n * sizeof(double);
    if (cudaMalloc(&dx, b) != cudaSuccess || cudaMalloc(&de, b) != cudaSuccess || cudaMalloc(&dl, b) != cudaSuccess) {
        cudaFree(dx); cudaFree(de); cudaFree(dl);
        lastErrEr = "erMathCheck: cudaMalloc";
        return false;
    }
    bool ok = cudaMemcpy(dx, x, b, cudaMemcpyHostToDevice) == cudaSuccess;
    if (ok) { er_math<<<1024, 256>>>(dx, n, de, dl); ok = cudaGetLastError() == cudaSuccess; }
    if (ok) ok = cudaMemcpy(e, de, b, cudaMemcpyDeviceToHost) == cudaSuccess &&
                 cudaMemcpy(l, dl, b, cudaMemcpyDeviceToHost) == cudaSuccess;
    cudaFree(dx); cudaFree(de); cudaFree(dl);
    if (!ok) lastErrEr = "erMathCheck: CUDA error";
    return ok;
}

}  // namespace gpu2
}  // namespace saige
