// gpu_er_fp32.cu — gpuPrecisionER: fp32. The same exact test as gpu_er.cu
// (same 2^k case/control assignments of the k carriers, same strata and
// lexicographic order, same two kernels: one thread per pair for k <= 12, one
// block per pair above), evaluated in fp32 in the log domain. Not bit-identical
// to the CPU; the fp64 path (gpu_er.cu) is.
//
// What changes against the fp64 kernels, and why each step is safe in fp32:
//
//  * Fisher weights. The CPU multiplies up to 20 odds (and divides for the
//    inverted strata); in fp32 such products underflow (odds 0.01, 20 cases:
//    1e-40). Here each assignment carries lr(S) = sum_{i in S} lo_i with
//    lo_i = log(mu_i) - log(1 - mu_i) centred on their mean (the centring is a
//    constant per stratum, so it cancels in the within-stratum normalisation
//    and keeps |lr| small, i.e. its absolute rounding small).
//  * Normalisation. Per stratum j, three online log-sum-exp accumulators with
//    Kahan-compensated sums: D_j over all assignments (the CPU's denomi[j]),
//    N_j over T >= T_obs, S_j over the ties. The three CPU passes collapse into
//    one: p = sum_j prob_j N_j / D_j / total, total = sum_j prob_j.
//  * Null distribution prob[j]. The hypergeometric leaves are accumulated as
//    log-sum-exp per j and normalised in the log domain.
//  * The O(N) loops of the fp64 kernel are gone. The non-carriers' mean of mu
//    is (sum of the trait's mu - carriers' mu) / (n - k) with the trait's sum
//    computed once at erCreate in double and kept as a float pair (hi, lo).
//    log C(n - k, ncase - j) enters only through differences between strata, so
//    it is built from log C(N, c - 1) - log C(N, c) = log c - log(N - c + 1)
//    over k terms instead of the CPU's ncase-term running sum (whose absolute
//    value, ~3e4 at N = 5e4, has an fp32 ulp of 4e-3). The CPU's convention
//    lCombinations(n, k > n) = 0 is reproduced on the branch where it can be
//    reached (ncase > n - k), with a Kahan sum.
//  * The statistic and its ties. T(S) = (G_S - c)^2 with G_S = sum_{i in S} g_i
//    and c = sum_i g_i mu_i; the CPU's Q - T is evaluated as
//    (G_obs - G_S)(G_obs + G_S - 2c), with G_S folded as a float-float pair, so
//    an exact tie gives exactly 0 and the |Q - T| <= 1e-6 test sees an error of
//    ~1e-7 relative instead of fp32's 1e-4 absolute at T ~ 1e3.
//  * Bin membership (10 bins of mu) uses the same double comparisons as the
//    fp64 kernel, so the groups are the CPU's.
//  * The result: p = exp(L_num - L_tot) - exp(L_same - L_tot) / 2 with the two
//    exps in double, so a p below fp32's range (1e-38) is still returned.
//
// fp64 arithmetic left in the pair prologue (O(k) per pair, none in the 2^k
// enumeration): 1 - mu at the carriers, and G_all, G_obs, 2c, each split into a
// float pair; plus the final two exps.

#include "gpu_er_fp32.cuh"

#include <cmath>

namespace saige {
namespace gpu2 {
namespace erfp32 {

namespace {

constexpr int KM = kErMaxCarriers;
constexpr int NBIN = 10;
constexpr int KSMALL = 12;          // == gpu_er.cu
constexpr int BT = 256;             // threads of the large-k kernel

// ---- compensated log-sum-exp accumulator: value = exp(m) * (s - c) ----
struct Acc { float m, s, c; };
__device__ __forceinline__ Acc accInit() { return Acc{-INFINITY, 0.f, 0.f}; }
__device__ __forceinline__ void accAdd(Acc& A, float x)
{
    if (x > A.m) {   // A.m = -inf: the factor is 0 and s, c are 0
        const float f = expf(A.m - x);
        A.s *= f; A.c *= f; A.m = x;
    }
    const float y = expf(x - A.m) - A.c;
    const float t = A.s + y;
    A.c = (t - A.s) - y;
    A.s = t;
}
__device__ __forceinline__ float accLog(const Acc& A)
{
    const float v = A.s - A.c;
    return (v > 0.f) ? A.m + logf(v) : -INFINITY;
}

// ---- float-float ----
struct FF { float h, l; };
__device__ __forceinline__ FF ffAdd(FF a, float bh, float bl)
{
    const float s = a.h + bh;
    const float bb = s - a.h;
    float e = (a.h - (s - bb)) + (bh - bb);
    e += a.l + bl;
    const float h = s + e;
    return FF{h, e - (h - s)};
}
__device__ __forceinline__ FF ffSplit(double x)
{
    const float h = (float)x;
    return FF{h, (float)(x - (double)h)};
}

// Everything the enumeration needs; the hand-off record of large-k pairs.
struct EnumIn {
    int k;
    float lo[KM];            // centred log odds
    float gh[KM], gl[KM];    // dosage as a float pair
    float sumLo;             // sum of lo (start of the inverted strata)
    FF gAll, gObs, c2;       // sum g, sum_{cases} g, 2 sum g mu
    float lprob[KM + 1];     // log prob[j]
};

// cnt consecutive r-subsets of [0, k) in lexicographic order from a[], into
// the stratum accumulators (D: all, N: T >= T_obs, S: ties). The prefix folds
// of lr and G are kept and only the changed suffix refolded, as in gpu_er.cu.
__device__ __forceinline__ void walk(const EnumIn& E, int j, int* a, long long cnt, Acc& D, Acc& N, Acc& S)
{
    const int k = E.k;
    const bool inv = !(j <= k / 2 + 1);
    const int r = inv ? (k - j) : j;
    float lr[KM + 1];
    FF G[KM + 1];
    lr[0] = inv ? E.sumLo : 0.f;
    G[0] = inv ? E.gAll : FF{0.f, 0.f};
    int from = 0;
    for (long long u = 0; u < cnt; u++) {
        for (int q = from; q < r; q++) {
            const int l = a[q];
            if (!inv) { lr[q + 1] = lr[q] + E.lo[l]; G[q + 1] = ffAdd(G[q],  E.gh[l],  E.gl[l]); }
            else      { lr[q + 1] = lr[q] - E.lo[l]; G[q + 1] = ffAdd(G[q], -E.gh[l], -E.gl[l]); }
        }
        const float x = lr[r];
        accAdd(D, x);
        const FF a1 = ffAdd(E.gObs, -G[r].h, -G[r].l);              // G_obs - G_S
        const FF a2 = ffAdd(ffAdd(E.gObs, G[r].h, G[r].l), -E.c2.h, -E.c2.l);   // G_obs + G_S - 2c
        float temp1 = (a1.h + a1.l) * (a2.h + a2.l);                 // Q - T
        if (fabsf(temp1) <= 1e-6f) temp1 = 0.f;   // 1e-6f is the largest float <= 1e-6
        if (temp1 <= 0.f) {
            accAdd(N, x);
            if (temp1 == 0.f) accAdd(S, x);
        }
        int q = r - 1;
        while (q >= 0 && a[q] == k - r + q) q--;
        if (q < 0) break;
        a[q]++;
        for (int v = q + 1; v < r; v++) a[v] = a[v - 1] + 1;
        from = q;
    }
}

// One stratum's three sums into the pair's p-value accumulators.
__device__ __forceinline__ void foldStratum(float lprob, float lD, float lN, float lS, Acc& Pn, Acc& Ps)
{
    if (lN > -INFINITY) accAdd(Pn, lprob + lN - lD);
    if (lS > -INFINITY) accAdd(Ps, lprob + lS - lD);
}

__device__ __forceinline__ double finalP(const Acc& Pn, const Acc& Ps, const Acc& T)
{
    const float lt = accLog(T);
    return exp((double)(accLog(Pn) - lt)) - 0.5 * exp((double)(accLog(Ps) - lt));
}

__global__ void __launch_bounds__(32)
er_pairs32(const double* __restrict__ mu, const long long* __restrict__ muOff,
           const int* __restrict__ nOf, const int* __restrict__ ncaseOf, const float* __restrict__ muSum,
           const ErPairIn* __restrict__ in, int nPairs,
           const uint32_t* __restrict__ carIdx, const double* __restrict__ carG,
           const unsigned char* __restrict__ carCase,
           ErPairOut* __restrict__ out, const int* __restrict__ bigSlot, EnumIn* __restrict__ prep)
{
    for (int pq = blockIdx.x * blockDim.x + threadIdx.x; pq < nPairs; pq += gridDim.x * blockDim.x) {
        const ErPairIn P = in[pq];
        const int t = P.trait, k = P.k;
        const int n = nOf[t], ncase = ncaseOf[t];
        const double* m = mu + muOff[t];
        EnumIn E;
        E.k = k;
        double pd[KM];
        bool cs[KM];
        float csum = 0.f, loMean = 0.f;
        double gAll = 0.0, gObs = 0.0, gmu = 0.0;
        for (int i = 0; i < k; i++) {
            pd[i] = m[carIdx[P.off + i]];
            const double g = carG[P.off + i];
            cs[i] = carCase[P.off + i] != 0;
            const float pf = (float)pd[i];
            csum += pf;
            E.lo[i] = logf(pf) - logf((float)(1.0 - pd[i]));
            loMean += E.lo[i];
            const FF gs = ffSplit(g);
            E.gh[i] = gs.h; E.gl[i] = gs.l;
            gAll += g; if (cs[i]) gObs += g; gmu += g * pd[i];
        }
        loMean /= (float)k;
        E.sumLo = 0.f;
        for (int i = 0; i < k; i++) { E.lo[i] -= loMean; E.sumLo += E.lo[i]; }
        E.gAll = ffSplit(gAll); E.gObs = ffSplit(gObs); E.c2 = ffSplit(2.0 * gmu);

        // ---- groups (SKATExactBin_ComputeProb_Group) ----
        const float p2 = ((muSum[2 * t] - csum) + muSum[2 * t + 1]) / (float)(n - k);
        const float lodd2 = logf(p2) - log1pf(-p2);
        float lw[NBIN + 1];
        int group[NBIN + 1];
        int ng = 0;
        for (int b = 0; b < NBIN; b++) {
            const double a1 = double(b) / NBIN;
            const double a2 = double(b + 1) / NBIN;
            float acc = 0.f;
            int cnt = 0;
            for (int i = 0; i < k; i++) {
                const double pc = (pd[i] >= 1) ? 0.999 : pd[i];
                const bool inBin = (pc >= a1) && ((b + 1) < NBIN ? (pc < a2) : (pc <= a2));
                if (!inBin) continue;
                acc += (float)pc;
                cnt++;
            }
            if (cnt > 0) {
                const float p1t = acc / (float)cnt;
                lw[ng] = (logf(p1t) - log1pf(-p1t)) - lodd2;
                group[ng] = cnt;
                ng++;
            }
        }
        lw[ng] = 0.f;   // the non-carriers' weight is p2odd / p2odd = 1
        group[ng] = n - k;
        const int ngroup = ng + 1;

        // ---- HyperGeo tables, log scale ----
        float tabC[KM + NBIN];
        int tabOff[NBIN];
        {
            int o = 0;
            for (int i = 0; i < ngroup - 1; i++) {
                tabOff[i] = o;
                float lch = 0.f;
                for (int j = 0; j <= group[i]; j++) {
                    tabC[o + j] = lch + (float)j * lw[i];
                    lch += logf((float)(group[i] - j)) - logf((float)(j + 1));
                }
                o += group[i] + 1;
            }
        }
        // last group: log C(nn0, ncase - j) up to a constant (lw = 0 there).
        float tabL[KM + 1];
        {
            const int nn0 = group[ngroup - 1];
            if (ncase <= nn0) {
                float rel = 0.f;
                for (int j = 0; j <= k; j++) {
                    const int c = ncase - j;
                    tabL[j] = (c >= 0) ? rel : 0.f;   // c < 0 is never a leaf
                    if (c >= 1) rel += logf((float)c) - logf((float)(nn0 - c + 1));
                }
            } else {
                // the CPU's lCombinations(nn0, c) with its 0 for c > nn0, absolute
                for (int j = 0; j <= k; j++) tabL[j] = 0.f;
                float r = 0.f, cmp = 0.f;
                for (int d = 1; d <= nn0; ++d) {
                    const float y = (logf((float)(nn0 - d + 1)) - logf((float)d)) - cmp;
                    const float s2 = r + y;
                    cmp = (s2 - r) - y;
                    r = s2;
                    const int j = ncase - d;
                    if (j >= 0 && j <= k) tabL[j] = r - cmp;
                }
            }
        }
        // Recursive(0, 0, 0), depth first; leaves into log-sum-exp per j
        Acc lk[KM + 1];
        for (int j = 0; j <= k; j++) lk[j] = accInit();
        {
            const int G = ngroup - 1;
            int sel[NBIN + 1];
            float ps[NBIN + 1];
            int nc[NBIN + 1];
            int lev = 0;
            ps[0] = 0.f; nc[0] = 0; sel[0] = -1;
            while (lev >= 0) {
                if (lev == G) {
                    accAdd(lk[nc[G]], ps[G] + tabL[nc[G]]);
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
        {
            Acc all = accInit();
            float lv[KM + 1];
            for (int j = 0; j <= k; j++) { lv[j] = accLog(lk[j]); if (lv[j] > -INFINITY) accAdd(all, lv[j]); }
            const float la = accLog(all);
            for (int j = 0; j <= k; j++) E.lprob[j] = lv[j] - la;
        }

        if (k > KSMALL) { prep[bigSlot[pq]] = E; continue; }

        // ---- the 2^k assignments, one pass ----
        Acc Pn = accInit(), Ps = accInit(), T = accInit();
        int a[KM];
        for (int j = 0; j <= k; j++) {
            if (!(E.lprob[j] > -INFINITY)) continue;   // prob[j] = 0: no contribution
            accAdd(T, E.lprob[j]);
            const bool inv = !(j <= k / 2 + 1);
            const int r = inv ? (k - j) : j;
            for (int q = 0; q < r; q++) a[q] = q;
            Acc D = accInit(), N = accInit(), S = accInit();
            walk(E, j, a, 1LL << 40, D, N, S);
            foldStratum(E.lprob[j], accLog(D), accLog(N), accLog(S), Pn, Ps);
        }
        out[pq].pval = finalP(Pn, Ps, T);
    }
}

__device__ __forceinline__ void reduceAcc(float* rm, float* rv, int tid)
{
    for (int s = BT / 2; s > 0; s >>= 1) {
        if (tid < s) {
            const float m1 = rm[tid], m2 = rm[tid + s];
            const float mx = fmaxf(m1, m2);
            if (mx == -INFINITY) { rm[tid] = -INFINITY; rv[tid] = 0.f; }
            else { rv[tid] = rv[tid] * expf(m1 - mx) + rv[tid + s] * expf(m2 - mx); rm[tid] = mx; }
        }
        __syncthreads();
    }
}

__global__ void __launch_bounds__(BT)
er_big32(const EnumIn* __restrict__ prep, const int* __restrict__ bigPair, int nBig, ErPairOut* __restrict__ out)
{
    const int b = blockIdx.x;
    if (b >= nBig) return;
    const int tid = threadIdx.x;
    __shared__ EnumIn E;
    __shared__ int Cb[KM + 1][KM + 1];
    __shared__ float rm[3][BT], rv[3][BT];
    if (tid == 0) {
        E = prep[b];
        for (int nn = 0; nn <= KM; nn++)
            for (int r = 0; r <= KM; r++)
                Cb[nn][r] = (r == 0) ? 1 : (nn == 0 ? 0 : Cb[nn - 1][r - 1] + Cb[nn - 1][r]);
    }
    __syncthreads();
    const int k = E.k;
    Acc Pn = accInit(), Ps = accInit(), T = accInit();   // thread 0's
    for (int j = 0; j <= k; j++) {
        if (!(E.lprob[j] > -INFINITY)) continue;   // uniform across the block
        const bool inv = !(j <= k / 2 + 1);
        const int r = inv ? (k - j) : j;
        const long long total = Cb[k][j];
        const long long per = (total + BT - 1) / BT;
        const long long c0 = (long long)tid * per;
        const long long cnt = (c0 < total) ? ((c0 + per < total) ? per : total - c0) : 0;
        Acc D = accInit(), N = accInit(), S = accInit();
        if (cnt > 0) {
            // unrank: the c0-th r-subset of [0, k) in lexicographic order
            int a[KM];
            long long mm = c0;
            int x = 0;
            for (int q = 0; q < r; q++) {
                while (true) {
                    const int c = Cb[k - x - 1][r - q - 1];
                    if (c <= mm) { mm -= c; x++; } else break;
                }
                a[q] = x++;
            }
            walk(E, j, a, cnt, D, N, S);
        }
        const Acc* A3[3] = {&D, &N, &S};
        for (int w = 0; w < 3; w++) { rm[w][tid] = A3[w]->m; rv[w][tid] = A3[w]->s - A3[w]->c; }
        __syncthreads();
        reduceAcc(rm[0], rv[0], tid);
        reduceAcc(rm[1], rv[1], tid);
        reduceAcc(rm[2], rv[2], tid);
        if (tid == 0) {
            accAdd(T, E.lprob[j]);
            auto lg = [](float m, float v) { return (v > 0.f) ? m + logf(v) : -INFINITY; };
            foldStratum(E.lprob[j], lg(rm[0][0], rv[0][0]), lg(rm[1][0], rv[1][0]), lg(rm[2][0], rv[2][0]), Pn, Ps);
        }
        __syncthreads();
    }
    if (tid == 0) out[bigPair[b]].pval = finalP(Pn, Ps, T);
}

}  // namespace

std::size_t prepBytes() { return sizeof(EnumIn); }

cudaError_t launch(cudaStream_t st,
                   const double* mu, const long long* muOff, const int* nOf, const int* ncaseOf,
                   const float* muSum,
                   const ErPairIn* in, int nPairs,
                   const uint32_t* carIdx, const double* carG, const unsigned char* carCase,
                   ErPairOut* out, const int* bigSlot, const int* bigPair, int nBig, void* prep)
{
    const int nt = 32;
    const int nb = (nPairs + nt - 1) / nt;
    er_pairs32<<<nb, nt, 0, st>>>(mu, muOff, nOf, ncaseOf, muSum, in, nPairs, carIdx, carG, carCase,
                                  out, bigSlot, static_cast<EnumIn*>(prep));
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess || nBig <= 0) return e;
    er_big32<<<nBig, BT, 0, st>>>(static_cast<const EnumIn*>(prep), bigPair, nBig, out);
    return cudaGetLastError();
}

}  // namespace erfp32
}  // namespace gpu2
}  // namespace saige
