// GPU 锁步 SPA（binary，SAIGE_MT_GPUSPA=1 时启用）
//
// 每个 block 处理一个回落对的一条尾（q 或 qinv）：块内归约算 K1/K2/Korg，
// 块内 Newton（init=0，tol、maxiter=1000、prevJump 二分保护）→ 鞍点概率。
// 公式与 SPA_binary.cpp 的 Korg/K1_adj/K2_fast_Binom 逐项一致（e^{-gt} 变形）。
// 非 fast 情形（p_com<0.5）等价于 NAmu=NAsigma=0 + 数组取全长 —— 同一 kernel。
//
// 只处理 logp=false 的对（pval_noadj==0 的极端对留在 CPU 线程安全路径）。
// pnorm 用 erfc 实现，与 Rmath 有末位差 —— 所以本路径不承诺逐字节一致，
// 验收看差异统计。
#include <cstdio>
#include <cmath>
#include <cuda_runtime.h>

#define SPAG_NTHR 256
#define SPAG_MAXIT 1000

// 上尾标准正态：P(Z > z)。erfc 版，双精度。
__device__ __forceinline__ double pnorm_upper(double z)
{
    return 0.5 * erfc(z / 1.4142135623730951);
}

struct SpaJob {          // 一条尾 = 一个 job
    int    off;          // 数组起点（扁平缓冲）
    int    len;          // 元素数
    double q;            // 该尾的 q（或 qinv）
    double namu, nasig;  // fast 形式的正态补偿；非 fast 传 0
    double gpos, gneg;   // Σg⁺ / Σg⁻（getroot 的早退判据，CPU 端算好）
};

__global__ void spa_gpu_kernel(const double* __restrict__ G,
                               const double* __restrict__ MU,
                               const SpaJob* __restrict__ jobs,
                               double* __restrict__ pout,     // 该尾的鞍点 p（带符号，同 CPU 语义）
                               unsigned char* __restrict__ okroot,
                               unsigned char* __restrict__ oksaddle,
                               double tol, int njob)
{
    __shared__ double s1[SPAG_NTHR], s2[SPAG_NTHR];
    __shared__ double t_sh, k1_sh, k2_sh;
    __shared__ int    state_sh;   // 0 迭代中 / 1 收敛 / 2 失败

    int jb = blockIdx.x;
    if (jb >= njob) return;
    SpaJob J = jobs[jb];
    const double *g  = G  + J.off;
    const double *mu = MU + J.off;
    int n = J.len;

    // getroot 早退：q 超出 K' 的可达范围
    if (J.q >= J.gpos || J.q <= J.gneg) {
        if (threadIdx.x == 0) { okroot[jb] = 1; oksaddle[jb] = 0; pout[jb] = 0; }
        return;   // root=inf：CPU 端 Get_Saddle_Prob(inf) 的 k1/k2 非有限 → isSaddle=false
    }

    // ---- Newton：K1(t) = Σ μg/((1−μ)e^{−gt}+μ) + NAmu + NAsigma·t − q ----
    auto K12 = [&](double t, double &K1, double &K2){
        double a1 = 0, a2 = 0;
        for (int i = threadIdx.x; i < n; i += SPAG_NTHR) {
            double gi = g[i], mi = mu[i];
            double e  = exp(-gi * t);
            double d  = (1.0 - mi) * e + mi;
            a1 += mi * gi / d;
            a2 += (1.0 - mi) * mi * gi * gi * e / (d * d);
        }
        s1[threadIdx.x] = a1; s2[threadIdx.x] = a2;
        __syncthreads();
        for (int w = SPAG_NTHR/2; w > 0; w >>= 1) {
            if (threadIdx.x < w) { s1[threadIdx.x] += s1[threadIdx.x+w]; s2[threadIdx.x] += s2[threadIdx.x+w]; }
            __syncthreads();
        }
        K1 = s1[0] + J.namu + J.nasig * t - J.q;
        K2 = s2[0] + J.nasig;
        __syncthreads();
    };

    if (threadIdx.x == 0) { t_sh = 0; state_sh = 0; }
    __syncthreads();

    double K1_eval, K2_eval;
    K12(0.0, K1_eval, K2_eval);
    double prevJump = INFINITY;
    double t = 0;

    for (int rep = 1; rep <= SPAG_MAXIT; ++rep) {
        double tnew = t - K1_eval / K2_eval;
        if (isnan(tnew)) { if (threadIdx.x == 0) state_sh = 2; __syncthreads(); break; }
        if (fabs(tnew - t) < tol) { if (threadIdx.x == 0) { t_sh = t; state_sh = 1; } __syncthreads(); break; }
        if (rep == SPAG_MAXIT) { if (threadIdx.x == 0) state_sh = 2; __syncthreads(); break; }
        double newK1, newK2;
        K12(tnew, newK1, newK2);
        if ((K1_eval * newK1) < 0) {
            if (fabs(tnew - t) > (prevJump - tol)) {
                double sgn = (newK1 - K1_eval) > 0 ? 1.0 : ((newK1 - K1_eval) < 0 ? -1.0 : 0.0);
                tnew = t + sgn * prevJump / 2;
                K12(tnew, newK1, newK2);
                prevJump = prevJump / 2;
            } else {
                prevJump = fabs(tnew - t);
            }
        }
        t = tnew;  K1_eval = newK1;  K2_eval = newK2;
    }
    __syncthreads();

    if (state_sh == 2) {
        if (threadIdx.x == 0) { okroot[jb] = 0; oksaddle[jb] = 0; pout[jb] = 0; }
        return;
    }
    t = t_sh;

    // ---- Get_Saddle_Prob：k1 = Korg(t)，k2 = K2(t) ----
    double a0 = 0;
    for (int i = threadIdx.x; i < n; i += SPAG_NTHR) {
        double gi = g[i], mi = mu[i];
        a0 += log(1.0 - mi + mi * exp(gi * t));
    }
    s1[threadIdx.x] = a0;
    __syncthreads();
    for (int w = SPAG_NTHR/2; w > 0; w >>= 1) {
        if (threadIdx.x < w) s1[threadIdx.x] += s1[threadIdx.x+w];
        __syncthreads();
    }
    if (threadIdx.x == 0) k1_sh = s1[0] + J.namu * t + 0.5 * J.nasig * t * t;
    __syncthreads();
    double dummyK1, k2v;
    K12(t, dummyK1, k2v);          // K2(t)；K1 值不用
    if (threadIdx.x == 0) k2_sh = k2v;
    __syncthreads();

    if (threadIdx.x == 0) {
        double k1 = k1_sh, k2 = k2_sh;
        double temp1 = t * J.q - k1;
        double pval = 0;  unsigned char isSaddle = 0;
        bool flagrun = false;
        double w = 0, v = 0;
        if (isfinite(k1) && isfinite(k2) && temp1 >= 0 && k2 >= 0) {
            double sgn = (t > 0) ? 1.0 : ((t < 0) ? -1.0 : 0.0);
            w = sgn * sqrt(2 * temp1);
            v = t * sqrt(k2);
            if (w != 0) flagrun = true;
        }
        if (flagrun) {
            double Ztest = w + (1 / w) * log(v / w);
            if (Ztest > 0) pval = pnorm_upper(Ztest);
            else           pval = -pnorm_upper(-Ztest);   // 下尾 = P(Z<z)=Φ̄(−z)，符号语义同 CPU
            isSaddle = 1;
        }
        okroot[jb] = 1;  oksaddle[jb] = isSaddle;  pout[jb] = pval;
    }
}

// ---- host 侧接口（Main.cpp 调）----
extern "C" {

int spa_gpu_available()
{
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n > 0 ? 1 : 0;
}

// jobs 打包好的一批尾；返回 0 成功。G/MU 是扁平缓冲（所有 job 共用）。
int spa_gpu_solve(int njob, long total_elems,
                  const double *G, const double *MU, const void *jobs_v,
                  double tol,
                  double *pout, unsigned char *okroot, unsigned char *oksaddle)
{
    const SpaJob *jobs = (const SpaJob*)jobs_v;
    static double *dG = nullptr, *dMU = nullptr, *dp = nullptr;
    static SpaJob *dJ = nullptr;
    static unsigned char *dr = nullptr, *ds = nullptr;
    static long   capE = 0;  static int capJ = 0;

    if (total_elems > capE) {
        if (dG)  cudaFree(dG);   if (dMU) cudaFree(dMU);
        if (cudaMalloc(&dG,  total_elems * 8) != cudaSuccess) return 1;
        if (cudaMalloc(&dMU, total_elems * 8) != cudaSuccess) return 1;
        capE = total_elems;
    }
    if (njob > capJ) {
        if (dJ) cudaFree(dJ);  if (dp) cudaFree(dp);
        if (dr) cudaFree(dr);  if (ds) cudaFree(ds);
        if (cudaMalloc(&dJ, njob * sizeof(SpaJob)) != cudaSuccess) return 1;
        if (cudaMalloc(&dp, njob * 8) != cudaSuccess) return 1;
        if (cudaMalloc(&dr, njob) != cudaSuccess) return 1;
        if (cudaMalloc(&ds, njob) != cudaSuccess) return 1;
        capJ = njob;
    }
    if (cudaMemcpy(dG,  G,  total_elems * 8, cudaMemcpyHostToDevice) != cudaSuccess) return 2;
    if (cudaMemcpy(dMU, MU, total_elems * 8, cudaMemcpyHostToDevice) != cudaSuccess) return 2;
    if (cudaMemcpy(dJ, jobs, njob * sizeof(SpaJob), cudaMemcpyHostToDevice) != cudaSuccess) return 2;

    spa_gpu_kernel<<<njob, SPAG_NTHR>>>(dG, dMU, dJ, dp, dr, ds, tol, njob);
    if (cudaDeviceSynchronize() != cudaSuccess) return 3;

    if (cudaMemcpy(pout, dp, njob * 8, cudaMemcpyDeviceToHost) != cudaSuccess) return 4;
    if (cudaMemcpy(okroot,  dr, njob, cudaMemcpyDeviceToHost) != cudaSuccess) return 4;
    if (cudaMemcpy(oksaddle, ds, njob, cudaMemcpyDeviceToHost) != cudaSuccess) return 4;
    return 0;
}

} // extern "C"
