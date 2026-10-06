// gpu_er_glibc.cuh — device ports of the exp and log the CPU build calls, so the
// device ER (gpu_er.cu) takes the same transcendental values bit for bit.
//
// The CPU's exp / log are glibc 2.35's (sysdeps/ieee754/dbl-64/e_exp.c,
// e_log.c). On x86-64 they are ifuncs; on a CPU with FMA + AVX2 -- every machine
// this runs on -- the resolver picks __exp_fma / __log_fma, the same C source
// built with -mfma -mavx2, where gcc contracted some a*b+c into fused
// multiply-adds. Which ones is not visible in the source, so the statements
// below follow the disassembly of those two functions in the installed
// libm.so.6 one instruction at a time (S2_ER_GPU.md §2): every vfmadd is an
// explicit fma() here and every other operation is a separately rounded
// __dmul_rn / __dadd_rn / __dsub_rn, which nvcc never contracts. The constants
// and tables are read out of that libm (gpu_er_glibc_tables.cuh).
//
// Covered: every path of both functions (the |x| < 2^-54, overflow and
// underflow paths of exp including the subnormal-result rounding, the
// near-1 polynomial and the subnormal-input path of log). errno and the
// floating-point exception flags are not reproduced -- nothing reads them.
// Checked against the host functions on random and edge inputs by
// tools/er_gpu_check (0 differing bits required).
#pragma once

#include <cstdint>
#include "gpu_er_glibc_tables.cuh"

namespace saige {
namespace gpu2 {
namespace glx {

__device__ __forceinline__ double asd(uint64_t u) { return __longlong_as_double((long long)u); }
__device__ __forceinline__ uint64_t asu(double d) { return (uint64_t)__double_as_longlong(d); }

// __exp_fma, glibc 2.35. N = 128.
__device__ double exp_(double x)
{
    const uint64_t ux = asu(x);
    uint32_t abstop = (uint32_t)(ux >> 52) & 0x7ffu;
    if (abstop - 0x3c9u >= 0x408u - 0x3c9u) {          // top12(0x1p-54), top12(512)
        if ((int32_t)(abstop - 0x3c9u) < 0)
            return __dadd_rn(x, 1.0);                   // |x| < 2^-54
        if (abstop >= 0x409u) {                          // |x| >= 1024
            if (ux == 0xfff0000000000000ULL) return 0.0;
            if (abstop >= 0x7ffu) return __dadd_rn(x, 1.0);
            if (ux >> 63) return 0.0;                    // __math_uflow(0): 0x1p-767 * 0x1p-767
            return asd(0x7ff0000000000000ULL);           // __math_oflow(0)
        }
        abstop = 0;                                      // large |x|: specialcase below
    }
    const double InvLn2N = asd(GLX_EXP_invln2N), Shift = asd(GLX_EXP_shift);
    const double NegLn2hiN = asd(GLX_EXP_negln2hiN), NegLn2loN = asd(GLX_EXP_negln2loN);
    const double C2 = asd(GLX_EXP_C2), C3 = asd(GLX_EXP_C3), C4 = asd(GLX_EXP_C4), C5 = asd(GLX_EXP_C5);

    double kd = fma(x, InvLn2N, Shift);                  // vfmadd132sd: z + Shift fused
    const uint64_t ki = asu(kd);
    kd = __dsub_rn(kd, Shift);
    double r = fma(kd, NegLn2hiN, x);
    r = fma(kd, NegLn2loN, r);
    const uint64_t idx = 2 * (ki % 128);
    const uint64_t top = ki << 45;                       // 52 - EXP_TABLE_BITS
    const double tail = asd(glx_exp_tab[idx]);
    uint64_t sbits = glx_exp_tab[idx + 1] + top;
    const double t23 = fma(r, C3, C2);
    const double rt  = __dadd_rn(r, tail);
    const double r2  = __dmul_rn(r, r);
    const double t45 = fma(r, C5, C4);
    const double a   = fma(t23, r2, rt);
    const double r4  = __dmul_rn(r2, r2);
    const double tmp = fma(r4, t45, a);
    if (abstop == 0) {
        // specialcase(tmp, sbits, ki)
        if ((ki & 0x80000000ULL) == 0) {
            sbits -= 1009ULL << 52;
            const double scale = asd(sbits);
            return __dmul_rn(fma(scale, tmp, scale), asd(GLX_K_p2_1009));
        }
        sbits += 1022ULL << 52;
        const double scale = asd(sbits);
        const double st = __dmul_rn(scale, tmp);
        double y = __dadd_rn(scale, st);
        if (1.0 > y) {
            const double hi = __dadd_rn(y, 1.0);
            double lo = __dsub_rn(scale, y);
            lo = __dadd_rn(lo, st);
            double s = __dsub_rn(1.0, hi);
            s = __dadd_rn(s, y);
            s = __dadd_rn(s, lo);
            s = __dadd_rn(s, hi);
            y = __dsub_rn(s, 1.0);
            if (y == 0.0) return 0.0;
        }
        return __dmul_rn(y, asd(GLX_K_p1022));
    }
    const double scale = asd(sbits);
    return fma(scale, tmp, scale);
}

// __log_fma, glibc 2.35. N = 128, OFF = 0x3fe6000000000000.
__device__ double log_(double x)
{
    uint64_t ix = asu(x);
    const uint32_t top = (uint32_t)(ix >> 48);
    if (ix - 0x3fee000000000000ULL < 0x3ff1090000000000ULL - 0x3fee000000000000ULL) {
        // x in [1 - 2^-4, 1 + 0x1.09p-4)
        if (ix == 0x3ff0000000000000ULL) return 0.0;
        const double B0 = asd(GLX_LOG_B0), B1 = asd(GLX_LOG_B1), B2 = asd(GLX_LOG_B2), B3 = asd(GLX_LOG_B3);
        const double B4 = asd(GLX_LOG_B4), B5 = asd(GLX_LOG_B5), B6 = asd(GLX_LOG_B6), B7 = asd(GLX_LOG_B7);
        const double B8 = asd(GLX_LOG_B8), B9 = asd(GLX_LOG_B9), B10 = asd(GLX_LOG_B10);
        const double two27 = asd(GLX_K_two27);
        const double r    = __dsub_rn(x, 1.0);
        const double b12  = fma(r, B2, B1);
        const double b45  = fma(r, B5, B4);
        const double r2   = __dmul_rn(r, r);
        const double b78  = fma(r, B8, B7);
        const double b123 = fma(r2, B3, b12);
        const double b456 = fma(r2, B6, b45);
        const double r3   = __dmul_rn(r, r2);
        const double b789 = fma(r2, B9, b78);
        const double b710 = fma(r3, B10, b789);
        const double q    = fma(b710, r3, b456);
        const double P    = fma(q, r3, b123);
        const double rw   = fma(r, two27, r);           // r + w, w = r * 2^27
        const double rhi  = fma(-two27, r, rw);         // vfnmadd: rw - 2^27 r
        const double rhi2 = __dmul_rn(rhi, rhi);
        const double rlo  = __dsub_rn(r, rhi);
        const double hi   = fma(rhi2, B0, r);
        const double t    = __dsub_rn(r, hi);
        const double s    = __dadd_rn(r, rhi);
        const double lo   = fma(rhi2, B0, t);
        const double b0rl = __dmul_rn(B0, rlo);
        const double lo2  = fma(b0rl, s, lo);
        const double y    = fma(P, r3, lo2);
        return __dadd_rn(hi, y);
    }
    if (top - 0x0010u >= 0x7ff0u - 0x0010u) {
        if (ix * 2 == 0) return asd(0xfff0000000000000ULL);          // __math_divzero(1)
        if (ix == 0x7ff0000000000000ULL) return x;                    // log(inf)
        if ((top & 0x8000u) || (top & 0x7ff0u) == 0x7ff0u)
            return __ddiv_rn(__dsub_rn(x, x), __dsub_rn(x, x));       // __math_invalid
        ix = asu(__dmul_rn(x, asd(GLX_K_two52)));                     // subnormal: normalise
        ix -= 52ULL << 52;
    }
    const uint64_t tmp = ix - 0x3fe6000000000000ULL;
    const int i = (int)((tmp >> 45) % 128);
    const int k = (int)((int64_t)tmp >> 52);
    const uint64_t iz = ix - (tmp & (0xfffULL << 52));
    const double invc = asd(glx_log_tab[2 * i]);
    const double logc = asd(glx_log_tab[2 * i + 1]);
    const double z = asd(iz);
    const double A0 = asd(GLX_LOG_A0), A1 = asd(GLX_LOG_A1), A2 = asd(GLX_LOG_A2);
    const double A3 = asd(GLX_LOG_A3), A4 = asd(GLX_LOG_A4);
    const double Ln2hi = asd(GLX_LOG_ln2hi), Ln2lo = asd(GLX_LOG_ln2lo);
    const double r   = fma(z, invc, asd(GLX_K_negone));
    const double kd  = (double)k;
    const double w   = fma(kd, Ln2hi, logc);
    const double hi  = __dadd_rn(r, w);
    double t         = __dsub_rn(w, hi);
    t                = __dadd_rn(t, r);
    const double lo  = fma(kd, Ln2lo, t);
    const double r2  = __dmul_rn(r, r);
    const double a12 = fma(r, A2, A1);
    const double a34 = fma(r, A4, A3);
    const double y0  = fma(r2, A0, lo);
    const double p   = fma(a34, r2, a12);
    const double rr2 = __dmul_rn(r, r2);
    const double y   = fma(rr2, p, y0);
    return __dadd_rn(y, hi);
}

}  // namespace glx
}  // namespace gpu2
}  // namespace saige
