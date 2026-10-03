// erfc_boost53.cuh -- the normal tail for the device, built from Boost.Math's
// 53-bit erf_imp (boost/math/special_functions/erf.hpp, Boost 1.85, Boost
// Software License 1.0), so that p = erfc(x/sqrt2)/2 agrees with the CPU's
// boost::math::cdf(normal) down to the last p the CPU can represent.
//
// What the CPU actually computes. Boost's default policy promotes double to
// long double (promote_double<true>): erfc(double z) is erf_imp<long double,
// 64-bit tag>(z) -- different rational approximations, 80-bit arithmetic, no
// underflow until z ~ 106 -- then narrowed to double. So the CPU's erfc is,
// to within its ~1e-19 relative error, the CORRECTLY ROUNDED double: a
// normal number for z < 26.55, a correctly rounded subnormal for
// 26.55 <= z < 27.226, and exactly 0 from z ~ 27.2261 (true erfc below half
// the smallest subnormal, 2^-1075). No double-only algorithm can be
// bit-identical to that (erfc_port_hosttest.cpp: the literal 53-bit port
// agrees bit for bit at ~60% of arguments, 1 ulp off elsewhere); what a
// double algorithm CAN reproduce is the correctly rounded subnormal and the
// zero cutoff, which are the two things that change routing (p == 0
// withdraws convergence; log p becomes -inf).
//
// erfImp53s (mode 1, default): Boost's 53-bit rational approximations, with
// the exp(-z^2) factor computed in a scaled form -- exp(-z^2) = 2^-k *
// exp(-z^2 + k ln2), k ln2 split Cody-Waite style so the shifted argument is
// exact -- and the scale applied last by one ldexp, which is the single
// rounding into the subnormal range. Within ~3 ulp of correctly rounded in
// the normal range; in the subnormal range the same value as the CPU except
// when the true value lies within that relative error of a rounding
// boundary; same zero cutoff to within 1e-16 relative in z. The "z >= 28 is
// zero" shortcut is kept (true value < 2^-1130 there).
//
// erfImp53 (mode 2): the literal 53-bit port, including Boost's own cutoff
// (result 0 for z >= 28 even though the true value is a subnormal from
// 26.55 to 27.23) and its unscaled exp(-z^2) that underflows to 0 well
// before 28. Kept for the measurement only.
//
// Both are plain host C++ when not compiled by nvcc (the CUDA qualifiers
// become empty), which is how erfc_port_hosttest.cpp compares them with
// Boost itself on the same libm.
#pragma once
#ifndef __CUDACC__
#  include <cmath>
#  define __device__
#  define __forceinline__ inline
#  define __constant__
#  define SPA_NAN (std::nan(""))
   using std::isnan; using std::floor; using std::ldexp; using std::frexp; using std::exp;
#else
#  include <math_constants.h>
#  define SPA_NAN CUDART_NAN
#endif
#include "erfc_boost53_consts.hpp"

namespace saige { namespace spa_gpu { namespace erfc53 {

// Polynomials as boost/math/tools/detail/polynomial_horner3_20.hpp evaluates
// them for GCC (BOOST_MATH_POLY_METHOD 3): same rounding sequence.
__device__ __forceinline__ double poly5(const double* a, double x)
{
    const double x2 = x * x;
    double t0 = a[4] * x2 + a[2];
    double t1 = a[3] * x2 + a[1];
    t0 *= x2; t0 += a[0];
    t1 *= x;
    return t0 + t1;
}
__device__ __forceinline__ double poly6(const double* a, double x)
{
    const double x2 = x * x;
    double t0 = a[5] * x2 + a[3];
    double t1 = a[4] * x2 + a[2];
    t0 *= x2; t1 *= x2;
    t0 += a[1]; t1 += a[0];
    t0 *= x;
    return t0 + t1;
}
__device__ __forceinline__ double poly7(const double* a, double x)
{
    const double x2 = x * x;
    double t0 = a[6] * x2 + a[4];
    double t1 = a[5] * x2 + a[3];
    t0 *= x2; t1 *= x2;
    t0 += a[2]; t1 += a[1];
    t0 *= x2; t0 += a[0];
    t1 *= x;
    return t0 + t1;
}

// The rational part R(z) = Y + P/Q of Boost's erfc branches for z >= 1.5
// (erfc(z) = R(z) * exp(-z^2) / z).
__device__ __forceinline__ double rationalPart(double z)
{
    if (z < 2.5)      return ERF_Y3 + poly6(ERF_P3, z - 1.5) / poly6(ERF_Q3, z - 1.5);
    else if (z < 4.5) return ERF_Y4 + poly6(ERF_P4, z - 3.5) / poly6(ERF_Q4, z - 3.5);
    else              return ERF_Y5 + poly7(ERF_P5, 1 / z) / poly7(ERF_Q5, 1 / z);
}

// Boost's split of z^2 into the double sq and the exact remainder err_sqr.
__device__ __forceinline__ void splitSquare(double z, double* sq, double* err_sqr)
{
    int expon;
    double hi = floor(ldexp(frexp(z, &expon), 26));
    hi = ldexp(hi, expon - 26);
    const double lo = z - hi;
    *sq = z * z;
    *err_sqr = ((hi * hi - *sq) + 2 * hi * lo) + lo * lo;
}

// --- mode 2: the literal 53-bit erf_imp -------------------------------------
__device__ double erfImp53(double z, bool invert)
{
    if (isnan(z)) return SPA_NAN;        // Boost raises a domain error here
    if (z < 0) {
        if (!invert)       return -erfImp53(-z, false);
        else if (z < -0.5) return 2 - erfImp53(-z, true);
        else               return 1 + erfImp53(-z, false);
    }
    double result;
    if (z < 0.5) {
        if (z < 1e-10) {
            result = (z == 0) ? 0.0 : (z * 1.125 + z * ERF_C_SMALL);
        } else {
            const double zz = z * z;
            result = z * (ERF_Y1 + poly5(ERF_P1, zz) / poly5(ERF_Q1, zz));
        }
    } else if (invert ? (z < 28) : (z < ERF_Z_ERF_LIMIT)) {
        invert = !invert;
        if (z < 1.5) {
            result = ERF_Y2 + poly6(ERF_P2, z - 0.5) / poly7(ERF_Q2, z - 0.5);
            result *= exp(-z * z) / z;
        } else {
            double sq, err_sqr;
            splitSquare(z, &sq, &err_sqr);
            result = rationalPart(z);
            result *= exp(-sq) * exp(-err_sqr) / z;
        }
    } else {
        // Any value of z larger than 28 will underflow to zero (Boost's words).
        result = 0;
        invert = !invert;
    }
    if (invert) result = 1 - result;
    return result;
}

// --- mode 1: the same approximations with the underflow deferred -----------
// ln 2 = LN2_HI + LN2_LO with LN2_HI carrying 32 significant bits, so k*LN2_HI
// is exact for k < 2^21 and (-sq) + k*LN2_HI is exact when the result is
// small (both operands are multiples of the same power of two).
constexpr double LN2_HI = 0.69314718036912382; // 0x1.62e42fep-1
constexpr double LN2_LO = 1.9082149292705877e-10;
constexpr double INV_LN2 = 1.4426950408889634;

__device__ double erfcScaled(double z)   // z >= 1.5
{
    double sq, err_sqr;
    splitSquare(z, &sq, &err_sqr);
    const int k = (int)(sq * INV_LN2) + 1;             // exp(-sq) = 2^-k exp(t) exp(k LN2_LO)
    const double t = (-sq) + (double)k * LN2_HI;       // exact; in (-0.7, 0.7]
    double m = rationalPart(z);
    m *= (exp(t) * exp((double)k * LN2_LO)) * exp(-err_sqr) / z;
    return ldexp(m, -k);                               // the one rounding into the subnormals
}

__device__ double erfImp53s(double z, bool invert)
{
    if (isnan(z)) return SPA_NAN;
    if (z < 0) {
        if (!invert)       return -erfImp53s(-z, false);
        else if (z < -0.5) return 2 - erfImp53s(-z, true);
        else               return 1 + erfImp53s(-z, false);
    }
    double result;
    if (z < 0.5) {
        if (z < 1e-10) {
            result = (z == 0) ? 0.0 : (z * 1.125 + z * ERF_C_SMALL);
        } else {
            const double zz = z * z;
            result = z * (ERF_Y1 + poly5(ERF_P1, zz) / poly5(ERF_Q1, zz));
        }
    } else if (invert ? (z < 28) : (z < ERF_Z_ERF_LIMIT)) {
        invert = !invert;
        if (z < 1.5) {
            result = ERF_Y2 + poly6(ERF_P2, z - 0.5) / poly7(ERF_Q2, z - 0.5);
            result *= exp(-z * z) / z;
        } else {
            result = erfcScaled(z);
        }
    } else {
        result = 0;
        invert = !invert;
    }
    if (invert) result = 1 - result;
    return result;
}

}}}  // namespace saige::spa_gpu::erfc53
