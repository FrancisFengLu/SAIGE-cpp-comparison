// gpu_firth_fp32.cuh — the fp32 variant of firth_pairs (gpuPrecisionFirth: fp32).
// Included by gpu_firth.cu inside its anonymous namespace, after the fp64
// kernel and its helpers (NT, NACC, blockSum, blockSum9, inv2x2).
//
// Same control flow as the fp64 kernel, statement for statement: pass A
// (b = XV g over carriers), pass B (g~ = g - XXVX_inv b), then
// fast_logistf_fit_simple's Newton loop with the step cap, maxit (which counts
// as converged, strict = 0), the singular-information exit and the stopping
// rule "max|delta| <= xconv and max|U*| <= gconv" with U* from before the
// step, with the caller's xconv and gconv unchanged.
//
// What is fp32 and what is not:
//   fp32  -- everything per sample and per iteration: y, offset, XV and
//            XXVX_inv are stored as float (converted once in firthCreate);
//            the dosage, g~ (scratch), the linear predictors, the logistic,
//            w, the residual, the Firth weight, and each thread's running
//            sums: FIRTH32_CHUNK samples added plainly, each chunk then added
//            Neumaier-compensated (pair s + c).
//   fp64  -- the cross-thread reductions (each thread's s + c widened to
//            double, then the fp64 kernel's fixed-shape block tree), the 2x2
//            inverse, U*, the step, alpha / beta / cov(beta) and the stopping
//            test. beta, alpha, se are handed back as double. The caller's
//            xconv and gconv are used unchanged.
//
// Why the residual is split around a reference point. Written directly,
// y_i - pi_i carries about one fp32 ulp of an O(1) number per sample, and that
// rounding changes from iteration to iteration as alpha / beta move, so the
// computed U* has a floor near eps32 sqrt(N) (2e-5 to 2e-4 at N = 50k here)
// that does not shrink at the root; for a rare variant F^-1 turns it into
// steps above xconv, and the first version of this kernel (tolerance raised
// to that floor) ran 899 of 2,245 bt_full fits to maxit where fp64 ran 7.
// So each iteration evaluates the logistic at a reference point
// (a_ref, b_ref) -- pi_ref, q_ref = 1 - pi_ref, r_ref = y - pi_ref -- and
// the step from it, dd = (alpha - a_ref) + (beta - b_ref) g~:
//   y - pi = r_ref - dpi,  dpi = pi - pi_ref = pi_ref q_ref (e^dd - 1) / (q_ref + pi_ref e^dd)
// with expm1f, so dpi is relatively accurate when small. While the fit is
// moving (last max|delta| > FIRTH32_FREEZE) the reference is the current
// iterate, rounded to float; once the steps are below that it stays put.
// Then pi_ref and r_ref are recomputed from identical inputs in an identical
// order every iteration, so the sums of r_ref and g~ r_ref (accumulated on
// their own) are bit-identical from one iteration to the next: their
// rounding is a fixed offset of the score, which moves the root by
// F^-1 x (that offset), not noise around it. Only dpi changes, and it is
// small, so the computed U* goes down to well below gconv like fp64's.
//
// Overflow-safe logistic, for x = eta_ref and for the step dd alike: with
// e = exp(-|x|) in (0, 1], pi = 1 / (1 + e), q = e / (1 + e) (swapped for
// x < 0); for the step, with E = exp(-|dd|) = expm1f(-|dd|) + 1,
//   dd >= 0: pi = pi_r / (q_r E + pi_r),  q = q_r E / (q_r E + pi_r),  dpi = -pi_r q_r expm1(-dd) / (q_r E + pi_r)
//   dd <  0: pi = pi_r E / (q_r + pi_r E), q = q_r / (q_r + pi_r E),   dpi =  pi_r q_r expm1(dd)  / (q_r + pi_r E)
// so nothing overflows and both pi and q = 1 - pi are relatively accurate;
// w = pi q, 0.5 - pi = (q - pi) / 2.
// Termination is the fp64 rule's: maxit (50) ends every fit, reported as fp64
// reports it (conv = 1, strict = 0); a non-finite sum fails inv2x2 and ends
// the fit as singular (conv = 0).
#pragma once

// reference point frozen once the last Newton step is below this (see above)
#define FIRTH32_FREEZE 1e-3
// per-thread sums: this many samples added plainly, then compensated
#define FIRTH32_CHUNK 4

__device__ __forceinline__ void nadd(float& s, float& c, float x)
{
    const float t = s + x;
    if (fabsf(s) >= fabsf(x)) c += (s - t) + x;
    else                      c += (x - t) + s;
    s = t;
}

__device__ __forceinline__ float dose32(const unsigned char* col, const float4& L, int i)
{
    const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
    return (c == 0) ? L.x : (c == 1 ? L.y : (c == 2 ? L.z : L.w));
}

template <bool OWN>
__global__ void __launch_bounds__(NT)
firth_pairs_f32(const unsigned char* __restrict__ packed, std::size_t bpv, const double* __restrict__ lut, int N,
                const double* __restrict__ plut, const uint64_t* __restrict__ masks, int maskWords,
                const float* __restrict__ Y, const float* __restrict__ OFF,
                const float* __restrict__ XV, const float* __restrict__ XX,
                const int* __restrict__ pOfTrait, std::size_t traitStride,
                const FirthPairIn* __restrict__ in, int nPairs, float* __restrict__ scratch,
                int maxit, double maxstep, double xconv, double gconv,
                FirthPairOut* __restrict__ out)
{
    __shared__ double sh[NACC * (NT / 32) + NACC];
    __shared__ float bsh[FIRTH_PMAX];
    float* gt = scratch + (std::size_t)blockIdx.x * N;
    for (int k = blockIdx.x; k < nPairs; k += gridDim.x) {
        const FirthPairIn pin = in[k];
        const unsigned char* col = packed + (std::size_t)pin.slot * bpv;
        float4 L;
        {
            const double* d = (OWN ? plut : lut) + 4 * (std::size_t)(OWN ? k : pin.slot);
            L = make_float4((float)d[0], (float)d[1], (float)d[2], (float)d[3]);
        }
        const int t = pin.trait;
        const uint64_t* mk = OWN ? masks + (std::size_t)t * maskWords : nullptr;
        auto in = [&](int i) -> bool { return !OWN || ((mk[i >> 6] >> (i & 63)) & 1ull); };
        const int p = pOfTrait[t];
        const float* y   = Y   + (std::size_t)t * N;
        const float* off = OFF + (std::size_t)t * N;
        const float* xv  = XV + (std::size_t)t * traitStride;   // sample i: xv[i*p + j]
        const float* xx  = XX + (std::size_t)t * traitStride;   // xx[j*N + i]

        // pass A: b = XV g over carriers, FIRTH_PTILE columns per pass
        for (int j0 = 0; j0 < p; j0 += FIRTH_PTILE) {
            const int pt = (p - j0 < FIRTH_PTILE) ? p - j0 : FIRTH_PTILE;
            float as[FIRTH_PTILE], ac[FIRTH_PTILE];
            for (int j = 0; j < FIRTH_PTILE; ++j) { as[j] = 0.f; ac[j] = 0.f; }
            for (int i = threadIdx.x; i < N; i += NT) {
                const float g = in(i) ? dose32(col, L, i) : 0.f;
                if (g == 0.f) continue;
                const float* x = xv + (std::size_t)i * p + j0;
                for (int j = 0; j < pt; ++j) nadd(as[j], ac[j], x[j] * g);
            }
            for (int j = 0; j < pt; ++j) {
                const double s = blockSum((double)as[j] + (double)ac[j], sh);
                if (threadIdx.x == 0) bsh[j0 + j] = (float)s;
            }
        }
        __syncthreads();
        // pass B: g~ = g - XXVX_inv b
        for (int i = threadIdx.x; i < N; i += NT) {
            const float g = in(i) ? dose32(col, L, i) : 0.f;
            float proj = 0.f;
            for (int j = 0; j < p; ++j) proj += xx[(std::size_t)j * N + i] * bsh[j];
            gt[i] = g - proj;
        }
        __syncthreads();

        // ---- fast_logistf_fit_simple(x = [1, g~], y, offset, firth, init = 0) ----
        double alpha = 0.0, beta = 0.0, cov11 = 0.0;
        double aRef = 0.0, bRef = 0.0, lastStep = HUGE_VAL;
        int iter = 0, flag = 0, strict = 0, singular = 0;
        while (iter <= maxit) {
            if (lastStep > FIRTH32_FREEZE) { aRef = (double)(float)alpha; bRef = (double)(float)beta; }
            const float ar = (float)aRef, br = (float)bRef;
            const float da = (float)(alpha - aRef), db = (float)(beta - bRef);
            float s[NACC + 2], c[NACC + 2];
            for (int q = 0; q < NACC + 2; q++) { s[q] = 0.f; c[q] = 0.f; }
            // FIRTH32_CHUNK samples summed plainly, then the chunk added compensated
            for (int i0 = threadIdx.x; i0 < N; i0 += NT * FIRTH32_CHUNK) {
                float u[NACC + 2];
                for (int q = 0; q < NACC + 2; q++) u[q] = 0.f;
                #pragma unroll
                for (int h = 0; h < FIRTH32_CHUNK; ++h) {
                    const int i = i0 + h * NT;
                    if (i >= N || !in(i)) continue;
                    const float g = gt[i], yi = y[i];
                    // the reference point
                    const float x = (ar + br * g) + off[i];
                    const float ex = expf(-fabsf(x)), rx = 1.f / (1.f + ex);
                    const float pr = (x >= 0.f) ? rx : ex * rx;
                    const float qr = (x >= 0.f) ? ex * rx : rx;
                    const float rr = yi * qr - (1.f - yi) * pr;            // y - pi_ref
                    // the step from it
                    const float dd = da + db * g;
                    const float em = expm1f(-fabsf(dd)), E = 1.f + em;
                    float pi, qq, dp;
                    if (dd >= 0.f) {
                        const float den = 1.f / (qr * E + pr);
                        pi = pr * den; qq = (qr * E) * den; dp = -((pr * qr) * em) * den;
                    } else {
                        const float den = 1.f / (qr + pr * E);
                        pi = (pr * E) * den; qq = qr * den; dp = ((pr * qr) * em) * den;
                    }
                    const float w = pi * qq;
                    const float ad = w * (0.5f * (qq - pi));               // w (0.5 - pi)
                    const float wg = w * g, adg = ad * g, adg2 = adg * g;
                    u[0] += w;   u[1] += wg;  u[2] += wg * g;
                    u[3] += dp;  u[4] += g * dp;
                    u[5] += ad;  u[6] += adg; u[7] += adg2; u[8] += adg2 * g;
                    u[9] += rr;  u[10] += g * rr;
                }
                for (int q = 0; q < NACC + 2; q++) nadd(s[q], c[q], u[q]);
            }
            double acc[NACC], ref[2];
            for (int q = 0; q < NACC; q++) acc[q] = (double)s[q] + (double)c[q];
            const double rs0 = blockSum((double)s[9] + (double)c[9], sh);
            const double rs1 = blockSum((double)s[10] + (double)c[10], sh);
            ref[0] = rs0; ref[1] = rs1;
            blockSum9(acc, sh);
            double A, Bc, C;
            if (!inv2x2(acc[0], acc[1], acc[2], A, Bc, C)) { singular = 1; break; }
            const double R0s = ref[0] - acc[3], R1s = ref[1] - acc[4];   // X' (y - pi)
            const double U0 = R0s + (A * acc[5] + 2.0 * Bc * acc[6] + C * acc[7]);
            const double U1 = R1s + (A * acc[6] + 2.0 * Bc * acc[7] + C * acc[8]);
            double d0 = A * U0 + Bc * U1, d1 = Bc * U0 + C * U1;
            double mxd = fmax(fabs(d0), fabs(d1));
            const double mx = mxd / maxstep;
            if (mx > 1.0) { d0 /= mx; d1 /= mx; mxd = fmax(fabs(d0), fabs(d1)); }
            cov11 = C;
            iter++;
            alpha += d0; beta += d1;
            lastStep = mxd;
            const bool small = (mxd <= xconv) && (fabs(U0) <= gconv) && (fabs(U1) <= gconv);
            if (iter == maxit || small) { flag = 1; strict = small ? 1 : 0; break; }
        }
        if (threadIdx.x == 0) {
            FirthPairOut o;
            o.beta = beta; o.alpha = alpha; o.se = sqrt(cov11);
            o.conv = flag; o.strict = strict; o.niter = iter; o.singular = singular;
            out[k] = o;
        }
        __syncthreads();
    }
}
