#!/usr/bin/env python3
"""Firth-corrected logistic regression for a batch of (marker, trait) pairs on the
GPU -- ONE CUDA BLOCK PER PAIR, fp64, the same parallel shape as the step-2 SPA
prototype (optimization/gpuassess/bench_spa_gpu.py).

What runs per block is the iteration of SAIGEClass::fast_logistf_fit_simple
(saige_test.cpp) for X = [1, gtilde], written in its two-parameter moment form:
every quantity of one Newton step -- the 2x2 Fisher information, the hat-matrix
diagonal h_i = w_i x_i' F^-1 x_i and the Firth-adjusted score -- is a sum over the
N samples of a function of (gtilde_i, y_i, offset_i, alpha, beta), so one pass
over N with nine accumulators replaces the N x 2 weighted matrix, the QR and the
hat vector. The step cap (max|delta| <= maxstep, maxstep = 15 in SAIGE, 1 in the
collaborator package), the stopping rule (max|delta| <= xconv and max|U*| <= gconv,
tested with the score from BEFORE the step), the maxit = 50 exit that SAIGE also
reports as "converged", and the singular-information exit are all kept verbatim,
so each block iterates exactly as many times as the CPU would for that pair.

gtilde = G - XXVX_inv * (XV * G) is formed on the device from the marker's 2-bit
.bed column through the per-marker code->dosage table (flip / mean imputation /
MAC-gated zeroing folded in, as in gpu_step2.cu) and the per-pair p-vector
b = XV * G; the per-trait N x p matrix XXVX_inv, y and offset are resident.
So the per-pair host->device payload is 2 ints + p doubles (32 B at p = 3) and
the per-marker payload is the N/4-byte packed column, shared by all traits.

Usage:
  firth_gpu.py validate <bed> <N> <pairdir> [--maxstep 15|1] [--blocks 320] [--rep R]
  firth_gpu.py synth <N> <npairs> [--maxstep 15]            # throughput only
"""
import sys, os, time, math, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arma_io import load_arma
from tsv import read_tsv

SRC = r'''
#define NT 256
#define NACC 9
__device__ __forceinline__ double warp_sum(double v) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffff, v, o);
    return v;
}
// Sum NACC values across the block; every thread gets the totals in acc[].
// Fixed shape (warp shuffle tree, then thread 0 adds the 8 warp partials in
// order) -> bit-reproducible run to run.
__device__ void block_sum9(double* acc, double* sh) {
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int q = 0; q < NACC; q++) {
        double v = warp_sum(acc[q]);
        if (lane == 0) sh[q * (NT / 32) + wid] = v;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        for (int q = 0; q < NACC; q++) {
            double s = 0.0;
            for (int w = 0; w < NT / 32; w++) s += sh[q * (NT / 32) + w];
            sh[NACC * (NT / 32) + q] = s;
        }
    }
    __syncthreads();
    for (int q = 0; q < NACC; q++) acc[q] = sh[NACC * (NT / 32) + q];
    __syncthreads();
}

extern "C" __global__ void __launch_bounds__(NT)
firth_pairs(const unsigned char* __restrict__ bedcols, long long Bbytes, int N,
            const double* __restrict__ lut,            // M x 4 code->dosage
            const double* __restrict__ XXVX,           // T x N x p (row-major)
            int p,
            const double* __restrict__ Y,              // T x N
            const double* __restrict__ OFF,            // T x N
            const int* __restrict__ pm, const int* __restrict__ pt,
            const double* __restrict__ pb,             // npairs x p  (b = XV * G)
            double* __restrict__ scratch,              // gridDim.x x N
            int npairs, int maxit, double maxstep, double xconv, double gconv,
            double* __restrict__ out,                  // npairs x 4: alpha, beta, se, cov11
            int* __restrict__ iout)                    // npairs x 4: flag, strict, iters, singular
{
    __shared__ double sh[NACC * (NT / 32) + NACC];
    double* gt = scratch + (size_t)blockIdx.x * N;
    const double eps = 2.220446049250313e-16;
    for (int k = blockIdx.x; k < npairs; k += gridDim.x) {
        const int m = pm[k], t = pt[k];
        const unsigned char* col = bedcols + (size_t)m * Bbytes;
        const double* L = lut + (size_t)m * 4;
        const double* Xt = XXVX + (size_t)t * N * p;
        const double* y = Y + (size_t)t * N;
        const double* off = OFF + (size_t)t * N;
        const double* b = pb + (size_t)k * p;
        // ---- gtilde = G - XXVX_inv * b  (getadjGFast) ----
        for (int i = threadIdx.x; i < N; i += NT) {
            const int code = (col[i >> 2] >> (2 * (i & 3))) & 3;
            double v = L[code];
            double xb = 0.0;
            for (int j = 0; j < p; j++) xb += Xt[(size_t)i * p + j] * b[j];
            gt[i] = v - xb;
        }
        __syncthreads();
        // ---- fast_logistf_fit_simple, X = [1, gtilde], init = 0 ----
        double alpha = 0.0, beta = 0.0, cov11 = 0.0;   // XX_covs starts as zeros
        int iter = 0, flag = 0, strict = 0, singular = 0;
        while (iter <= maxit) {
            double acc[NACC];
            for (int q = 0; q < NACC; q++) acc[q] = 0.0;
            for (int i = threadIdx.x; i < N; i += NT) {
                const double g = gt[i];
                // legacy: pi = 1 / (exp(-x*beta - offset) + 1)
                const double pi = 1.0 / (exp(-(alpha + beta * g) - off[i]) + 1.0);
                const double w = pi * (1.0 - pi);
                const double W2 = sqrt(w), c = g * W2;        // XW2 columns
                const double r = y[i] - pi;
                const double a = w * (0.5 - pi);               // Firth term weight
                acc[0] += W2 * W2;  acc[1] += W2 * c;  acc[2] += c * c;   // Fisher = XW2' XW2
                acc[3] += r;        acc[4] += g * r;                      // X' (y - pi)
                acc[5] += a;        acc[6] += a * g;  acc[7] += a * g * g;  acc[8] += a * g * g * g;
            }
            block_sum9(acc, sh);
            const double s0 = acc[0], s1 = acc[1], s2 = acc[2];
            const double det = s0 * s2 - s1 * s1;
            // inv_sympd 2x2 (armadillo apply_tiny_2x2): a<=0, det<eps, det>1/eps or NaN
            // fail the closed form; det>1/eps then succeeds in LAPACK, the others do not.
            if (!(s0 > 0.0) || !(det >= eps)) { singular = 1; break; }
            const double A = s2 / det, Bc = -s1 / det, C = s0 / det;   // F^-1
            // U* = X' ((y - pi) + h (0.5 - pi)),  h_i = w_i (A + 2 Bc g_i + C g_i^2)
            const double U0 = acc[3] + (A * acc[5] + 2.0 * Bc * acc[6] + C * acc[7]);
            const double U1 = acc[4] + (A * acc[6] + 2.0 * Bc * acc[7] + C * acc[8]);
            double d0 = A * U0 + Bc * U1, d1 = Bc * U0 + C * U1;      // delta = XX_covs * U*
            double mxd = fmax(fabs(d0), fabs(d1));
            const double mx = mxd / maxstep;
            if (mx > 1.0) { d0 /= mx; d1 /= mx; mxd = fmax(fabs(d0), fabs(d1)); }
            cov11 = C;
            iter++;
            alpha += d0; beta += d1;
            const bool small = (mxd <= xconv) && (fabs(U0) <= gconv) && (fabs(U1) <= gconv);
            if (iter == maxit || small) { flag = 1; strict = small ? 1 : 0; break; }
        }
        if (threadIdx.x == 0) {
            out[4 * k] = alpha; out[4 * k + 1] = beta; out[4 * k + 2] = sqrt(cov11); out[4 * k + 3] = cov11;
            iout[4 * k] = flag; iout[4 * k + 1] = strict; iout[4 * k + 2] = iter; iout[4 * k + 3] = singular;
        }
        __syncthreads();
    }
}
'''

def load_set(bed, N, pairdir, rep=1):
    """Host-side inputs for one pair directory. Returns dict of numpy arrays."""
    B = (N + 3) // 4
    raw = np.fromfile(bed, dtype=np.uint8)
    assert raw[:3].tolist() == [0x6c, 0x1b, 0x01]
    lut = np.load(os.path.join(pairdir, 'lut.npy'))
    M = lut.shape[0]; assert raw.shape[0] == 3 + M * B
    P = read_tsv(os.path.join(pairdir, 'pairs.tsv'))
    pm = P['marker'].astype(np.int32); pt = P['trait'].astype(np.int32)
    traits = [l.rstrip('\n').split('\t') for l in open(os.path.join(pairdir, 'traits.tsv'))]
    T = len(traits); p = None
    Y = []; OFF = []; XXVX = []; XV = []
    for ti, name, mdir in traits:
        Y.append(load_arma(mdir + '/y.arma')); OFF.append(load_arma(mdir + '/offset.arma'))
        xx = load_arma(mdir + '/XXVX_inv.arma'); XXVX.append(np.ascontiguousarray(xx)); XV.append(load_arma(mdir + '/XV.arma'))
        p = xx.shape[1]
    # markers actually used: decode once, b = XV * G per pair (sparse in the CPU; dense here, host side)
    cols = raw[3:].reshape(M, B)
    idx = np.arange(N); byte = idx >> 2; shift = (2 * (idx & 3)).astype(np.uint8)
    pb = np.empty((len(pm), p))
    cache = {}
    for k in range(len(pm)):
        m = int(pm[k]); t = int(pt[k])
        if (m) not in cache:
            code = (cols[m][byte] >> shift) & 3
            cache[m] = lut[m][code]
        G = cache[m]
        # same order as getadjGFast: m_XVG += XV.col(i) * G(i) over the carriers, in index order
        nz = np.nonzero(G)[0]
        pb[k] = np.cumsum(XV[t][:, nz] * G[nz], axis=1)[:, -1] if len(nz) else 0.0
    if rep > 1:
        pm = np.tile(pm, rep); pt = np.tile(pt, rep); pb = np.tile(pb, (rep, 1))
    return dict(bedcols=raw[3:], B=B, N=N, lut=lut, XXVX=np.stack(XXVX), p=p, Y=np.stack(Y), OFF=np.stack(OFF),
                pm=pm, pt=pt, pb=pb, pairs=P, ntraits=T)

def run_kernel(d, maxstep=15.0, maxit=50, xconv=1e-5, gconv=1e-5, nblocks=320, timing_runs=3):
    import cupy as cp
    kern = cp.RawKernel(SRC, 'firth_pairs')
    npairs = len(d['pm']); N = d['N']
    dev = dict(bedcols=cp.asarray(d['bedcols']), lut=cp.asarray(d['lut']), XXVX=cp.asarray(d['XXVX']),
               Y=cp.asarray(d['Y']), OFF=cp.asarray(d['OFF']), pm=cp.asarray(d['pm']), pt=cp.asarray(d['pt']), pb=cp.asarray(d['pb']))
    nb = min(nblocks, npairs)
    scratch = cp.empty(nb * N, np.float64)
    out = cp.empty(npairs * 4, np.float64); iout = cp.empty(npairs * 4, np.int32)
    args = (dev['bedcols'], np.int64(d['B']), np.int32(N), dev['lut'], dev['XXVX'], np.int32(d['p']), dev['Y'], dev['OFF'],
            dev['pm'], dev['pt'], dev['pb'], scratch, np.int32(npairs), np.int32(maxit), np.float64(maxstep),
            np.float64(xconv), np.float64(gconv), out, iout)
    kern((nb,), (256,), args); cp.cuda.Device().synchronize()
    ts = []
    for _ in range(timing_runs):
        cp.cuda.Device().synchronize(); t0 = time.perf_counter(); kern((nb,), (256,), args); cp.cuda.Device().synchronize()
        ts.append(time.perf_counter() - t0)
    o = out.get().reshape(-1, 4); io = iout.get().reshape(-1, 4)
    return dict(alpha=o[:, 0], beta=o[:, 1], se=o[:, 2], flag=io[:, 0], strict=io[:, 1], iters=io[:, 2], singular=io[:, 3],
                t=min(ts), ts=ts, nblocks=nb)

def qval(p):
    """|Phi^-1(p/2)| as saige_test.cpp uses it for the Firth SE back-calculation (boost quantile)."""
    from math import sqrt
    try:
        from scipy.special import ndtri
        return np.abs(ndtri(p / 2.0))
    except Exception:
        import statistics
        nd = statistics.NormalDist()
        return np.array([abs(nd.inv_cdf(v / 2.0)) if np.isfinite(v) and 0 < v < 1 else np.nan for v in p])

def validate(a):
    d = load_set(a.bed, a.N, a.pairdir, rep=a.rep)
    cpu = read_tsv(os.path.join(a.pairdir, a.cpu))
    ms = a.maxstep
    r = run_kernel(d, maxstep=ms, nblocks=a.blocks)
    npairs0 = len(cpu['pair']); n = len(d['pm'])
    print(f'[{a.pairdir}] N={a.N} pairs={n} (distinct {npairs0}, rep {a.rep}) traits={d["ntraits"]} maxstep={ms} blocks={r["nblocks"]}')
    print(f'  GPU kernel: {r["t"]*1e3:.1f} ms -> {r["t"]/n*1e6:.2f} us/pair  (runs: {", ".join(f"{t*1e3:.1f}" for t in r["ts"])} ms)')
    passes = (r['iters'] + 1).sum()
    print(f'  iterations: mean {r["iters"].mean():.2f} max {r["iters"].max()}  -> {passes*a.N/r["t"]/1e9:.1f} G sample-passes/s (decode pass + one per iteration)')
    # compare the first npairs0 (the CPU file) against the matching tag
    tag = 'leg' if ms == 15 else ('m1' if ms == 1 else None)
    if tag is None:
        print('  (no CPU column for this maxstep)'); return
    cb = cpu['beta_' + tag][:npairs0]; gb = r['beta'][:npairs0]
    cflag = cpu['flag_' + tag][:npairs0]; gflag = r['flag'][:npairs0]
    citer = cpu['iter_' + tag][:npairs0]; giter = r['iters'][:npairs0]
    if tag == 'leg':
        cstrict = cpu['strict_leg'][:npairs0]; csing = cpu['sing_leg'][:npairs0]; cse = cpu['se_leg'][:npairs0]
    else:
        cstrict = cpu['conv_m1'][:npairs0]; csing = (cpu['fail_m1'][:npairs0] == 2).astype(float); cse = cpu['se_m1'][:npairs0]
    gse = r['se'][:npairs0]
    ok = np.isfinite(cb) & np.isfinite(gb)
    absd = np.abs(gb - cb); reld = absd / np.maximum(np.abs(cb), 1e-300)
    conv_both = ok & (cstrict == 1) & (r['strict'][:npairs0] == 1)
    print(f'  flag identical: {int((cflag == gflag).sum())}/{npairs0}   strict identical: {int((cstrict == r["strict"][:npairs0]).sum())}/{npairs0}   '
          f'singular identical: {int((csing == r["singular"][:npairs0]).sum())}/{npairs0}   iteration count identical: {int((citer == giter).sum())}/{npairs0}')
    print(f'  CPU: strict {int(cstrict.sum())}, hit maxit {int(((cflag==1)&(cstrict==0)).sum())}, singular {int(csing.sum())};  '
          f'GPU: strict {int(r["strict"][:npairs0].sum())}, hit maxit {int(((gflag==1)&(r["strict"][:npairs0]==0)).sum())}, singular {int(r["singular"][:npairs0].sum())}')
    if conv_both.any():
        print(f'  beta, strictly converged on both ({int(conv_both.sum())}): max |d| {absd[conv_both].max():.2e}  max rel {reld[conv_both].max():.2e}  '
              f'median rel {np.median(reld[conv_both]):.1e}')
        sd = np.abs(gse - cse)[conv_both]
        print(f'  se (sqrt cov11, the fit SE before SAIGE overwrites it): max |d| {sd.max():.2e}  max rel {(sd/np.abs(cse[conv_both])).max():.2e}')
    nc = ok & ~conv_both
    if nc.any():
        print(f'  beta, not strictly converged ({int(nc.sum())}): max |d| {absd[nc].max():.2e}  (these are the oscillating fits; see FIRTH_GPU.md)')
    # SAIGE's reported SE = |beta| / qnorm(p/2); p is untouched by Firth
    P = d['pairs']; pp = P['p'][:npairs0]
    if np.isfinite(pp).any():
        q = qval(pp); seC = np.abs(cb) / q; seG = np.abs(gb) / q
        m = np.isfinite(seC) & conv_both
        print(f'  reported SE = |beta|/qnorm(p/2) (p unchanged by Firth): max |d| {np.abs(seG-seC)[m].max():.2e}  max rel {(np.abs(seG-seC)/seC)[m].max():.2e}')
        if 'beta_out' in P and np.isfinite(P['beta_out'][:npairs0]).any():
            # the fit sees the flipped genotype; main.cpp writes Beta * (1 - 2*flip)
            sgn = 1.0 - 2.0 * P['flip'][:npairs0]
            bo = P['beta_out'][:npairs0]; mo = np.isfinite(bo) & (bo != 0)
            print(f'  vs production BETA (6 printed digits, flip sign applied): max rel {(np.abs(gb*sgn-bo)/np.abs(bo))[mo].max():.2e} over {int(mo.sum())} pairs')
    if a.dump:
        with open(a.dump, 'w') as f:
            f.write('pair\tbeta_gpu\tse_gpu\tflag_gpu\tstrict_gpu\titer_gpu\tsing_gpu\talpha_gpu\n')
            for k in range(npairs0):
                f.write(f'{int(cpu["pair"][k])}\t{float(gb[k])!r}\t{float(gse[k])!r}\t{gflag[k]}\t{r["strict"][k]}\t{giter[k]}\t{r["singular"][k]}\t{float(r["alpha"][k])!r}\n')

def synth(a):
    """Throughput at a biobank N: random 2-bit columns (no missing), y ~ 5% cases, offset ~ N(-3, 0.7)."""
    rng = np.random.default_rng(3); N = a.N; B = (N + 3) // 4; Mk = a.markers; p = 3; T = 4
    maf = np.exp(rng.uniform(math.log(0.001), math.log(0.5), Mk))
    G = rng.binomial(2, maf[:, None], size=(Mk, N)).astype(np.uint8)
    code = np.array([3, 2, 0], np.uint8)[G]
    cols = np.zeros((Mk, B), np.uint8)
    for j in range(4): cols[:, : (N - j + 3) // 4] |= (code[:, j::4] << (2 * j)).astype(np.uint8)
    lut = np.tile(np.array([2.0, 0.0, 1.0, 0.0]), (Mk, 1))
    X = np.column_stack([np.ones(N), rng.normal(size=N), rng.normal(size=N)])
    Y = np.empty((T, N)); OFF = np.empty((T, N)); XXVX = np.empty((T, N, p)); XV = np.empty((T, p, N))
    for t in range(T):
        off = rng.normal(0, 0.7, N) + X @ rng.normal(0, 0.3, p); mu = 1 / (1 + np.exp(-(off + math.log(a.case_rate / (1 - a.case_rate)))))
        Y[t] = rng.random(N) < mu; OFF[t] = off; w = mu * (1 - mu)
        XtWX = (X.T * w) @ X; XV[t] = np.linalg.solve(XtWX, X.T * w); XXVX[t] = X
    pm = rng.integers(0, Mk, a.npairs).astype(np.int32); pt = rng.integers(0, T, a.npairs).astype(np.int32)
    pb = np.empty((a.npairs, p))
    for k in range(a.npairs): pb[k] = XV[pt[k]] @ G[pm[k]]
    d = dict(bedcols=cols.reshape(-1), B=B, N=N, lut=lut, XXVX=XXVX, p=p, Y=Y, OFF=OFF, pm=pm, pt=pt, pb=pb, ntraits=T)
    r = run_kernel(d, maxstep=a.maxstep, nblocks=a.blocks)
    n = a.npairs
    print(f'[synth] N={N} pairs={n} markers={Mk} case_rate={a.case_rate} maxstep={a.maxstep} blocks={r["nblocks"]}')
    print(f'  GPU kernel: {r["t"]*1e3:.1f} ms -> {r["t"]/n*1e6:.2f} us/pair  (runs: {", ".join(f"{t*1e3:.1f}" for t in r["ts"])} ms)')
    print(f'  iterations: mean {r["iters"].mean():.2f} p99 {np.percentile(r["iters"],99):.0f} max {r["iters"].max()}  strict {int(r["strict"].sum())}/{n}  hit maxit {int(((r["flag"]==1)&(r["strict"]==0)).sum())}  singular {int(r["singular"].sum())}'
          f'  -> {(r["iters"]+1).sum()*N/r["t"]/1e9:.1f} G sample-passes/s')
    print(f'  device-resident per trait: {2*N*8/1e6:.1f} MB (y, offset) + {N*p*8/1e6:.1f} MB (XXVX_inv); per marker {B/1e3:.1f} KB packed; per pair {8+8*p} B')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd', required=True)
    v = sub.add_parser('validate'); v.add_argument('bed'); v.add_argument('N', type=int); v.add_argument('pairdir')
    v.add_argument('--cpu', default='cpu_ref.tsv'); v.add_argument('--maxstep', type=float, default=15.0)
    v.add_argument('--blocks', type=int, default=320); v.add_argument('--rep', type=int, default=1); v.add_argument('--dump', default=None)
    s = sub.add_parser('synth'); s.add_argument('N', type=int); s.add_argument('npairs', type=int); s.add_argument('--markers', type=int, default=2000)
    s.add_argument('--case_rate', type=float, default=0.05); s.add_argument('--maxstep', type=float, default=15.0); s.add_argument('--blocks', type=int, default=320)
    a = ap.parse_args()
    import cupy as cp
    print(cp.cuda.runtime.getDeviceProperties(0)['name'].decode())
    validate(a) if a.cmd == 'validate' else synth(a)
