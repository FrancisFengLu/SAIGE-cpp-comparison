#!/usr/bin/env python3
"""
Region-test memory budget planner for SAIGE-GENE+ C++.

Model (validated 2026-07-29 against local N=1000 runs and the two UKB
measurements in SAIGE_GENE_STATUS.md §4):

    peak_GB(N, T, C) = base_GB + T * (2 * C * N * 8 / 2**30)

where C is the effective chunk row-count used to size P1Mat/P2Mat.

IMPORTANT: before the max_markers_region fix, C was NOT
markers_per_chunk_in_groupTest — main.cpp resized P1Mat/P2Mat to
max_markers_region (default 100000) regardless of the chunk knob. Use
--legacy to model that behavior. See BENCH_REPORT.md.

Validation:
    N=1000,   C=100000, T=1  -> predicted 1.60 GB, observed 1.54 GB
    N=165582, C=100000, T=1  -> predicted 246.7 GB, observed ~242 GB (T2D)
    N=158598, C=100000, T=1  -> predicted 236.3 GB, observed ~230 GB (LDL)

Usage:
    plan_budget.py --N 165000 --budget-gb 200
    plan_budget.py --N 165000 --budget-gb 200 --legacy   # pre-fix behavior
"""
import argparse, sys

BASE_GB = 0.05          # observed floor at N=1000 post-fix; grows slowly with N
LEGACY_CHUNK = 100000   # max_markers_region default, used pre-fix

def per_thread_gb(N, chunk):
    """P1Mat (chunk x N) + P2Mat (N x chunk), doubles."""
    return 2 * chunk * N * 8 / (1024 ** 3)

def peak_gb(N, T, chunk):
    return BASE_GB + T * per_thread_gb(N, chunk)

def plan(N, budget_gb, tmin, tmax, chunks):
    cand = []
    for T in range(tmin, tmax + 1):
        for C in chunks:
            gb = peak_gb(N, T, C)
            if gb <= budget_gb:
                # prefer more threads, then larger chunk (fewer spills -> stays
                # on the single-chunk, R-validated code path)
                cand.append((-T * 1_000_000 - C, T, C, gb))
    cand.sort()
    return [(T, C, gb) for _, T, C, gb in cand[:10]]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True, help="sample size")
    ap.add_argument("--budget-gb", type=float, required=True)
    ap.add_argument("--min-threads", type=int, default=1)
    ap.add_argument("--max-threads", type=int, default=32)
    ap.add_argument("--max-gene-markers", type=int, default=2000,
                    help="largest expected passing-marker count in any gene; "
                         "chunks below this force nchunks>1, which is NOT yet "
                         "validated against R (see BENCH_REPORT.md §4)")
    ap.add_argument("--legacy", action="store_true",
                    help="model pre-fix behavior (chunk pinned to 100000)")
    args = ap.parse_args()

    if args.legacy:
        gb = peak_gb(args.N, 1, LEGACY_CHUNK)
        print(f"# LEGACY (pre-fix): P1/P2 pinned to max_markers_region={LEGACY_CHUNK}")
        print(f"# per-thread = {per_thread_gb(args.N, LEGACY_CHUNK):.1f} GB")
        print(f"# T=1 peak   = {gb:.1f} GB   (budget {args.budget_gb} GB)")
        if gb > args.budget_gb:
            print("# -> WILL OOM even single-threaded. Apply the fix.")
        sys.exit(0)

    chunks = [c for c in (500, 1000, 2000, 3000, 5000, 10000, 20000)
              if c >= args.max_gene_markers]
    if not chunks:
        chunks = [args.max_gene_markers]

    print(f"# N={args.N}, budget={args.budget_gb} GB")
    print(f"# only chunk >= {args.max_gene_markers} considered "
          f"(keeps nchunks==1, the R-validated path)")
    print(f"# per-thread @ chunk={chunks[0]}: "
          f"{per_thread_gb(args.N, chunks[0]):.2f} GB")
    print()

    top = plan(args.N, args.budget_gb, args.min_threads, args.max_threads, chunks)
    if not top:
        need = peak_gb(args.N, 1, chunks[0])
        print(f"!! nothing fits: even T=1 chunk={chunks[0]} needs {need:.1f} GB",
              file=sys.stderr)
        sys.exit(1)

    print(f"{'nThreads':>8}  {'markers_per_chunk':>18}  {'est_peak_GB':>11}")
    for T, C, gb in top:
        print(f"{T:>8}  {C:>18}  {gb:>11.2f}")
    print()
    T, C, gb = top[0]
    print(f"# suggested: nThreads={T}, markers_per_chunk_in_groupTest={C}")
    print(f"# est peak ~ {gb:.1f} GB (< {args.budget_gb} GB budget)")

if __name__ == "__main__":
    main()
