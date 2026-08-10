#!/usr/bin/env bash
# Driver: generate the benchmark dataset (if absent) and run both benchmark
# phases.  This is the single entry point -- re-run it on a bigger box to
# reproduce/extend the numbers in BENCH_CPP_VS_R.md.
#
#   run_bench_all.sh [data_dir] [out_dir]
#
# Phase A (OPENBLAS_NUM_THREADS=1) is run first and writes timings.csv
# incrementally, so a partial run is still usable.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
DATA="${1:-$REPO/bench_data/N10000_g300}"
OUT="${2:-$HERE/bench_timing_20260810}"
mkdir -p "$OUT"

if [[ ! -f "$DATA/wes.bed" ]]; then
  python3 "$HERE/gen_bench_data.py" --out "$DATA" --n 10000 --grm-markers 10000 \
          --genes 300 --min-var 10 --max-var 100 --seed 21 2>&1 | tee "$OUT/gen_data.log"
fi

# Step 1 is only required single- vs multi-threaded, so it uses {4,1}; the
# step-2 stages get the full {4,2,1} sweep.
STAGES=step1            THREADS_A="4 1"   REPS=3   "$HERE/bench_cpp_vs_r.sh" "$DATA" "$OUT" A
SKIP_PREP=1 STAGES="region single" THREADS_A="4 2 1" REPS=3 "$HERE/bench_cpp_vs_r.sh" "$DATA" "$OUT" A
SKIP_PREP=1 STAGES=step1            THREADS_B="4"     REPS_B=2 "$HERE/bench_cpp_vs_r.sh" "$DATA" "$OUT" B
SKIP_PREP=1 STAGES="region single"  THREADS_B="4 1"   REPS_B=2 "$HERE/bench_cpp_vs_r.sh" "$DATA" "$OUT" B

python3 "$HERE/summarize_bench.py" "$OUT/timings.csv" > "$OUT/tables.md"

# --- equivalence check: same null model in, do the p-values match? -----------
{
  echo "== region  (C++ nThreads=1 rep1  vs  R) =="
  python3 "$HERE/cmp_bench_out.py" region "$OUT/results/A_region_cpp_T1_rep1.out" \
                                          "$OUT/results/A_region_r_T1_rep1.out"
  echo
  echo "== single-variant  (C++ nThreads=1 rep1  vs  R) =="
  python3 "$HERE/cmp_bench_out.py" single "$OUT/results/A_single_cpp_T1_rep1.out" \
                                          "$OUT/results/A_single_r_T1_rep1.out"
  echo
  echo "== C++ nThreads=4 vs nThreads=1 (threading must not change results) =="
  python3 "$HERE/cmp_bench_out.py" region "$OUT/results/A_region_cpp_T4_rep1.out" \
                                          "$OUT/results/A_region_cpp_T1_rep1.out"
} > "$OUT/equivalence.txt" 2>&1

echo "wrote $OUT/tables.md and $OUT/equivalence.txt"
