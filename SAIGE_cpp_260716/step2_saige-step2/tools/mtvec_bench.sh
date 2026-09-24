#!/bin/bash
# mtvec_bench.sh -- wall clock A/B for mtVecQuantStats, and the emit split.
#
# The covariate fold (mtFoldQuantProj) is ON on both sides: with it off the two
# wide 3P GEMMs cost 3.9 s/trait and bury the tail.
#
# env: DATA (default the 10^6-marker g1m), REPS, PLIST, PROFBIN
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA=${DATA:-/opt/saige/logs/tg2_step2/data/g1m}
REPS=${REPS:-2}
PLIST=${PLIST:-"8 32"}
PROFBIN=${PROFBIN:-/opt/saige/logs/mtvec/bin/saige-step2.prof}
export RUNROOT=${RUNROOT:-/opt/saige/logs/mtvec/runs}
export FOLD=1

wall() { grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K.*' "$1"; }
secs() { awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }'; }

echo "===== wall clock, $DATA, fold ON, text output ====="
for P in $PLIST; do
  for V in 0 1; do
    for R in $(seq 1 "$REPS"); do
      T=bench_P${P}_v${V}_r${R}
      bash "$HERE/mtvec_run_one.sh" "$T" "$DATA" "$V" q:$P > "$RUNROOT/$T.run.log" 2>&1 \
        || { echo "FAILED $T"; tail -5 "$RUNROOT/$T.run.log"; exit 1; }
      W=$(wall "$RUNROOT/$T/time.txt")
      echo "P=$P vec=$V rep=$R  $W  ($(echo "$W" | secs) s)"
      rm -rf "$RUNROOT/$T/out"
    done
  done
done

echo
echo "===== emit split (MTVEC_PROF build, cpu-s) ====="
for P in $PLIST; do
  for V in 0 1; do
    T=prof_P${P}_v${V}
    BIN=$PROFBIN bash "$HERE/mtvec_run_one.sh" "$T" "$DATA" "$V" q:$P \
        > "$RUNROOT/$T.run.log" 2>&1 || { echo "FAILED $T"; exit 1; }
    echo "P=$P vec=$V  $(wall "$RUNROOT/$T/time.txt")"
    grep -E "mtvec prof" "$RUNROOT/$T/log.txt" | sed 's/^/    /'
    rm -rf "$RUNROOT/$T/out"
  done
  T=prof_P${P}_v1_nofmt
  MTVEC_PROF_NOFMT=1 BIN=$PROFBIN bash "$HERE/mtvec_run_one.sh" "$T" "$DATA" 1 q:$P \
      > "$RUNROOT/$T.run.log" 2>&1 || { echo "FAILED $T"; exit 1; }
  echo "P=$P vec=1 NOFMT  $(wall "$RUNROOT/$T/time.txt")"
  grep -E "mtvec prof" "$RUNROOT/$T/log.txt" | sed 's/^/    /'
  rm -rf "$RUNROOT/$T/out"
done
df -h /opt/saige | tail -1
