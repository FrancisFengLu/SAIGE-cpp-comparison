#!/bin/bash
# run_gate.sh -- wait until the box is free (no TIMING_LOCK, no step-2 /
# step-1 / benchmark job), then run spa_gpu_test with the given arguments.
# Usage: bash run_gate.sh LOGFILE [spa_gpu_test args...]
set -u
LOG=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
PAT='saige-step2|saige-null|benchmark|spa_gpu_check'
waited=0
while [ -e /opt/saige/logs/TIMING_LOCK ] || pgrep -f "$PAT" > /dev/null; do
    sleep 30; waited=$((waited + 30))
    if [ $((waited % 600)) -eq 0 ]; then echo "$(date +%T) still waiting (${waited}s)" >> "$LOG"; fi
done
echo "$(date +%T) box free after ${waited}s; running: spa_gpu_test $*" >> "$LOG"
"$HERE/spa_gpu_test" "$@" >> "$LOG" 2>&1
echo "exit=$? $(date +%T)" >> "$LOG"
