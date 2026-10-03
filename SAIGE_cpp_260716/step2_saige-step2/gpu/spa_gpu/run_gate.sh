#!/bin/bash
# run_gate.sh -- run spa_gpu_test under the box's sharing protocol, then exit.
#
#   bash run_gate.sh LOGFILE [spa_gpu_test args...]
#   SPA_GPU_TEST_BIN=./spa_gpu_test_integ bash run_gate.sh LOG --impl integ ...
#
# Protocol (coordinator, 2026-10-03). The integrator's timing queue holds
# /opt/saige/logs/TIMING_LOCK for hours. A short job may interleave only when
# that file contains "PAUSE supported":
#   1. create /opt/saige/logs/TIMING_PAUSE with our tag and start time
#      (the queue then starts no new cell);
#   2. wait until no saige-step2 / saige-null / sgs2txt process runs;
#   3. run the job;
#   4. delete TIMING_PAUSE promptly.
# With no lock at all the same process check applies and no pause file is
# written. Any other lock content: keep waiting.
set -u
LOG=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)
LOCK=/opt/saige/logs/TIMING_LOCK
PAUSE=/opt/saige/logs/TIMING_PAUSE
PAT='saige-step2|saige-null|sgs2txt|spa_gpu_check'
TAG="spa-gpu-lib"
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
BIN=${SPA_GPU_TEST_BIN:-$HERE/spa_gpu_test}
wrote_pause=0

cleanup() { if [ "$wrote_pause" = 1 ]; then rm -f "$PAUSE"; echo "$(date +%T) removed TIMING_PAUSE" >> "$LOG"; fi; }
trap cleanup EXIT INT TERM

waited=0
while :; do
    if [ ! -e "$LOCK" ]; then mode=free; break; fi
    if grep -q "PAUSE supported" "$LOCK" 2>/dev/null && [ ! -e "$PAUSE" ]; then mode=pause; break; fi
    sleep 30; waited=$((waited + 30))
    if [ $((waited % 1800)) -eq 0 ]; then echo "$(date +%T) still waiting for the lock protocol (${waited}s)" >> "$LOG"; fi
done
if [ "$mode" = pause ]; then
    echo "$TAG $(date '+%Y-%m-%d %H:%M:%S') short job: $(basename "$BIN") $*" > "$PAUSE"
    wrote_pause=1
    echo "$(date +%T) lock supports PAUSE; wrote TIMING_PAUSE" >> "$LOG"
fi
pw=0
while pgrep -f "$PAT" > /dev/null; do sleep 10; pw=$((pw + 10)); done
echo "$(date +%T) box free (mode=$mode, waited ${waited}s + ${pw}s for processes); running: $(basename "$BIN") $*" >> "$LOG"
"$BIN" "$@" >> "$LOG" 2>&1
echo "exit=$? $(date +%T)" >> "$LOG"
