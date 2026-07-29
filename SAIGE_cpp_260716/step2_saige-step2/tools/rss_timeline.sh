#!/usr/bin/env bash
# Sample /proc/<pid>/status VmRSS while running saige-step2, and interleave the
# samples with the program's own stdout lines so you can see WHICH PHASE grows.
#
# Usage: rss_timeline.sh <config.yaml> <out_prefix> [sample_ms]
set -uo pipefail

CFG="${1:?usage: rss_timeline.sh <config.yaml> <out_prefix> [sample_ms]}"
PFX="${2:?}"
MS="${3:-20}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN="${BENCH_SAIGE_BIN:-$SCRIPT_DIR/../saige-step2}"

TS0=$(date +%s.%N)
# process substitution keeps $! = the binary's PID while timestamping each line
stdbuf -oL -eL "$BIN" "$CFG" 2>&1 \
  > >(stdbuf -oL awk -v t0="$TS0" \
        '{ "date +%s.%N" | getline now; close("date +%s.%N");
           printf "%8.3f  %s\n", now-t0, $0; fflush() }' > "$PFX.stdout") &
PID=$!

: > "$PFX.rss"
START=$(date +%s.%N)
while kill -0 "$PID" 2>/dev/null; do
  if [[ -r /proc/$PID/status ]]; then
    rss=$(grep -m1 '^VmRSS:' /proc/$PID/status 2>/dev/null | awk '{print $2}')
    hwm=$(grep -m1 '^VmHWM:' /proc/$PID/status 2>/dev/null | awk '{print $2}')
    now=$(date +%s.%N)
    t=$(awk -v a="$now" -v b="$START" 'BEGIN{printf "%.3f", a-b}')
    [[ -n "${rss:-}" ]] && echo "$t $rss $hwm" >> "$PFX.rss"
  fi
  sleep "$(awk -v m="$MS" 'BEGIN{print m/1000}')"
done
wait "$PID"; EC=$?

echo "exit=$EC" >&2
echo "--- RSS timeline (GB) ---"
awk '{printf "%7.3fs  rss=%7.3f GB  hwm=%7.3f GB\n", $1, $2/1048576, $3/1048576}' "$PFX.rss" \
  | awk 'NR==1||NR%5==0||/hwm/' | tail -60
echo
echo "--- peak ---"
awk 'BEGIN{m=0} {if($3>m)m=$3} END{printf "VmHWM = %.3f GB\n", m/1048576}' "$PFX.rss"
