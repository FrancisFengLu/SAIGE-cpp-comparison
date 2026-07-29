#!/usr/bin/env bash
# Region-test memory + wall-time benchmark harness.
#
# Sweeps (nThreads, markers_per_chunk_in_groupTest) for the region test,
# captures GNU-time -v Maximum-RSS, wall/user/sys, and writes a CSV.
#
# Usage:
#   bench_region_mem.sh <template.yaml> <out_dir> [reps]
#     template.yaml must have: nThreads / markers_per_chunk_in_groupTest keys
#                              and outputFile pointing anywhere (rewritten per-cell)
#
# Env overrides:
#   THREADS_LIST="1 2 4 8"           (default)
#   CHUNK_LIST="100 500 1000"        (default)
#   BENCH_SAIGE_BIN=/path/to/saige-step2  (default: sibling of this script's step2 dir)
#   BENCH_R_HOME=<R.home>            (only needed if binary embeds R; step2 does not)

set -euo pipefail

TPL="${1:?usage: bench_region_mem.sh <template.yaml> <out_dir> [reps]}"
OUT="${2:?}"
REPS="${3:-1}"
THREADS_LIST="${THREADS_LIST:-1 2 4 8}"
CHUNK_LIST="${CHUNK_LIST:-100 500 1000}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN="${BENCH_SAIGE_BIN:-$SCRIPT_DIR/../saige-step2}"
[[ -x "$BIN" ]] || { echo "ERR: saige-step2 not executable at $BIN" >&2; exit 2; }
[[ -f "$TPL" ]] || { echo "ERR: template not found: $TPL" >&2; exit 2; }
mkdir -p "$OUT"

CSV="$OUT/bench.csv"
echo "rep,nThreads,markers_per_chunk,exit,wall_s,user_s,sys_s,max_rss_kb,max_rss_gb,n_regions" > "$CSV"

for rep in $(seq 1 "$REPS"); do
  for T in $THREADS_LIST; do
    for C in $CHUNK_LIST; do
      cell="T${T}_C${C}_rep${rep}"
      cfg="$OUT/cfg_${cell}.yaml"
      log="$OUT/log_${cell}.txt"
      tim="$OUT/time_${cell}.txt"
      outfile="$OUT/out_${cell}.txt"

      # Rewrite nThreads, markers_per_chunk, outputFile in template
      sed -e "s|^nThreads:.*|nThreads: ${T}|" \
          -e "s|^markers_per_chunk_in_groupTest:.*|markers_per_chunk_in_groupTest: ${C}|" \
          -e "s|^outputFile:.*|outputFile: ${outfile}|" \
          "$TPL" > "$cfg"

      echo "[bench] $cell ..." >&2
      set +e
      /usr/bin/time -v "$BIN" "$cfg" > "$log" 2> "$tim"
      ec=$?
      set -e

      wall=$(grep -oP 'Elapsed \(wall clock\).*: \K[0-9:.]+' "$tim" | tail -1)
      user=$(grep -oP 'User time \(seconds\): \K[0-9.]+'      "$tim" | tail -1)
      sys=$( grep -oP 'System time \(seconds\): \K[0-9.]+'    "$tim" | tail -1)
      rss=$( grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$tim" | tail -1)

      # wall as float seconds
      if [[ "$wall" == *:*:* ]]; then
        wall_s=$(awk -F: '{print $1*3600+$2*60+$3}' <<<"$wall")
      elif [[ "$wall" == *:* ]]; then
        wall_s=$(awk -F: '{print $1*60+$2}' <<<"$wall")
      else
        wall_s="$wall"
      fi
      rss_gb=$(awk -v k="${rss:-0}" 'BEGIN{printf "%.3f", k/1024/1024}')
      nreg=$(grep -oP 'Total regions processed:\s*\K[0-9]+' "$log" | tail -1 || echo 0)

      echo "$rep,$T,$C,$ec,${wall_s:-NA},${user:-NA},${sys:-NA},${rss:-NA},${rss_gb},${nreg:-0}" >> "$CSV"
    done
  done
done

echo "[bench] wrote $CSV" >&2
column -s, -t "$CSV"
