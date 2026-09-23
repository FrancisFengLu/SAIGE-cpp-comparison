#!/bin/bash
# run_one.sh TAG PLINK FOLD SPEC...
set -euo pipefail
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
BIN=${BIN:-/opt/saige/logs/mtfold/bin/saige-step2}
TAG=$1; PLINK=$2; FOLD=$3; shift 3
D=${RUNROOT:-/opt/saige/logs/mtfold/runs}/$TAG
rm -rf "$D"; mkdir -p "$D/out"
python3 "$(dirname "$0")"/mtfold_gen_cfg.py "$D" "$PLINK" "$FOLD" "$@" >/dev/null
/usr/bin/time -v -o "$D/time.txt" "$BIN" "$D/cfg.yaml" > "$D/log.txt" 2>&1 || { echo FAIL; tail -20 "$D/log.txt"; exit 1; }
grep -E "Elapsed \(wall|Maximum resident" "$D/time.txt"
grep -E "mtFoldQuantProj|MT batch coverage|Sigma p" "$D/log.txt" || true
