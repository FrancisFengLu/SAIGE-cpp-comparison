#!/bin/bash
# mtvec_run_one.sh TAG PLINK VEC(0|1) SPEC...
# env: BIN, RUNROOT, FOLD, SGS, NTHREADS  (see mtvec_gen_cfg.py)
set -euo pipefail
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
BIN=${BIN:-/opt/saige/logs/mtvec/bin/saige-step2}
TAG=$1; PLINK=$2; VEC=$3; shift 3
D=${RUNROOT:-/opt/saige/logs/mtvec/runs}/$TAG
rm -rf "$D"; mkdir -p "$D/out"
python3 "$(dirname "$0")"/mtvec_gen_cfg.py "$D" "$PLINK" "$VEC" "$@" >/dev/null
/usr/bin/time -v -o "$D/time.txt" "$BIN" "$D/cfg.yaml" > "$D/log.txt" 2>&1 \
  || { echo FAIL; tail -20 "$D/log.txt"; exit 1; }
grep -E "Elapsed \(wall|Maximum resident" "$D/time.txt"
grep -E "mtVecQuantStats|mtFoldQuantProj|MT batch coverage|mtvec prof" "$D/log.txt" || true
