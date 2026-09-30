#!/bin/bash
# s2fuse_run_one.sh TAG PLINK SPEC...     switches via env FOLD/FUSE/PCAF (0/1)
set -euo pipefail
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
BIN=${BIN:-/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/saige-step2}
TAG=$1; PLINK=$2; shift 2
D=${RUNROOT:-/opt/saige/logs/s2fuse/runs}/$TAG
rm -rf "$D"; mkdir -p "$D/out"
python3 "$(dirname "$0")"/s2fuse_gen_cfg.py "$D" "$PLINK" "$@" >/dev/null
echo "bin $(md5sum "$BIN" | cut -c1-12)  FOLD=${FOLD:-0} FUSE=${FUSE:-0} PCAF=${PCAF:-0} $*" > "$D/meta.txt"
/usr/bin/time -v -o "$D/time.txt" "$BIN" "$D/cfg.yaml" > "$D/log.txt" 2>&1 || { echo FAIL; tail -20 "$D/log.txt"; exit 1; }
grep -E "Elapsed \(wall|Maximum resident" "$D/time.txt"
grep -E "mtFoldQuantProj|mtFuseGemm|mtPopcountAF|MT batch coverage|Sigma p" "$D/log.txt" || true
