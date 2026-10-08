#!/bin/bash
# run_precision_gate.sh BIN_DIR WORKDIR KEY=VALUE [KEY=VALUE ...]
#
# Accuracy gate of one per-stage GPU precision setting (gpu/gpu_precision.hpp)
# against the all-fp64 GPU run of the same build. For each data set below it
# runs BIN_DIR/saige-step2 twice -- useGPU with the defaults (all fp64), and
# with the given keys on top -- writing text output, then
# tools/precision_compare.py prints per-trait and overall differences and
# writes WORKDIR/<set>.json.
#
#   bt_full     bingpu_test bt, dense models, Firth on, first 8 traits
#   bt_sparse   bingpu_test bt, sparse GRM, fast test off, Firth on, 8 traits
#   g200k       g200k (200,000 markers), 8 binary traits, Firth on
#
# SETS="bt_full g200k" restricts the sets. Simulated data only. Example:
#   tests/run_precision_gate.sh ../bin /opt/saige/logs/x/gate_spa gpuPrecisionSPA=fp32
# A mode that is not implemented stops the run; its log says
# "useGPU: <stage> precision <mode> is not implemented yet".
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(cd "${1:?usage: run_precision_gate.sh BIN_DIR WORKDIR KEY=VALUE ...}" && pwd)
W=$(mkdir -p "${2:?WORKDIR}" && cd "$2" && pwd)
shift 2
[ $# -ge 1 ] || { echo "give at least one KEY=VALUE (e.g. gpuPrecisionSPA=fp32)"; exit 2; }
CMP=$HERE/../tools/precision_compare.py
BT=/opt/saige/data/bingpu_test/ref
declare -A SRC=(
  [bt_full]=$BT/full_f1_text/cfg.yaml
  [bt_sparse]=$BT/sparse_nofast_f1_text/cfg.yaml
  [g200k]=/opt/saige/logs/s2-pgen/runs/gate_v1/g200k_8_f1/gpu_bed/cfg.yaml
)
SETS=${SETS:-"bt_full bt_sparse g200k"}

mkcfg() {   # SRC OUTDIR KEY=VALUE...
  python3 - "$@" <<'EOF'
import sys, os, yaml
src, outdir = sys.argv[1], sys.argv[2]
c = yaml.safe_load(open(src))
c['models'] = c['models'][:8]
os.makedirs(os.path.join(outdir, 'out'), exist_ok=True)
for m in c['models']:
    m['outputFile'] = os.path.join(outdir, 'out', m['traitName'] + '.txt')
c['useGPU'] = True
c['outputFormat'] = 'text'
for kv in sys.argv[3:]:
    k, v = kv.split('=', 1)
    c[k] = yaml.safe_load(v)
yaml.safe_dump(c, open(os.path.join(outdir, 'cfg.yaml'), 'w'), sort_keys=False)
EOF
}
run() {     # DIR -> rc
  "$BIN/saige-step2" "$1/cfg.yaml" > "$1/log.txt" 2>&1
  local rc=$?
  echo $rc > "$1/rc"
  return $rc
}

nfail=0
for s in $SETS; do
  src=${SRC[$s]:-}
  [ -n "$src" ] && [ -f "$src" ] || { echo "== $s: no source config ($src)"; nfail=$((nfail+1)); continue; }
  rm -rf "$W/$s"
  mkcfg "$src" "$W/$s/fp64"
  mkcfg "$src" "$W/$s/test" "$@"
  echo "== $s ($*)"
  run "$W/$s/fp64" || { echo "   fp64 reference run failed: $(grep -m1 ERROR "$W/$s/fp64/log.txt")"; nfail=$((nfail+1)); continue; }
  if ! run "$W/$s/test"; then
    echo "   test run failed: $(grep -m1 ERROR "$W/$s/test/log.txt")"; nfail=$((nfail+1)); continue
  fi
  grep -h "^GPU precision" "$W/$s/test/log.txt" | sed 's/^/   /'
  grep -h "useGPU: refused" "$W/$s/test/log.txt" | sed 's/^/   /'
  python3 "$CMP" "$W/$s/fp64" "$W/$s/test" --quiet --json "$W/$s.json" | sed 's/^/   /'
done
[ "$nfail" = 0 ] || { echo "$nfail set(s) did not run"; exit 1; }
