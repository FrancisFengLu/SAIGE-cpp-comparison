#!/bin/bash
# run_precision_gate.sh BIN_DIR WORKDIR KEY=VALUE [KEY=VALUE ...]
#
# Accuracy gate of one per-stage GPU precision setting (gpu/gpu_precision.hpp)
# against the all-fp64 GPU run of the same build. For each data set below it
# runs BIN_DIR/saige-step2 twice -- useGPU with the defaults (all fp64), and
# with the given keys on top -- writing text output and route dumps
# (SAIGE_STEP2_ROUTE_DUMP), then tools/precision_compare.py prints the overall
# differences (p.value / BETA / SE, crossings of 5e-8 and 1e-5, Is.SPA, Firth
# route and non-converged changes) and writes WORKDIR/<set>.<TAG>.json.
#
#   bt_full     bingpu_test bt, dense models, Firth on, first 8 traits
#   bt_sparse   bingpu_test bt, sparse GRM, fast test off, Firth on, 8 traits
#   bm_full     bingpu_test bm, dense models, Firth on, 8 traits
#   bm_sparse   bingpu_test bm, sparse GRM, fast test on, Firth on, 8 traits
#   rare520     bt dense models on the rare520 genotypes, MACCutoffforER 20,
#               Firth on, 8 traits
#   g200k       g200k (200,000 markers), 8 binary traits, Firth on
#
# Environment:
#   SETS="bt_full g200k"  restrict the sets
#   TAG=name              name of the test run directory (default test); the
#                         all-fp64 reference (WORKDIR/<set>/<REF_TAG>) is reused
#                         when it already ran with rc 0, so several settings can
#                         share one reference
#   REF_TAG / REF_KEYS    a different reference, e.g. REF_TAG=fp64_own
#                         REF_KEYS="gpuSpaImpl=own" for the own SPA kernel
#   MAIN_BIN=DIR          also run DIR/saige-step2 (e.g. an origin/main build)
#                         with the reference config and report whether its text
#                         output and route dumps are byte-identical to the
#                         reference's
#
# Simulated data only. Example:
#   tests/run_precision_gate.sh ../bin /opt/saige/logs/x/gate gpuPrecisionSPA=fp32
# A mode that is not implemented stops the run; its log says
# "useGPU: <stage> precision <mode> is not implemented yet".
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(cd "${1:?usage: run_precision_gate.sh BIN_DIR WORKDIR KEY=VALUE ...}" && pwd)
W=$(mkdir -p "${2:?WORKDIR}" && cd "$2" && pwd)
shift 2
CMP=$HERE/../tools/precision_compare.py
TAG=${TAG:-test}
REF_TAG=${REF_TAG:-fp64}
read -r -a REFK <<< "${REF_KEYS:-}"
BT=/opt/saige/data/bingpu_test/ref
declare -A SRC=(
  [bt_full]=$BT/full_f1_text/cfg.yaml
  [bt_sparse]=$BT/sparse_nofast_f1_text/cfg.yaml
  [bm_full]=$BT/bm_full_f1_text/cfg.yaml
  [bm_sparse]=$BT/bm_sparse_fast_f1_text/cfg.yaml
  [rare520]=$BT/full_f1_text/cfg.yaml
  [g200k]=/opt/saige/logs/s2-pgen/runs/gate_v1/g200k_8_f1/gpu_bed/cfg.yaml
)
declare -A SETKEYS=(
  [rare520]="plinkFile=/opt/saige/logs/gpuassess/rare520 MACCutoffforER=20"
)
SETS=${SETS:-"bt_full bt_sparse bm_full bm_sparse rare520 g200k"}

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
run() {     # BIN DIR -> rc
  rm -rf "$2/routes"; mkdir -p "$2/routes"
  ( cd "$2" && SAIGE_STEP2_ROUTE_DUMP="$2/routes" /usr/bin/time -f '%e' -o wall \
      "$1/saige-step2" cfg.yaml > log.txt 2>&1 )
  local rc=$?
  echo $rc > "$2/rc"
  return $rc
}
same_tree() {   # A B SUBDIR -> "n files, k differ"
  local n=0 k=0 f
  for f in "$1/$3"/*; do
    n=$((n+1)); cmp -s "$f" "$2/$3/$(basename "$f")" || k=$((k+1))
  done
  echo "$n files, $k differ"
}

nfail=0
for s in $SETS; do
  src=${SRC[$s]:-}
  [ -n "$src" ] && [ -f "$src" ] || { echo "== $s: no source config ($src)"; nfail=$((nfail+1)); continue; }
  read -r -a SK <<< "${SETKEYS[$s]:-}"
  R=$W/$s/$REF_TAG
  echo "== $s [$TAG] ($*)"
  if [ "$(cat "$R/rc" 2>/dev/null)" != 0 ]; then
    rm -rf "$R"; mkcfg "$src" "$R" "${SK[@]}" "${REFK[@]}"
    run "$BIN" "$R" || { echo "   reference run failed: $(grep -m1 ERROR "$R/log.txt")"; nfail=$((nfail+1)); continue; }
  fi
  if [ -n "${MAIN_BIN:-}" ]; then
    M=$W/$s/main_$REF_TAG
    rm -rf "$M"; mkcfg "$src" "$M" "${SK[@]}" "${REFK[@]}"
    if run "$(cd "$MAIN_BIN" && pwd)" "$M"; then
      echo "   MAIN_BIN vs $REF_TAG: text $(same_tree "$R" "$M" out); routes $(same_tree "$R" "$M" routes)"
    else
      echo "   MAIN_BIN run failed: $(grep -m1 ERROR "$M/log.txt")"; nfail=$((nfail+1))
    fi
  fi
  [ $# -ge 1 ] || continue
  T=$W/$s/$TAG
  rm -rf "$T"; mkcfg "$src" "$T" "${SK[@]}" "${REFK[@]}" "$@"
  if ! run "$BIN" "$T"; then
    echo "   test run failed: $(grep -m1 ERROR "$T/log.txt")"; nfail=$((nfail+1)); continue
  fi
  grep -h "^GPU precision" "$T/log.txt" | sed 's/^/   /'
  grep -h "useGPU: refused" "$T/log.txt" | sed 's/^/   /'
  python3 "$CMP" "$R" "$T" --quiet --json "$W/$s.$TAG.json" | sed 's/^/   /'
done
[ "$nfail" = 0 ] || { echo "$nfail set(s) did not run"; exit 1; }
