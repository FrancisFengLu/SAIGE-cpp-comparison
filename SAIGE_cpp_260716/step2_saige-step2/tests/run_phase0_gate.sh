#!/bin/bash
# run_phase0_gate.sh -- byte-identity gate for the two CPU-side binary-trait
# switches added for S2_BINARY_GPU phase 0:
#
#   spaScratch               gpos/gneg once per pair + thread_local SPA buffers
#   mtPopcountCtrlFromTotal  control code counts = marker counts - case counts
#
# Three runs per case: the baseline binary (switches do not exist), the new
# binary with the switches absent, the new binary with both on. Every output
# file must be byte-identical across the three. The cases cover the SPA fast
# variant (MAF < 0.29) and the full-N variant, mean-imputed missing calls (the
# popcount replay path), Firth, isMoreOutput, ER on MAC <= 4, MAC 5-20, and
# case rates of 1 / 5 / 10 / ~52 %.
#
# usage: tests/run_phase0_gate.sh [-n NEW_BIN] [-b BASE_BIN] [-w WORKDIR] [-c "CASE ..."]
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEW="$HERE/../saige-step2"
BASE="/opt/saige/logs/s2bingpu/bin/saige-step2.base"
WORK="/opt/saige/logs/s2bingpu/gates/phase0"
CASES="mid_b8 g50k_imb7 rare520_imb grare500_imb"
while getopts "n:b:w:c:" o; do case $o in
  n) NEW=$OPTARG;; b) BASE=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;;
esac; done
NEW="$(cd "$(dirname "$NEW")" && pwd)/$(basename "$NEW")"
[ -x "$NEW" ] || { echo "FAIL: no binary $NEW"; exit 1; }
[ -x "$BASE" ] || { echo "FAIL: no baseline $BASE"; exit 1; }
BAL="/opt/saige/logs/binsplit/models/spa"            # 32 binary models, ~52% cases, mid samples
IMB="/opt/saige/logs/spagaps/step1/out"              # 7 binary models, 1/5/10% cases, mid samples
IMBY="c01_1 c01_2 c05_1 c05_2 c10_1 c10_2 c10_causal"
mkdir -p "$WORK"
FAIL=0
echo "workdir  $WORK"; echo "new      $NEW"; echo "baseline $BASE"; echo "cases    $CASES"; echo

head_yaml() {   # head_yaml <bed> <extra keys...>
  local BED=$1; shift
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "mtPopcountAF: true"
  for L in "$@"; do echo "$L"; done
  echo "models:"
}
one_model() {   # one_model <name> <modeldir> <vrfile> <outdir>
  echo "  - traitName: $1"; echo "    modelFile: $2"
  echo "    varianceRatioFile: $3"; echo "    outputFile: $4/$1.txt"
}
bal_models() { local P=$1 OD=$2 k; for k in $(seq 1 "$P"); do one_model "y$k" "$BAL/m/y$k" "$BAL/mvr_y$k.varianceRatio.txt" "$OD"; done; }
imb_models() { local OD=$1 y; for y in $IMBY; do one_model "$y" "$IMB/m/$y" "$IMB/mvr_$y.varianceRatio.txt" "$OD"; done; }

cfg_for() {     # cfg_for <case> <outdir> <switches: on|off>
  local C=$1 OD=$2 SW=$3 X=()
  [ "$SW" = on ] && X=("spaScratch: true" "mtPopcountCtrlFromTotal: true")
  case "$C" in
    mid_b8)      head_yaml /opt/saige/data/mid "isMoreOutput: true" "isFirth: true" "is_Firth_beta: true" "pCutoffforFirth: 0.05" "${X[@]}"
                 bal_models 8 "$OD" ;;
    g50k_imb7)   head_yaml /opt/saige/logs/binsplit/data/g50k "${X[@]}"
                 imb_models "$OD" ;;
    rare520_imb) head_yaml /opt/saige/logs/gpuassess/rare520 "isFirth: true" "is_Firth_beta: true" "pCutoffforFirth: 0.05" "isMoreOutput: true" "${X[@]}"
                 imb_models "$OD" ;;
    grare500_imb) head_yaml /opt/saige/logs/binsplit/data/grare500 "isFirth: true" "is_Firth_beta: true" "pCutoffforFirth: 0.05" "${X[@]}"
                 imb_models "$OD" ;;
    *) echo "unknown case $C" >&2; return 1 ;;
  esac
}

run_one() {     # run_one <bin> <case> <variant> <switches>
  local BIN=$1 C=$2 V=$3 SW=$4 D="$WORK/$2/$3"
  rm -rf "$D"; mkdir -p "$D/out"
  cfg_for "$C" "$D/out" "$SW" > "$D/cfg.yaml" || return 1
  ( cd "$D" && /usr/bin/time -f '%e' -o wall "$BIN" cfg.yaml > log.txt 2>&1 ); local rc=$?
  echo "    $V rc=$rc wall=$(cat "$D/wall" 2>/dev/null)s files=$(ls "$D/out" | wc -l)"
  return $rc
}

for C in $CASES; do
  echo "== $C"
  run_one "$BASE" "$C" base off || FAIL=1
  run_one "$NEW"  "$C" off  off || FAIL=1
  run_one "$NEW"  "$C" on   on  || FAIL=1
  grep -q 'spaScratch: on' "$WORK/$C/on/log.txt" || { echo "    FAIL: spaScratch not reported on"; FAIL=1; }
  grep -q 'mtPopcountCtrlFromTotal: on' "$WORK/$C/on/log.txt" || { echo "    FAIL: mtPopcountCtrlFromTotal not reported on"; FAIL=1; }
  nd=0
  for f in "$WORK/$C/base/out/"*.txt; do
    b=$(basename "$f")
    cmp -s "$f" "$WORK/$C/off/out/$b" || { echo "    DIFF base vs off: $b"; nd=$((nd+1)); }
    cmp -s "$f" "$WORK/$C/on/out/$b"  || { echo "    DIFF base vs on:  $b"; nd=$((nd+1)); }
  done
  [ "$nd" = 0 ] && echo "    byte-identical: $(ls "$WORK/$C/base/out" | wc -l) files x 3 runs" || FAIL=1
  grep -h 'mtPopcountAF:' "$WORK/$C/on/log.txt" | tail -1 | sed 's/^/    /'
  grep -h 'Firth approx' "$WORK/$C/on/log.txt" | head -2 | sed 's/^/    /'
done
echo
[ "$FAIL" = 0 ] && echo "PHASE0 GATE: PASS" || echo "PHASE0 GATE: FAIL"
exit $FAIL
