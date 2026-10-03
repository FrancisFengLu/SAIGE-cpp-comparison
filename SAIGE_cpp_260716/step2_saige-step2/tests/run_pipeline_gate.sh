#!/bin/bash
# run_pipeline_gate.sh -- byte-identity gates for the two step-2 pipeline
# switches (S2_PIPELINE.md): gpuPrefetch (read / compute overlap on the GPU
# path) and parallelModelLoad (concurrent null-model loading, both paths).
#
# Per case, the same configuration is run as
#   base   the baseline binary (the commit the branch started from), switches absent
#   off    the new binary, switches absent
#   on     the new binary, gpuPrefetch: true + parallelModelLoad: true
#   on2    the new binary, gpuPrefetch: true (2 sets, 1 reader thread) + parallelModelLoad: true
# with SAIGE_STEP2_ROUTE_DUMP set, and then every output file (text or sgs,
# the markers.sgs included) and every route file must be byte-identical across
# the four runs; the loader section of the log (between "Loading null model"
# and "Setting global variables") must be identical between off and on; a
# bingpu_test case's text output (for an sgs case: its sgs2txt conversion)
# must also equal the harness's CPU reference text.
#
# Cases:
#   bt:<ref>          a /opt/saige/data/bingpu_test reference run (its cfg.yaml, with
#                     nThreads 8 and the GPU keys added; the GPU path refuses bm_* and
#                     sparse_nofast, which then exercise the CPU path's model loading)
#   g200k:<P>:<fmt>   g200k x P balanced models (the phase-3 cells' models), fmt text|sgs;
#                     list a P's text case before its sgs case to get the sgs -> text
#                     round trip between them
#
# usage: tests/run_pipeline_gate.sh -B BASE_GPU_BIN -G NEW_GPU_BIN -X NEW_SGS2TXT \
#            -w WORKDIR [-c "CASE ..."] [-k]      (-k keeps the outputs)
set -u
BASE=""; NEW=""; CVT=""; WORK=""; KEEP=0
CASES="bt:full_f0_text bt:full_f1_text bt:sparse_fast_f0_text bt:sparse_fast_f1_text bt:full_f0_sgs bt:bm_full_f0_text g200k:8:text g200k:8:sgs g200k:32:text g200k:32:sgs g200k:128:sgs"
while getopts "B:G:X:w:c:k" o; do case $o in
  B) BASE=$OPTARG;; G) NEW=$OPTARG;; X) CVT=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;; k) KEEP=1;;
esac; done
for b in "$BASE" "$NEW" "$CVT"; do [ -x "$b" ] || { echo "FAIL: no executable '$b'"; exit 1; }; done
[ -n "$WORK" ] || { echo "FAIL: -w WORKDIR"; exit 1; }
BT=/opt/saige/data/bingpu_test
BAL=/opt/saige/logs/binsplit/models/spa
mkdir -p "$WORK"
FAIL=0
echo "workdir $WORK"; echo "base    $BASE"; echo "new     $NEW"; echo "sgs2txt $CVT"; echo "cases   $CASES"; echo

ONKEYS=("gpuPrefetch: true" "parallelModelLoad: true")
ON2KEYS=("gpuPrefetch: true" "gpuPrefetchSets: 2" "gpuPrefetchThreads: 1" "parallelModelLoad: true")
GPUKEYS=("useGPU: true" "gpuBinary: true" "gpuSpa: true" "gpuSpaImpl: lib")

# cfg_bt <ref> <outdir> <extra keys...>: the reference's cfg.yaml with the
# output paths moved, 8 threads, and the GPU keys (as the lib gate ran it).
cfg_bt() {
  local R=$1 OD=$2; shift 2
  local CFG="$BT/ref/$R/cfg.yaml"
  sed -n '1,/^models:/p' "$CFG" | sed '$d' | sed 's/^nThreads: .*/nThreads: 8/'
  for L in "${GPUKEYS[@]}" "$@"; do echo "$L"; done
  echo "models:"
  sed -n '/^models:/,$p' "$CFG" | sed 1d | sed "s#$BT/ref/$R/out/#$OD/#"
}
# cfg_g200k <P> <fmt> <outdir> <extra keys...>: the phase-3 cell configuration.
cfg_g200k() {
  local P=$1 FMT=$2 OD=$3; shift 3
  echo "genoType: plink"; echo "plinkFile: /opt/saige/logs/binsplit/data/g200k"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "mtPopcountAF: true"; echo "spaScratch: true"; echo "mtPopcountCtrlFromTotal: true"
  echo "outputFormat: $FMT"
  for L in "${GPUKEYS[@]}" "$@"; do echo "$L"; done
  echo "models:"
  local k
  for k in $(seq 1 "$P"); do
    local m=$(( (k-1) % 32 + 1 ))
    echo "  - traitName: y$k"; echo "    modelFile: $BAL/m/y$m"
    echo "    varianceRatioFile: $BAL/mvr_y$m.varianceRatio.txt"; echo "    outputFile: $OD/y$k.txt"
  done
}
run_one() {   # run_one <bin> <casedir> <variant> <cfg-generator> <generator args...>
  local BIN=$1 CD=$2 V=$3 GEN=$4; shift 4
  local D="$CD/$V"
  rm -rf "$D"; mkdir -p "$D/out" "$D/routes"
  $GEN "$@" > "$D/cfg.yaml" || return 1
  ( cd "$D" && SAIGE_STEP2_ROUTE_DUMP="$D/routes" /usr/bin/time -f '%e %P' -o wall "$BIN" cfg.yaml > log.txt 2>&1 ); local rc=$?
  echo "    $V rc=$rc wall=$(cut -d' ' -f1 "$D/wall" 2>/dev/null)s cpu=$(cut -d' ' -f2 "$D/wall" 2>/dev/null) $(grep -h 'useGPU: refused' "$D/log.txt" | head -1) $(grep -h '\[gpu pipeline\]' "$D/log.txt" | sed 's/^ *//')"
  [ $rc = 0 ] || FAIL=1
  return $rc
}
cmp_dirs() {  # cmp_dirs <dirA> <dirB> <label>: the same file names, every file byte-identical
  local A=$1 B=$2 LBL=$3 n=0 d=0 f
  for f in "$A"/*; do [ -e "$f" ] || continue; n=$((n+1)); cmp -s "$f" "$B/$(basename "$f")" || { d=$((d+1)); echo "    DIFF $LBL: $(basename "$f")"; }; done
  for f in "$B"/*; do [ -e "$f" ] || continue; [ -e "$A/$(basename "$f")" ] || { d=$((d+1)); echo "    DIFF $LBL: $(basename "$f") only in B"; }; done
  [ $n -gt 0 ] || { d=$((d+1)); echo "    DIFF $LBL: no files in $A"; }
  if [ $d = 0 ]; then echo "    $LBL: $n files byte-identical"; else FAIL=1; fi
}
loader_log() { sed -n '/===== Loading null model/,/===== Setting global variables/p' "$1" | grep -v '^\[TIMING\]'; }
sgs2txt_dir() {  # sgs2txt_dir <sgsdir> <txtdir>
  rm -rf "$2"; mkdir -p "$2"
  local f b
  for f in "$1"/*.txt.sgs; do [[ "$f" == *".markers.sgs" ]] && continue; b=$(basename "$f" .sgs); "$CVT" -o "$2/$b" "$f" > /dev/null 2>&1 || { echo "    sgs2txt failed: $b"; FAIL=1; }; done
}

for C in $CASES; do
  echo "== $C"
  CD="$WORK/$(echo "$C" | tr ':' '_')"; mkdir -p "$CD"
  case "$C" in
    bt:*)
      R=${C#bt:}
      run_one "$BASE" "$CD" base cfg_bt "$R" "$CD/base/out"
      run_one "$NEW"  "$CD" off  cfg_bt "$R" "$CD/off/out"
      run_one "$NEW"  "$CD" on   cfg_bt "$R" "$CD/on/out"  "${ONKEYS[@]}"
      run_one "$NEW"  "$CD" on2  cfg_bt "$R" "$CD/on2/out" "${ON2KEYS[@]}"
      ;;
    g200k:*)
      P=$(echo "$C" | cut -d: -f2); FMT=$(echo "$C" | cut -d: -f3)
      run_one "$BASE" "$CD" base cfg_g200k "$P" "$FMT" "$CD/base/out"
      run_one "$NEW"  "$CD" off  cfg_g200k "$P" "$FMT" "$CD/off/out"
      run_one "$NEW"  "$CD" on   cfg_g200k "$P" "$FMT" "$CD/on/out"  "${ONKEYS[@]}"
      run_one "$NEW"  "$CD" on2  cfg_g200k "$P" "$FMT" "$CD/on2/out" "${ON2KEYS[@]}"
      ;;
    *) echo "    unknown case $C"; FAIL=1; continue;;
  esac
  grep -h '^  gate:\|^  GPU coverage\|device SPA:' "$CD/on/log.txt" | sed 's/^/    /'
  for V in off on on2; do
    grep -q 'gpuPrefetch: on' "$CD/$V/log.txt" && PF=1 || PF=0
    grep -q 'parallelModelLoad: on' "$CD/$V/log.txt" && PL=1 || PL=0
    echo "    $V: gpuPrefetch=$PF parallelModelLoad=$PL $(grep -h 'TIMING\] parallelModelLoad' "$CD/$V/log.txt" | sed 's/.*parallelModelLoad: //') $(grep -h 'TIMING\] 20_null_model_loaded' "$CD/$V/log.txt" | sed 's/.*loaded//')"
  done
  cmp_dirs "$CD/base/out"    "$CD/off/out"    "base vs off (output)"
  cmp_dirs "$CD/base/routes" "$CD/off/routes" "base vs off (routes)"
  cmp_dirs "$CD/off/out"     "$CD/on/out"     "off vs on (output)"
  cmp_dirs "$CD/off/routes"  "$CD/on/routes"  "off vs on (routes)"
  cmp_dirs "$CD/off/out"     "$CD/on2/out"    "off vs on2 (output)"
  cmp_dirs "$CD/off/routes"  "$CD/on2/routes" "off vs on2 (routes)"
  if diff -q <(loader_log "$CD/off/log.txt") <(loader_log "$CD/on/log.txt") > /dev/null; then
    echo "    loader log off vs on: identical ($(loader_log "$CD/off/log.txt" | wc -l) lines)"
  else
    echo "    DIFF loader log off vs on"; diff <(loader_log "$CD/off/log.txt") <(loader_log "$CD/on/log.txt") | head -20 | sed 's/^/      /'; FAIL=1
  fi
  case "$C" in
    bt:*_text) cmp_dirs "$BT/ref/$R/out" "$CD/on/out" "reference vs on (text)";;
    bt:*_sgs)  sgs2txt_dir "$CD/on/out" "$CD/on/txt"
               cmp_dirs "$BT/ref/${R%_sgs}_text/out" "$CD/on/txt" "reference text vs sgs2txt(on sgs)"
               rm -rf "$CD/on/txt";;
  esac
  # sgs -> text round trip of the on run against this P's text case (kept for it)
  if [[ "$C" == g200k:*:sgs ]]; then
    TXT="$WORK/g200k_${P}_text/on/out"
    if [ -d "$TXT" ] && ls "$TXT"/*.txt > /dev/null 2>&1; then
      sgs2txt_dir "$CD/on/out" "$CD/on/txt"
      cmp_dirs "$TXT" "$CD/on/txt" "sgs2txt(on sgs) vs on text"
      rm -rf "$CD/on/txt"
      [ $KEEP = 0 ] && rm -rf "$TXT"
    fi
  fi
  if [ $KEEP = 0 ]; then
    for V in base off on2; do rm -rf "$CD/$V/out"; done
    [[ "$C" == g200k:*:text ]] || rm -rf "$CD/on/out"
  fi
done
echo
[ "$FAIL" = 0 ] && echo "PIPELINE GATE: PASS" || echo "PIPELINE GATE: FAIL"
exit $FAIL
