#!/bin/bash
# run_devstats_gate.sh -B BASE_BIN -G NEW_BIN -S NEW_SGS2TXT -w WORKDIR [-c "CASES"]
# Byte-identity gates for gpuDeviceStats (A) and sgsRawDouble (B), S2_DEVSTATS.md.
# Per case every variant writes to the SAME output path <case>/out (moved to
# <case>/<variant>/out after the run: an sgs file's header names the markers
# file's path), with SAIGE_STEP2_ROUTE_DUMP. Keys: useGPU: true only (the
# defaults of origin/main), so gpuBinary / gpuSpa / gpuFirth / gpuOverlap /
# gpuER are all on. Variants and comparisons:
#   base   origin/main binary                      | base vs off: switches absent == origin/main
#   off    new binary, switches absent             |
#   A      + gpuDeviceStats: true                  | off vs A: output + routes byte-identical
#   Ap32   + gpuDeviceStats + gpuPrecisionScan fp32, against offp32 (selected cases)
#   sgs cases only:
#   B      + sgsRawDouble: true                    | off vs B / A vs AB: routes; sgs differ by design
#   AB     + both                                  |
#   Atext  A with outputFormat: text               | sgs2txt(off | B | AB) == Atext, byte for byte
# Cases: bt:<ref> (bingpu_test/ref/<ref>, or the s2-sparse-missing template for
# bm_sparse_nofast_* / mix_*), g200k:<P>:f<0|1>:<text|sgs>. Caller holds the
# pause / lock.
set -u
BASE=""; NEW=""; S2T=""; WORK=""
CASES="bt:full_f0_text bt:full_f1_text bt:sparse_fast_f0_text bt:sparse_fast_f1_text bt:sparse_nofast_f1_text bt:bm_full_f1_text bt:bm_sparse_fast_f1_text bt:mix_full_f1_text bt:full_f1_sgs bt:sparse_fast_f1_sgs bt:bm_sparse_fast_f1_sgs g200k:8:f0:text g200k:8:f1:text g200k:8:f0:sgs g200k:128:f0:text g200k:128:f1:text g200k:128:f0:sgs"
while getopts "B:G:S:w:c:" o; do case $o in B) BASE=$OPTARG;; G) NEW=$OPTARG;; S) S2T=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;; esac; done
for b in "$BASE" "$NEW" "$S2T"; do [ -x "$b" ] || { echo "FAIL: no executable '$b'"; exit 1; }; done
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
BT=/opt/saige/data/bingpu_test
TMPL=/opt/saige/logs/s2-sparse-missing/tmpl
BAL=/opt/saige/logs/binsplit/models/spa
mkdir -p "$WORK"; FAIL=0
echo "base $BASE"; echo "new  $NEW"; echo "sgs2txt $S2T"; echo "cases $CASES"; echo
cfg_bt() {   # cfg_bt <ref> <outdir> keys...
  local R=$1 OD=$2; shift 2
  local CFG="$BT/ref/$R/cfg.yaml"; [ -e "$TMPL/$R/cfg.yaml" ] && CFG="$TMPL/$R/cfg.yaml"
  [ -e "$BT/qm/cfgs/$R/cfg.yaml" ] && CFG="$BT/qm/cfgs/$R/cfg.yaml"
  sed -n '1,/^models:/p' "$CFG" | sed '$d' | sed 's/^nThreads: .*/nThreads: 8/'
  for L in "$@"; do echo "$L"; done
  echo "models:"
  sed -n '/^models:/,$p' "$CFG" | sed 1d | sed -e "s#$BT/ref/[a-z_]*_f[01]_[a-z]*/out/#$OD/#" -e "s#$BT/qm/ref/[a-z0-9_]*/out/#$OD/#"
}
cfg_g200k() {   # cfg_g200k <P> <firth 0|1> <text|sgs> <outdir> keys...
  local P=$1 F=$2 FMT=$3 OD=$4; shift 4
  echo "genoType: plink"; echo "plinkFile: /opt/saige/logs/binsplit/data/g200k"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"
  if [ "$F" = 1 ]; then echo "isFirth: true"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.05"
  else echo "isFirth: false"; echo "is_Firth_beta: false"; fi
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "outputFormat: $FMT"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  local k
  for k in $(seq 1 "$P"); do
    local m=$(( (k-1) % 32 + 1 ))
    echo "  - traitName: y$k"; echo "    modelFile: $BAL/m/y$m"
    echo "    varianceRatioFile: $BAL/mvr_y$m.varianceRatio.txt"; echo "    outputFile: $OD/y$k.txt"
  done
}
run_one() {   # run_one <bin> <casedir> <variant> <gen> <gen args...>
  local BIN=$1 CD=$2 V=$3 GEN=$4; shift 4
  local D="$CD/$V"
  rm -rf "$D" "$CD/out"; mkdir -p "$D/routes" "$CD/out"
  $GEN "$@" > "$D/cfg.yaml" || return 1
  ( cd "$D" && SAIGE_STEP2_ROUTE_DUMP="$D/routes" /usr/bin/time -f '%e %P %MkB' -o wall "$BIN" cfg.yaml > log.txt 2>&1 ); local rc=$?
  mv "$CD/out" "$D/out"
  echo "    $V rc=$rc $(cat "$D/wall" 2>/dev/null) $(grep -h 'useGPU: refused' "$D/log.txt" | head -1) $(grep -h '^  gpuDeviceStats: [0-9]' "$D/log.txt" | sed 's/^ *//' | cut -c1-230)"
  grep -h 'gpuDeviceStats: off\|gpuDeviceStats: on -- S\|self-check' "$D/log.txt" | sed 's/^ */      /' | cut -c1-330
  [ $rc = 0 ] || FAIL=1
}
cmp_dirs() {
  local A=$1 B=$2 LBL=$3 n=0 d=0 f
  for f in "$A"/*; do [ -e "$f" ] || continue; n=$((n+1)); cmp -s "$f" "$B/$(basename "$f")" || { d=$((d+1)); echo "    DIFF $LBL: $(basename "$f")"; }; done
  for f in "$B"/*; do [ -e "$f" ] || continue; [ -e "$A/$(basename "$f")" ] || { d=$((d+1)); echo "    DIFF $LBL: $(basename "$f") only in B"; }; done
  [ $n -gt 0 ] || { d=$((d+1)); echo "    DIFF $LBL: no files in $A"; }
  if [ $d = 0 ]; then echo "    $LBL: $n files byte-identical"; else FAIL=1; fi
}
pair() {  # pair <CD> <A> <B>
  cmp_dirs "$1/$2/out" "$1/$3/out" "$2 vs $3 (output)"; cmp_dirs "$1/$2/routes" "$1/$3/routes" "$2 vs $3 (routes)"
}
to_text() {  # to_text <CD> <variant>: sgs2txt of <variant>/out into <variant>/txt (headers name <CD>/out/*.txt)
  local CD=$1 V=$2 D="$1/$2"
  rm -rf "$CD/out"; mkdir -p "$CD/out"
  local M; M=$(ls "$D"/out/*.markers.sgs | head -1)
  local t0; t0=$(date +%s.%N)
  "$S2T" -j 8 -m "$M" "$D"/out/*.txt.sgs > "$D/sgs2txt.log" 2>&1; local rc=$?
  local t1; t1=$(date +%s.%N)
  rm -rf "$D/txt"; mv "$CD/out" "$D/txt"
  echo "    sgs2txt $V rc=$rc $(echo "$t1 - $t0" | bc) s"
  [ $rc = 0 ] || FAIL=1
}
for C in $CASES; do
  echo "== $C"
  CD="$WORK/$(echo "$C" | tr ':' '_')"; mkdir -p "$CD"
  SGS=0
  case "$C" in
    bt:*) R=${C#bt:}; GEN=cfg_bt; GA=("$R" "$CD/out"); F=0; [[ "$R" == *_f1_* ]] && F=1; [[ "$R" == *_sgs ]] && SGS=1;;
    g200k:*) P=$(echo "$C" | cut -d: -f2); F=$(echo "$C" | cut -d: -f3 | tr -d f); FMT=$(echo "$C" | cut -d: -f4); GEN=cfg_g200k; GA=("$P" "$F" "$FMT" "$CD/out"); [ "$FMT" = sgs ] && SGS=1;;
    *) echo "    unknown case"; FAIL=1; continue;;
  esac
  K=("useGPU: true")
  run_one "$BASE" "$CD" base $GEN "${GA[@]}" "${K[@]}"
  run_one "$NEW"  "$CD" off  $GEN "${GA[@]}" "${K[@]}"
  run_one "$NEW"  "$CD" A    $GEN "${GA[@]}" "${K[@]}" "gpuDeviceStats: true"
  grep -h '^  gate:\|device SPA:\|device Firth:\|GPU coverage' "$CD/A/log.txt" | sed 's/^/    /'
  pair "$CD" base off; pair "$CD" off A
  case "$C" in bt:full_f1_text|bt:sparse_fast_f1_text|g200k:8:f1:text)
    run_one "$NEW" "$CD" offp32 $GEN "${GA[@]}" "${K[@]}" "gpuPrecisionScan: fp32"
    run_one "$NEW" "$CD" Ap32   $GEN "${GA[@]}" "${K[@]}" "gpuPrecisionScan: fp32" "gpuDeviceStats: true"
    pair "$CD" offp32 Ap32
    run_one "$NEW" "$CD" offi8 $GEN "${GA[@]}" "${K[@]}" "gpuPrecisionScan: int8"
    run_one "$NEW" "$CD" Ai8   $GEN "${GA[@]}" "${K[@]}" "gpuPrecisionScan: int8" "gpuDeviceStats: true"
    pair "$CD" offi8 Ai8;;
  esac
  if [ "$SGS" = 1 ]; then
    run_one "$NEW" "$CD" B  $GEN "${GA[@]}" "${K[@]}" "sgsRawDouble: true"
    run_one "$NEW" "$CD" AB $GEN "${GA[@]}" "${K[@]}" "sgsRawDouble: true" "gpuDeviceStats: true"
    case "$C" in
      bt:*)    TR=${R%_sgs}_text; run_one "$NEW" "$CD" Atext cfg_bt "$TR" "$CD/out" "${K[@]}" "gpuDeviceStats: true";;
      g200k:*) run_one "$NEW" "$CD" Atext cfg_g200k "$P" "$F" text "$CD/out" "${K[@]}" "gpuDeviceStats: true";;
    esac
    cmp_dirs "$CD/off/routes" "$CD/B/routes" "off vs B (routes)"
    cmp_dirs "$CD/A/routes" "$CD/AB/routes" "A vs AB (routes)"
    cmp_dirs "$CD/A/routes" "$CD/Atext/routes" "A vs Atext (routes)"
    to_text "$CD" off; to_text "$CD" B; to_text "$CD" AB
    cmp_dirs "$CD/off/txt" "$CD/Atext/out" "sgs2txt(off) vs Atext"
    cmp_dirs "$CD/B/txt" "$CD/Atext/out" "sgs2txt(B) vs Atext"
    cmp_dirs "$CD/AB/txt" "$CD/Atext/out" "sgs2txt(AB) vs Atext"
    # the raw p-values are doubles, so an .sgs under B holds full-precision p (E_PVALRAW); size for the record
    echo "    sgs bytes: off $(du -sb "$CD/off/out" | cut -f1), B $(du -sb "$CD/B/out" | cut -f1), AB $(du -sb "$CD/AB/out" | cut -f1)"
  fi
  # keep logs, drop bulky outputs
  rm -rf "$CD"/*/out "$CD"/*/txt
done
echo
[ "$FAIL" = 0 ] && echo "DEVSTATS GATE: PASS" || echo "DEVSTATS GATE: FAIL"
exit $FAIL
