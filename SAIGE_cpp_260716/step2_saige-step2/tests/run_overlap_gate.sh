#!/bin/bash
# gate.sh -B BASE_BIN -G NEW_BIN -w WORKDIR [-c "CASES"]
# Byte-identity gates for gpuOverlap. Per case, every variant writes to the SAME
# output path <case>/out (moved to <case>/<variant>/out after the run: an sgs
# file's header names the markers file's path), with SAIGE_STEP2_ROUTE_DUMP.
# Then every output file (text or sgs incl. markers.sgs) and every route file
# must be byte-identical between: base vs off (switch absent: new binary == origin/main),
# off vs on, off vs on2, and each extra pair (offX vs onX). Caller holds the lock.
#   ALL keys: useGPU gpuBinary gpuSpa(lib) gpuOwnSampleSets gpuSparse gpuPrefetch(4 readers)
#             parallelModelLoad + the five kernel-tune switches; + gpuFirth when Firth is on
#   on  = ALL + gpuOverlap
#   on2 = ALL - gpuPrefetch + gpuOverlap, 3 sets, lag 1 (scan worker reads; tightest ring)
#   _nf = gpuFirth off (Firth cases), _min = useGPU/gpuBinary/gpuSpa with every other GPU switch written off,
#   _def = only useGPU: true (the defaults of 69a2c3b8), _own = gpuSpaImpl own,
#   _sb = gpuBlockSize 1024 (many small superblocks; on_sb with 4 sets, lag 3)
# Cases: bt:<ref> (bingpu_test/ref/<ref>, or the s2-sparse-missing template for bm_sparse_nofast_*),
#        g200k:<P>:f<0|1> (sgs)
set -u
BASE=""; NEW=""; WORK=""
CASES="bt:full_f0_text bt:full_f1_text bt:full_f0_sgs bt:sparse_fast_f0_text bt:sparse_fast_f1_text bt:sparse_nofast_f0_text bt:sparse_nofast_f1_text bt:bm_full_f0_text bt:bm_full_f1_text bt:bm_sparse_fast_f0_text bt:bm_sparse_fast_f1_text bt:bm_sparse_nofast_f1_text bt:bm_sparse_fast_f1_sgs bt:mix_full_f1_text bt:mix_sparse_fast_f1_text bt:qm_sparse_fast_f0_text bt:qm_full_f0_sgs g200k:8:f0 g200k:8:f1 g200k:32:f0 g200k:32:f1 g200k:128:f0 g200k:128:f1"
while getopts "B:G:w:c:" o; do case $o in B) BASE=$OPTARG;; G) NEW=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;; esac; done
for b in "$BASE" "$NEW"; do [ -x "$b" ] || { echo "FAIL: no executable '$b'"; exit 1; }; done
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
BT=/opt/saige/data/bingpu_test
TMPL=/opt/saige/logs/s2-sparse-missing/tmpl
BAL=/opt/saige/logs/binsplit/models/spa
mkdir -p "$WORK"; FAIL=0
echo "base $BASE"; echo "new  $NEW"; echo "cases $CASES"; echo
ALL=("useGPU: true" "gpuBinary: true" "gpuSpa: true" "gpuSpaImpl: lib" "gpuOwnSampleSets: true" "gpuSparse: true"
     "gpuPrefetch: true" "gpuPrefetchThreads: 4" "parallelModelLoad: true"
     "gpuSpaOrder: trait" "gpuSpaDynamic: true" "gpuSpaFused: true" "gpuSpaMinBlocks: 3" "gpuDecodeX2: true")
cfg_bt() {   # cfg_bt <ref> <outdir> keys...
  local R=$1 OD=$2; shift 2
  local CFG="$BT/ref/$R/cfg.yaml"; [ -e "$TMPL/$R/cfg.yaml" ] && CFG="$TMPL/$R/cfg.yaml"
  [ -e "$BT/qm/cfgs/$R/cfg.yaml" ] && CFG="$BT/qm/cfgs/$R/cfg.yaml"
  sed -n '1,/^models:/p' "$CFG" | sed '$d' | sed 's/^nThreads: .*/nThreads: 8/'
  for L in "$@"; do echo "$L"; done
  echo "models:"
  sed -n '/^models:/,$p' "$CFG" | sed 1d | sed -e "s#$BT/ref/[a-z_]*_f[01]_[a-z]*/out/#$OD/#" -e "s#$BT/qm/ref/[a-z0-9_]*/out/#$OD/#"
}
cfg_g200k() {   # cfg_g200k <P> <firth 0|1> <outdir> keys...
  local P=$1 F=$2 OD=$3; shift 3
  echo "genoType: plink"; echo "plinkFile: /opt/saige/logs/binsplit/data/g200k"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"
  if [ "$F" = 1 ]; then echo "isFirth: true"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01"
  else echo "isFirth: false"; echo "is_Firth_beta: false"; fi
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "mtPopcountAF: true"; echo "spaScratch: true"; echo "mtPopcountCtrlFromTotal: true"
  echo "outputFormat: sgs"
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
  echo "    $V rc=$rc $(cat "$D/wall" 2>/dev/null) $(grep -h 'useGPU: refused' "$D/log.txt" | head -1) $(grep -h '\[gpu overlap\]' "$D/log.txt" | sed 's/^ *//' | cut -c1-160)"
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
without() { local x; for x in "${ALL[@]}"; do [[ "$x" == "$1"* ]] || echo "$x"; done; }
for C in $CASES; do
  echo "== $C"
  CD="$WORK/$(echo "$C" | tr ':' '_')"; mkdir -p "$CD"
  case "$C" in
    bt:*) R=${C#bt:}; GEN=cfg_bt; GA=("$R" "$CD/out"); F=0; [[ "$R" == *_f1_* ]] && F=1;;
    g200k:*) P=$(echo "$C" | cut -d: -f2); F=$(echo "$C" | cut -d: -f3 | tr -d f); GEN=cfg_g200k; GA=("$P" "$F" "$CD/out");;
    *) echo "    unknown case"; FAIL=1; continue;;
  esac
  K=("${ALL[@]}"); [ "$F" = 1 ] && K+=("gpuFirth: true")
  mapfile -t KNP < <(without gpuPrefetch); [ "$F" = 1 ] && KNP+=("gpuFirth: true")
  run_one "$BASE" "$CD" base $GEN "${GA[@]}" "${K[@]}"
  run_one "$NEW"  "$CD" off  $GEN "${GA[@]}" "${K[@]}"
  run_one "$NEW"  "$CD" on   $GEN "${GA[@]}" "${K[@]}" "gpuOverlap: true"
  run_one "$NEW"  "$CD" on2  $GEN "${GA[@]}" "${KNP[@]}" "gpuOverlap: true" "gpuOverlapSets: 3" "gpuOverlapLag: 1"
  grep -h '^  gate:\|device SPA:\|device Firth:\|GPU coverage' "$CD/on/log.txt" | sed 's/^/    /'
  pair "$CD" base off; pair "$CD" off on; pair "$CD" off on2
  EXTRA=""
  [ "$F" = 1 ] && [[ "$C" != g200k:128:* ]] && [[ "$C" != g200k:32:* ]] && EXTRA="nf"   # CPU Firth at P>=32: ~1000 s
  case "$C" in bt:full_f0_text|bt:sparse_fast_f1_text|g200k:32:f1) EXTRA="$EXTRA min";; esac
  case "$C" in bt:full_f1_text|g200k:8:f0) EXTRA="$EXTRA own";; esac
  case "$C" in bt:full_f1_text|bt:bm_sparse_fast_f1_text|bt:mix_full_f1_text|g200k:32:f1) EXTRA="$EXTRA def";; esac
  case "$C" in bt:full_f1_text|bt:bm_sparse_fast_f1_text|bt:sparse_nofast_f1_text|bt:bm_full_f0_text) EXTRA="$EXTRA sb";; esac
  for X in $EXTRA; do
    case $X in
      nf)  KX=("${ALL[@]}");;
      min) KX=("useGPU: true" "gpuBinary: true" "gpuSpa: true" "gpuFirth: false" "gpuSparse: false" "gpuOwnSampleSets: false"
               "gpuPrefetch: false" "parallelModelLoad: false" "gpuSpaOrder: marker" "gpuSpaDynamic: false"
               "gpuSpaFused: false" "gpuSpaMinBlocks: 0" "gpuDecodeX2: false");;
      def) KX=("useGPU: true");;
      own) mapfile -t KX < <(without gpuSpaImpl); KX+=("gpuSpaImpl: own"); [ "$F" = 1 ] && KX+=("gpuFirth: true");;
      sb)  KX=("${K[@]}" "gpuBlockSize: 1024");;
    esac
    OX=("gpuOverlap: true"); [ $X = sb ] && OX+=("gpuOverlapSets: 4" "gpuOverlapLag: 3")
    run_one "$NEW" "$CD" off_$X $GEN "${GA[@]}" "${KX[@]}"
    run_one "$NEW" "$CD" on_$X  $GEN "${GA[@]}" "${KX[@]}" "${OX[@]}"
    pair "$CD" off_$X on_$X
  done
  # keep logs, drop bulky outputs
  rm -rf "$CD"/*/out
done
echo
[ "$FAIL" = 0 ] && echo "OVERLAP GATE: PASS" || echo "OVERLAP GATE: FAIL"
exit $FAIL
