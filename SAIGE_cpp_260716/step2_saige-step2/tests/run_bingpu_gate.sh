#!/bin/bash
# run_bingpu_gate.sh -- accuracy gates for the binary-trait GPU path
# (useGPU + gpuBinary [+ gpuSpa]; S2_BINARY_GPU.md).
#
# Per case, four runs of the same configuration, all with the CPU-side
# switches on (spaScratch, mtPopcountAF, mtPopcountCtrlFromTotal):
#   cpu_text   CPU multi-trait loop, text output, route dump
#   gpu_text   GPU path, text output, route dump
#   cpu_sgs    CPU loop, outputFormat: sgs          (full-precision fields)
#   gpu_sgs    GPU path, outputFormat: sgs
# and then:
#   1. tools/sgs2txt(gpu_sgs) must be byte-identical to gpu_text  (the
#      results-level cut: what the CPU box reconstructs is what the GPU box
#      would have printed);
#   2. tools/s2bingpu_cmp.py cpu_text vs gpu_text: exact columns identical,
#      routing identical pair for pair (route dumps), per-field max |delta|,
#      max |delta(-log10 p)|, 5e-8 / 1e-5 crossings; plus the fp64 fields from
#      the two sgs runs.
#
# usage: tests/run_bingpu_gate.sh [-n CPU_BIN] [-g GPU_BIN] [-w WORKDIR]
#                                 [-c "CASE ..."] [-S]   (-S adds gpuSpa: true)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPU="$HERE/../saige-step2"
GPU="$HERE/../saige-step2.cuda"
CVT="$HERE/../tools/sgs2txt"
CMP="$HERE/../tools/s2bingpu_cmp.py"
WORK="/opt/saige/logs/s2bingpu/gates/phase1"
CASES="bal8_g200k imb7_g200k mid_b8 rare520_imb grare500_imb mixed_g50k"
SPA=0
while getopts "n:g:w:c:S" o; do case $o in
  n) CPU=$OPTARG;; g) GPU=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;; S) SPA=1;;
esac; done
for b in "$CPU" "$GPU" "$CVT"; do [ -x "$b" ] || { echo "FAIL: no executable $b"; exit 1; }; done
BAL="/opt/saige/logs/binsplit/models/spa"
IMB="/opt/saige/logs/spagaps/step1/out"
IMBY="c01_1 c01_2 c05_1 c05_2 c10_1 c10_2 c10_causal"
QNT="/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out"
mkdir -p "$WORK"
FAIL=0
echo "workdir $WORK"; echo "cpu     $CPU"; echo "gpu     $GPU"; echo "cases   $CASES"; echo "gpuSpa  $SPA"; echo

# head_yaml <bed> <extra keys...>. yaml-cpp keeps the FIRST of two equal keys,
# so the Firth / isMoreOutput trio is emitted once, from FIRTH / MORE: FIRTH=1
# means isFirth + is_Firth_beta true with pCutoffforFirth 0.05 (the models
# carry 0.01), MORE=1 means isMoreOutput.
FIRTH=0; MORE=0
head_yaml() {
  local BED=$1; shift
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: $([ "$MORE" = 1 ] && echo true || echo false)"
  if [ "$FIRTH" = 1 ]; then echo "isFirth: true"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.05"
  else echo "isFirth: false"; echo "is_Firth_beta: false"; fi
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "mtPopcountAF: true"; echo "spaScratch: true"; echo "mtPopcountCtrlFromTotal: true"
  for L in "$@"; do echo "$L"; done
  echo "models:"
}
one_model() { echo "  - traitName: $1"; echo "    modelFile: $2"; echo "    varianceRatioFile: $3"; echo "    outputFile: $4/$1.txt"; }
bal_models() { local P=$1 OD=$2 k; for k in $(seq 1 "$P"); do one_model "y$k" "$BAL/m/y$k" "$BAL/mvr_y$k.varianceRatio.txt" "$OD"; done; }
imb_models() { local OD=$1 y; for y in $IMBY; do one_model "$y" "$IMB/m/$y" "$IMB/mvr_$y.varianceRatio.txt" "$OD"; done; }

cfg_for() {     # cfg_for <case> <outdir> <extra...>
  local C=$1 OD=$2; shift 2
  case "$C" in
    bal8_g200k)   head_yaml /opt/saige/logs/binsplit/data/g200k "$@"; bal_models 8 "$OD" ;;
    imb7_g200k)   head_yaml /opt/saige/logs/binsplit/data/g200k "$@"; imb_models "$OD" ;;
    imb7_g50k)    head_yaml /opt/saige/logs/binsplit/data/g50k "$@"; imb_models "$OD" ;;
    mid_b8)       FIRTH=1 MORE=1 head_yaml /opt/saige/data/mid "$@"; bal_models 8 "$OD" ;;
    mid_b8_nofirth) FIRTH=0 MORE=1 head_yaml /opt/saige/data/mid "$@"; bal_models 8 "$OD" ;;
    rare520_imb)  FIRTH=1 MORE=1 head_yaml /opt/saige/logs/gpuassess/rare520 "$@"; imb_models "$OD" ;;
    grare500_imb) FIRTH=1 MORE=0 head_yaml /opt/saige/logs/binsplit/data/grare500 "$@"; imb_models "$OD" ;;
    imb7_g50k_firth) FIRTH=1 MORE=0 head_yaml /opt/saige/logs/binsplit/data/g50k "$@"; imb_models "$OD" ;;
    mixed_g50k)   head_yaml /opt/saige/logs/binsplit/data/g50k "$@"
                  local k; for k in 1 2 3 4; do
                    one_model "b$k" "$BAL/m/y$k" "$BAL/mvr_y$k.varianceRatio.txt" "$OD"
                    one_model "q$k" "$QNT/m/y$k" "$QNT/mvr_y$k.varianceRatio.txt" "$OD"
                  done ;;
    *) echo "unknown case $C" >&2; return 1 ;;
  esac
}

run_one() {     # run_one <bin> <case> <variant> <extra yaml...>
  local BIN=$1 C=$2 V=$3; shift 3
  local D="$WORK/$C/$V"
  rm -rf "$D"; mkdir -p "$D/out" "$D/routes"
  cfg_for "$C" "$D/out" "$@" > "$D/cfg.yaml" || return 1
  ( cd "$D" && SAIGE_STEP2_ROUTE_DUMP="$D/routes" /usr/bin/time -f '%e' -o wall "$BIN" cfg.yaml > log.txt 2>&1 ); local rc=$?
  echo "    $V rc=$rc wall=$(cat "$D/wall" 2>/dev/null)s $(grep -h 'useGPU: refused' "$D/log.txt" | head -1)"
  return $rc
}

GPUKEYS=("useGPU: true" "gpuBinary: true")
[ "$SPA" = 1 ] && GPUKEYS+=("gpuSpa: true")

for C in $CASES; do
  echo "== $C"
  run_one "$CPU" "$C" cpu_text || FAIL=1
  run_one "$GPU" "$C" gpu_text "${GPUKEYS[@]}" || FAIL=1
  run_one "$CPU" "$C" cpu_sgs "outputFormat: sgs" || FAIL=1
  run_one "$GPU" "$C" gpu_sgs "outputFormat: sgs" "${GPUKEYS[@]}" || FAIL=1
  if grep -q 'useGPU: refused' "$WORK/$C/gpu_text/log.txt"; then echo "    FAIL: GPU path refused"; FAIL=1; fi
  # A "Firth on" case must have exercised Firth on both paths (see the FIRTH
  # note above head_yaml): the per-trait "[trait] Firth approx was applied to N
  # markers" lines must sum to > 0 (the 2026-10-03 03:21-04:34 runs summed the
  # word "to" -- field 6 -- and printed a false FAIL for every Firth case; the
  # numbers were rechecked from the kept logs), and the GPU gate line must show needFirth.
  if grep -q '^is_Firth_beta: true' "$WORK/$C/cpu_text/cfg.yaml"; then
    for V in cpu_text gpu_text; do
      nf=$(grep -h 'Firth approx was applied to' "$WORK/$C/$V/log.txt" | grep -oP 'applied to \K[0-9]+' | awk '{s+=$1} END{print s+0}')
      [ "$nf" -gt 0 ] && echo "    $V: Firth applied to $nf (marker, trait) pairs" || { echo "    FAIL: $V ran with the Firth path never taken"; FAIL=1; }
    done
    grep -h '^  gate:' "$WORK/$C/gpu_text/log.txt" | grep -q 'needFirth [1-9]' || grep -h 'device SPA:' "$WORK/$C/gpu_text/log.txt" | grep -q '[1-9][0-9]* asked for Firth' || echo "    note: the GPU gate flagged no pair for Firth before SPA (Firth decided after SPA only)"
  fi
  grep -h '^  gate:\|^  AF_case/AF_ctrl:\|^  GPU coverage\|^  gpuSpa:\|device SPA' "$WORK/$C/gpu_text/log.txt" | sed 's/^/    /'
  # 1. sgs -> text round trip of the GPU run
  rm -rf "$WORK/$C/gpu_sgs/txt"; mkdir -p "$WORK/$C/gpu_sgs/txt"
  for f in "$WORK/$C/gpu_sgs/out/"*.txt.sgs; do
    [[ "$f" == *".markers.sgs" ]] && continue
    b=$(basename "$f" .sgs)
    "$CVT" -o "$WORK/$C/gpu_sgs/txt/$b" "$f" > /dev/null 2>&1 || { echo "    sgs2txt failed: $b"; FAIL=1; }
  done
  nd=0; nf=0
  for f in "$WORK/$C/gpu_text/out/"*.txt; do
    b=$(basename "$f"); nf=$((nf+1))
    cmp -s "$f" "$WORK/$C/gpu_sgs/txt/$b" || { echo "    DIFF sgs2txt vs text: $b"; nd=$((nd+1)); }
  done
  [ "$nd" = 0 ] && echo "    sgs -> text round trip: $nf files byte-identical" || FAIL=1
  # 2. CPU vs GPU
  python3 "$CMP" --cpu "$WORK/$C/cpu_text/out" --gpu "$WORK/$C/gpu_text/out" \
          --routes-cpu "$WORK/$C/cpu_text/routes" --routes-gpu "$WORK/$C/gpu_text/routes" \
          --sgs-cpu "$WORK/$C/cpu_sgs/out" --sgs-gpu "$WORK/$C/gpu_sgs/out" --tag "$C" \
          > "$WORK/$C/cmp.txt" 2>&1 || FAIL=1
  sed 's/^/    /' "$WORK/$C/cmp.txt"
done
echo
[ "$FAIL" = 0 ] && echo "BINGPU GATE: PASS" || echo "BINGPU GATE: FAIL"
exit $FAIL
