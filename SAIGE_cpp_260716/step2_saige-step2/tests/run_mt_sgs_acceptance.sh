#!/bin/bash
# run_mt_sgs_acceptance.sh -- outputFormat: sgs on the CPU multi-trait path
# (mainMarkerMT). The GPU path has had it since 2026-09-18; this is the same
# sink driven from the other loop, so the question is only whether that loop
# hands it the same columns.
#
# Two gates:
#
#   1. ROUND TRIP. For every case, run it twice -- once with the default text
#      writer, once with outputFormat: sgs -- and convert the .sgs back with
#      tools/sgs2txt. Every trait's text must come back BYTE-IDENTICAL.
#
#   2. TEXT UNTOUCHED. The same text runs, done with a baseline binary, must
#      produce the same bytes as the binary under test. `outputFormat` absent
#      has to mean nothing changed.
#
# The cases deliberately include everything the GPU path cannot reach, because
# that is the part of the sgs writer the GPU acceptance run never exercised:
# binary traits (p.value.NA / Is.SPA / AF_case / N_case / the isMoreOutput
# quartet), mixed binary+quantitative in one run (traits are reordered
# binary-first internally, so the marker block and the trait blocks have to stay
# in step), the scalar fallback (mtBatch: false), different sample sets per
# model (per-trait AC / AF / MissingRate, which is what the format's
# F_OVERRIDE_* flags exist for), missing calls and QC-dropped markers (the
# present mask), and rare variants with ER.
#
# usage: tests/run_mt_sgs_acceptance.sh [-n NEW_BIN] [-b BASE_BIN] [-w WORKDIR]
#                                       [-c "CASE ..."] [-1|-2]
#   -1 gate 1 only   -2 gate 2 only
# Exit status is 0 only if every comparison passed.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEW="$HERE/../saige-step2"
BASE="/opt/saige/logs/mtsgs/bin/saige-step2.base"
WORK="/opt/saige/logs/mtsgs/accept"
DATA="/opt/saige/logs/tg2_step2"
SUBNULL="${MS_NULL:-/opt/saige/logs/missing_mt/step2/null}"
MID2K="/opt/saige/data/mid2k"
RARE="/opt/saige/data/rare"
CASES="q8 q32 qnobatch qmiss b8 bmiss mixed differ_block differ_indep rare_er"
G1=1; G2=1
while getopts "n:b:w:c:12" o; do case $o in
  n) NEW=$OPTARG;; b) BASE=$OPTARG;; w) WORK=$OPTARG;; c) CASES=$OPTARG;;
  1) G2=0;; 2) G1=0;;
esac; done
NEW="$(cd "$(dirname "$NEW")" && pwd)/$(basename "$NEW")"
CVT="$HERE/../tools/sgs2txt"
[ -x "$NEW" ] || { echo "FAIL: no binary $NEW"; exit 1; }
[ -x "$CVT" ] || { echo "FAIL: no converter $CVT (make tools/sgs2txt)"; exit 1; }
MODELS="$DATA/runs/s1_cpp_P128/out"     # 128 quantitative models, n = 50,000
BINMODELS="$DATA/bin_models/spa"        # 8 binary models, same samples
mkdir -p "$WORK"
FAIL=0
echo "workdir  $WORK"
echo "new      $NEW"
echo "baseline $BASE"
echo "cases    $CASES"
echo

# ------------------------------------------------ a small file with MISSING
# Same construction as tests/run_gpu_acceptance.sh, same seed, so the two
# acceptance runs exercise the same bed.
MISS="$WORK/g5kmiss"
if [ ! -f "$MISS.bed" ]; then
  echo "building $MISS.bed (g5k with injected missing calls) ..."
  python3 - "$DATA/data/g5k" "$MISS" <<'PY'
import numpy as np, shutil, sys
src, dst = sys.argv[1], sys.argv[2]
N, M = 50000, 5000
bpv = (N + 3) // 4
raw = np.fromfile(src + '.bed', dtype=np.uint8)
assert raw[0] == 0x6c and raw[1] == 0x1b and raw[2] == 0x01
body = raw[3:].reshape(M, bpv).copy()
codes = np.zeros((M, N), dtype=np.uint8)
for k in range(4):
    codes[:, k::4] = (body >> (2 * k)) & 3
rng = np.random.default_rng(20260918)
mask = rng.random((M, N)) < 0.02
hi = rng.choice(M, 50, replace=False)
mask[hi] |= rng.random((50, N)) < 0.30
codes[mask] = 1
out = np.zeros((M, bpv), dtype=np.uint8)
for k in range(4):
    out |= (codes[:, k::4].astype(np.uint8) << (2 * k))
with open(dst + '.bed', 'wb') as f:
    f.write(bytes([0x6c, 0x1b, 0x01])); out.tofile(f)
shutil.copy(src + '.bim', dst + '.bim'); shutil.copy(src + '.fam', dst + '.fam')
mr = (codes == 1).mean(axis=1)
print(f"  overall missing {(codes==1).mean():.4f}; {int((mr>0.15).sum())} markers over maxMissRate")
PY
fi

# ---------------------------------------------------------------- config gen
# head <bed> <extra keys...>
head_yaml() {
  local BED=$1; shift
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  for L in "$@"; do echo "$L"; done
  echo "models:"
}
# one_model <name> <modeldir> <vrfile> <outdir>
one_model() {
  echo "  - traitName: $1"; echo "    modelFile: $2"
  echo "    varianceRatioFile: $3"; echo "    outputFile: $4/$1.txt"
}
# seq_models <P> <modelroot> <prefix> <outdir>
seq_models() {
  local P=$1 MR=$2 PF=$3 OD=$4 k
  for k in $(seq 1 "$P"); do
    one_model "$PF$k" "$MR/m/$PF$k" "$MR/mvr_$PF$k.varianceRatio.txt" "$OD"
  done
}

# cfg_for <case> <outdir>  -> prints the whole yaml
cfg_for() {
  local C=$1 OD=$2
  case "$C" in
    q8)       head_yaml "$DATA/data/g5k";                    seq_models 8  "$MODELS"    y "$OD" ;;
    q32)      head_yaml "$DATA/data/g5k";                    seq_models 32 "$MODELS"    y "$OD" ;;
    qnobatch) head_yaml "$DATA/data/g5k" "mtBatch: false";   seq_models 8  "$MODELS"    y "$OD" ;;
    qmiss)    head_yaml "$MISS";                             seq_models 8  "$MODELS"    y "$OD" ;;
    b8)       head_yaml "$DATA/data/g5k" "isMoreOutput: true" "isFirth: true" \
                        "is_Firth_beta: true" "pCutoffforFirth: 0.05"
              seq_models 8 "$BINMODELS" y "$OD" ;;
    bmiss)    head_yaml "$MISS" "isMoreOutput: true";        seq_models 8  "$BINMODELS" y "$OD" ;;
    mixed)    # 4 binary + 4 quantitative in one run, interleaved in the config
              head_yaml "$DATA/data/g5k" "isMoreOutput: true"
              local k
              for k in 1 2 3 4; do
                one_model "b$k" "$BINMODELS/m/y$k" "$BINMODELS/mvr_y$k.varianceRatio.txt" "$OD"
                one_model "q$k" "$MODELS/m/y$k"    "$MODELS/mvr_y$k.varianceRatio.txt"    "$OD"
              done ;;
    differ_block)  # two sample sets: y1-y8 the union, y9-y16 a 60% subset
              head_yaml "$MID2K"
              for k in 1 2 3 4 9 10 11 12; do
                one_model "y$k" "$SUBNULL/block16/m/y$k" \
                          "$SUBNULL/block16/mvr_y$k.varianceRatio.txt" "$OD"
              done ;;
    differ_indep)  # 8 different sample sets, each missing an independent 5%
              head_yaml "$MID2K" "isMoreOutput: true"
              for k in 1 2 3 4 5 6 7 8; do
                one_model "y$k" "$SUBNULL/indep16/m/y$k" \
                          "$SUBNULL/indep16/mvr_y$k.varianceRatio.txt" "$OD"
              done ;;
    rare_er)  # rare variants: ER at MAC <= 4, the dosage-zeroing gate, missing
              head_yaml "$RARE" "isMoreOutput: true" "MACCutoffforER: 20"
              local k
              for k in 1 2 3 4; do
                one_model "b$k" "$SUBNULL/indep16/m/y$k" \
                          "$SUBNULL/indep16/mvr_y$k.varianceRatio.txt" "$OD"
                one_model "q$k" "$SUBNULL/qmiss/m/qm$k" \
                          "$SUBNULL/qmiss/mvr_qm$k.varianceRatio.txt" "$OD"
              done ;;
    *) echo "unknown case $C" >&2; return 1 ;;
  esac
}

have_case() {  # the differ / rare cases need assets that may not be installed
  case "$1" in
    differ_*|rare_er) [ -f "$SUBNULL/indep16/m/y1/nullmodel.json" ] && [ -f "$MID2K.bed" ] ;;
    *) return 0 ;;
  esac
}

run() {  # run <bin> <cfg> <log>
  "$1" "$2" > "$3" 2>&1
  local rc=$?
  [ $rc -ne 0 ] && { echo "    RUN FAILED rc=$rc, see $3"; tail -8 "$3"; }
  return $rc
}

# ==========================================================================
# GATE 1 -- sgs -> sgs2txt is byte-identical to the text run
# ==========================================================================
if [ $G1 -eq 1 ]; then
echo "=== gate 1: sgs round trip (CPU multi-trait path) ==="
for C in $CASES; do
  have_case "$C" || { printf "  %-14s SKIP (assets missing)\n" "$C"; continue; }
  TD="$WORK/g1_${C}_txt"; SD="$WORK/g1_${C}_sgs"
  rm -rf "$TD" "$SD"; mkdir -p "$TD" "$SD"
  cfg_for "$C" "$TD"                     > "$WORK/g1_$C.txt.yaml" || { FAIL=1; continue; }
  { cfg_for "$C" "$SD"; echo "outputFormat: sgs"; } > "$WORK/g1_$C.sgs.yaml"
  run "$NEW" "$WORK/g1_$C.txt.yaml" "$WORK/g1_$C.txt.log" || { FAIL=1; continue; }
  run "$NEW" "$WORK/g1_$C.sgs.yaml" "$WORK/g1_$C.sgs.log" || { FAIL=1; continue; }
  grep -q "useGPU" "$WORK/g1_$C.sgs.log" && echo "    NOTE: a GPU line appeared in $C"
  ls "$SD"/*.sgs >/dev/null 2>&1 || { echo "    $C: no .sgs written"; FAIL=1; continue; }
  "$CVT" -j 8 "$SD"/*.sgs > "$WORK/g1_$C.cvt.log" 2>&1 || {
      echo "    $C: sgs2txt failed"; tail -5 "$WORK/g1_$C.cvt.log"; FAIL=1; continue; }
  n=0; bad=0; rows=0
  for f in "$TD"/*.txt; do
    g="$SD/$(basename "$f")"; n=$((n+1))
    rows=$((rows + $(wc -l < "$f") - 1))
    cmp -s "$f" "$g" || { bad=$((bad+1)); echo "    DIFFERS: $(basename "$f")"; }
  done
  tb=$(du -sb "$TD" | cut -f1)
  sb=$(du -sb --exclude='*.txt' "$SD" | cut -f1)
  printf "  %-14s %2d traits, %8d rows, %d differ   text %sB -> sgs %sB (%.2fx)\n" \
         "$C" "$n" "$rows" "$bad" "$tb" "$sb" \
         "$(awk -v a=$tb -v b=$sb 'BEGIN{print (b>0? a/b : 0)}')"
  [ "$bad" -ne 0 ] && FAIL=1
done
echo
fi

# ==========================================================================
# GATE 2 -- outputFormat absent: byte-identical to the baseline binary
# ==========================================================================
if [ $G2 -eq 1 ]; then
if [ -x "$BASE" ]; then
echo "=== gate 2: text output unchanged against $BASE ==="
NCMP=0; NBAD=0
for C in $CASES; do
  have_case "$C" || { printf "  %-14s SKIP (assets missing)\n" "$C"; continue; }
  OB="$WORK/g2_${C}_base"; ON="$WORK/g2_${C}_new"
  rm -rf "$OB" "$ON"; mkdir -p "$OB" "$ON"
  cfg_for "$C" "$OB" > "$WORK/g2_$C.base.yaml"
  cfg_for "$C" "$ON" > "$WORK/g2_$C.new.yaml"
  run "$BASE" "$WORK/g2_$C.base.yaml" "$WORK/g2_$C.base.log" || { FAIL=1; continue; }
  run "$NEW"  "$WORK/g2_$C.new.yaml"  "$WORK/g2_$C.new.log"  || { FAIL=1; continue; }
  if diff -r -q "$OB" "$ON" > "$WORK/g2_$C.diff" 2>&1; then
    d=0
  else
    d=$(wc -l < "$WORK/g2_$C.diff"); cat "$WORK/g2_$C.diff" | sed 's/^/    /'
  fi
  n=$(ls "$OB" | wc -l)
  NCMP=$((NCMP+n)); NBAD=$((NBAD+d))
  printf "  %-14s %2d files, diff -r: %s\n" "$C" "$n" \
         "$([ "$d" -eq 0 ] && echo clean || echo "$d differences")"
  [ "$d" -ne 0 ] && FAIL=1
done
echo "  ---- gate 2: $NCMP files compared, $NBAD differences ----"
echo
else
echo "=== gate 2: skipped, no baseline binary at $BASE ==="; echo
fi
fi

# ==========================================================================
# GATE 3 -- the paths that do NOT write sgs must say so, not write nothing
# ==========================================================================
echo "=== refusals ==="
RD="$WORK/refuse"; rm -rf "$RD"; mkdir -p "$RD"
{ head_yaml "$DATA/data/g5k" "outputFormat: sgs"
  one_model y1 "$MODELS/m/y1" "$MODELS/mvr_y1.varianceRatio.txt" "$RD"; } > "$RD/p1.yaml"
if "$NEW" "$RD/p1.yaml" > "$RD/p1.log" 2>&1; then
  echo "  P=1, useGPU off: ACCEPTED -- should have refused"; FAIL=1
else
  grep -q "single-trait path" "$RD/p1.log" \
    && echo "  P=1, useGPU off: refused, message names the single-trait path" \
    || { echo "  P=1, useGPU off: refused with the wrong message:"; \
         grep -i "sgs" "$RD/p1.log" | head -3; FAIL=1; }
fi
GRP="$WORK/g5k.group"
[ -f "$GRP" ] || python3 "$HERE/make_group_file.py" "$DATA/data/g5k.bim" "$GRP" all >/dev/null
{ head_yaml "$DATA/data/g5k" "outputFormat: sgs" "groupFile: $GRP" \
      "annotationList:" "  - \"null\"" "maxMAFList:" "  - 0.5" "r_corr: 0"
  one_model y1 "$MODELS/m/y1" "$MODELS/mvr_y1.varianceRatio.txt" "$RD"; } > "$RD/region.yaml"
if "$NEW" "$RD/region.yaml" > "$RD/region.log" 2>&1; then
  echo "  region test:      ACCEPTED -- should have refused"; FAIL=1
else
  grep -q "region / group testing" "$RD/region.log" \
    && echo "  region test:      refused, message names region testing" \
    || { echo "  region test:      refused with the wrong message:"; \
         grep -i "sgs" "$RD/region.log" | head -3; FAIL=1; }
fi
# A one-model run with useGPU: true goes through mainMarkerMT (and falls back to
# it when no device is present), so it MUST be accepted.
SD1="$WORK/refuse/p1gpu"; mkdir -p "$SD1"
{ head_yaml "$DATA/data/g5k" "outputFormat: sgs" "useGPU: true"
  one_model y1 "$MODELS/m/y1" "$MODELS/mvr_y1.varianceRatio.txt" "$SD1"; } > "$RD/p1gpu.yaml"
if "$NEW" "$RD/p1gpu.yaml" > "$RD/p1gpu.log" 2>&1 && [ -f "$SD1/y1.txt.sgs" ]; then
  echo "  P=1, useGPU on:   accepted, wrote $SD1/y1.txt.sgs"
else
  echo "  P=1, useGPU on:   FAILED -- should have written sgs"; tail -5 "$RD/p1gpu.log"; FAIL=1
fi
echo
[ $FAIL -eq 0 ] && echo "ALL CHECKS PASSED" || echo "SOME CHECKS FAILED"
exit $FAIL
