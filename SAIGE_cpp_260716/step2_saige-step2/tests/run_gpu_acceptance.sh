#!/bin/bash
# run_gpu_acceptance.sh — the repeatable check for the step-2 GPU path
# (config key useGPU; gpu/gpu_step2.hpp).
#
# Two questions, and they need different answers:
#
#   A. FLAG OFF: byte identity. With useGPU absent the binary must produce the
#      same bytes as a pre-GPU build, for every trait type and every test type
#      -- including the ones the GPU path would never touch (binary traits with
#      SPA, region/set tests). Anything else means the plumbing leaked.
#
#   B. FLAG ON: agreement within a stated tolerance. The CPU path reduces in
#      double with OpenBLAS's dgemm; the GPU reduces with cuBLAS in double
#      (gpuPrecision: fp64, the default) or float (fp32). Neither association
#      order is the other's, so byte identity is not available even in fp64 --
#      but the EXACT quantities must still be exact:
#
#        same rows, same order                      (the QC decisions)
#        CHR POS MarkerID Allele1 Allele2           (the marker)
#        AC_Allele2 AF_Allele2 MissingRate N        (counts and imputation)
#
#      and only BETA / SE / Tstat / var / p.value may move. Those are reported
#      as max and median |delta(-log10 p)| and |delta Tstat| / sd(Tstat).
#
# Usage:
#   tests/run_gpu_acceptance.sh [-b BASELINE_BIN] [-n NEW_BIN] [-d DATADIR]
#                               [-w WORKDIR] [-q] [-A|-B]
#     -b  pre-GPU binary for part A (default: tests/saige-step2.pre-gpu if present)
#     -n  binary under test        (default: ./saige-step2)
#     -q  quick: 5,000-marker file only; skip the 10^6-marker tolerance run
#     -A  part A only     -B  part B only
#
# Exit status is 0 only if every comparison passed.
set -u

BASE=""
NEW="./saige-step2"
DATA="/opt/saige/logs/tg2_step2"
WORK="/opt/saige/logs/gpu_step2/accept"
QUICK=0
DOA=1; DOB=1
while getopts "b:n:d:w:qAB" o; do case $o in
  b) BASE=$OPTARG;; n) NEW=$OPTARG;; d) DATA=$OPTARG;; w) WORK=$OPTARG;;
  q) QUICK=1;; A) DOB=0;; B) DOA=0;;
esac; done

HERE="$(cd "$(dirname "$0")" && pwd)"
[ -z "$BASE" ] && [ -f "$HERE/saige-step2.pre-gpu" ] && BASE="$HERE/saige-step2.pre-gpu"
NEW="$(cd "$(dirname "$NEW")" && pwd)/$(basename "$NEW")"

MODELS="$DATA/runs/s1_cpp_P128/out"          # 128 quantitative null models
BINMODELS="$DATA/bin_models/spa"             # 8 binary null models
mkdir -p "$WORK"
FAIL=0
echo "workdir  $WORK"
echo "new      $NEW"
echo "baseline ${BASE:-<none: part A skipped>}"
echo

# ---------------------------------------------------------------- config gen
# $1 P, $2 outdir, $3 plink prefix, $4 model root, rest: extra top-level keys
gen() {
  local P=$1 OD=$2 BED=$3 MR=$4; shift 4
  mkdir -p "$OD"
  echo "genoType: plink"
  echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  for k in $(seq 1 "$P"); do
    echo "  - traitName: y$k"
    echo "    modelFile: $MR/m/y$k"
    echo "    varianceRatioFile: $MR/mvr_y$k.varianceRatio.txt"
    echo "    outputFile: $OD/y$k.txt"
  done
}

run() {  # run <bin> <cfg> <log>
  "$1" "$2" > "$3" 2>&1
  local rc=$?
  [ $rc -ne 0 ] && { echo "    RUN FAILED rc=$rc, see $3"; tail -5 "$3"; }
  return $rc
}

# ------------------------------------------------ a small file with MISSING
# The benchmark .bed has no missing calls at all, and missing genotypes are the
# one thing the GPU path handles differently in kind (a per-marker entry in the
# code->dosage table instead of a per-cell branch). So part A and part B both
# run against a copy with ~2% missing plus 50 markers over the maxMissRate
# cutoff, to exercise both the imputation and the QC drop.
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
mask[hi] |= rng.random((50, N)) < 0.30      # over maxMissRate: must be dropped
codes[mask] = 1                             # PLINK code 01 == missing
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

# ==========================================================================
# PART A — flag off must be byte-identical to a pre-GPU build
# ==========================================================================
if [ $DOA -eq 1 ] && [ -n "$BASE" ]; then
  echo "=== A. flag off: byte identity against $BASE ==="
  NCMP=0; NBAD=0
  caseA() {  # caseA <name> <cfg-file>
    local name=$1 cfg=$2
    local ob="$WORK/A_${name}_base" on="$WORK/A_${name}_new"
    rm -rf "$ob" "$on"; mkdir -p "$ob" "$on"
    sed "s|$WORK/A_${name}_out|$ob|" "$cfg" > "$cfg.base"
    sed "s|$WORK/A_${name}_out|$on|" "$cfg" > "$cfg.new"
    run "$BASE" "$cfg.base" "$WORK/A_${name}_base.log" || { FAIL=1; return; }
    run "$NEW"  "$cfg.new"  "$WORK/A_${name}_new.log"  || { FAIL=1; return; }
    local n=0 bad=0
    for f in "$ob"/*; do
      local g="$on/$(basename "$f")"
      n=$((n+1)); NCMP=$((NCMP+1))
      cmp -s "$f" "$g" || { bad=$((bad+1)); NBAD=$((NBAD+1)); echo "    DIFFERS: $(basename "$f")"; }
    done
    printf "  %-28s %2d files, %d differ\n" "$name" "$n" "$bad"
    [ $bad -ne 0 ] && FAIL=1
  }

  gen 1  "$WORK/A_q1_out"    "$DATA/data/g5k" "$MODELS"    > "$WORK/A_q1.yaml"
  gen 8  "$WORK/A_q8_out"    "$DATA/data/g5k" "$MODELS"    > "$WORK/A_q8.yaml"
  gen 32 "$WORK/A_q32_out"   "$DATA/data/g5k" "$MODELS"    > "$WORK/A_q32.yaml"
  gen 8  "$WORK/A_qmiss_out" "$MISS"          "$MODELS"    > "$WORK/A_qmiss.yaml"
  gen 8  "$WORK/A_qnobatch_out" "$DATA/data/g5k" "$MODELS" "mtBatch: false" > "$WORK/A_qnobatch.yaml"
  # Binary traits: SPA on, Firth on, more output. None of this is reachable
  # from the GPU path -- which is exactly why it is worth comparing.
  gen 1  "$WORK/A_b1_out" "$DATA/data/g5k" "$BINMODELS" "isMoreOutput: true" "isFirth: true" "is_Firth_beta: true" "pCutoffforFirth: 0.05" > "$WORK/A_b1.yaml"
  gen 8  "$WORK/A_b8_out" "$DATA/data/g5k" "$BINMODELS" "isMoreOutput: true" "isFirth: true" "is_Firth_beta: true" "pCutoffforFirth: 0.05" > "$WORK/A_b8.yaml"
  gen 8  "$WORK/A_bmiss_out" "$MISS"       "$BINMODELS" "isMoreOutput: true" > "$WORK/A_bmiss.yaml"

  for c in q1 q8 q32 qmiss qnobatch b1 b8 bmiss; do caseA "$c" "$WORK/A_$c.yaml"; done

  # Region / set-based test: single trait, group file, SKAT-O + burden.
  # Never GPU, always CPU -- which is exactly why it is worth comparing.
  GRP="$WORK/g5k.group"
  [ -f "$GRP" ] || python3 "$HERE/make_group_file.py" "$DATA/data/g5k.bim" "$GRP" all
  gen 1 "$WORK/A_region_out" "$DATA/data/g5k" "$MODELS" \
      "groupFile: $GRP" \
      "annotationList:" "  - \"null\"" \
      "maxMAFList:" "  - 0.01" "  - 0.5" \
      "r_corr: 0" "MACCutoff_to_CollapseUltraRare: 10" \
      "is_single_in_groupTest: true" \
      "markers_per_chunk_in_groupTest: 500" > "$WORK/A_region.yaml"
  caseA region "$WORK/A_region.yaml"

  echo "  ---- part A: $NCMP files compared, $NBAD differ ----"
  echo
elif [ $DOA -eq 1 ]; then
  echo "=== A. skipped: no baseline binary (-b) ==="; echo
fi

# ==========================================================================
# PART B — flag on: same rows, exact QC columns, statistic within tolerance
# ==========================================================================
if [ $DOB -eq 1 ]; then
  echo "=== B. flag on: CPU vs GPU ==="
  # The CPU side of a (bed, P) pair is run ONCE and reused by both precisions;
  # at 10^6 markers and P=128 that run is 25 minutes and there is no reason to
  # pay for it twice.
  caseB() {  # caseB <name> <P> <bed> <prec>
    local name=$1 P=$2 BED=$3 PREC=$4
    local key; key="$(basename "$BED")_P${P}"
    local oc="$WORK/B_cpu_$key" og="$WORK/B_${name}_gpu"
    if [ ! -f "$oc/y${P}.txt" ] || [ "$NEW" -nt "$oc/y${P}.txt" ]; then
      rm -rf "$oc"
      gen "$P" "$oc" "$BED" "$MODELS" > "$WORK/B_cpu_$key.yaml"
      run "$NEW" "$WORK/B_cpu_$key.yaml" "$WORK/B_cpu_$key.log" || { FAIL=1; return; }
    fi
    rm -rf "$og"
    gen "$P" "$og" "$BED" "$MODELS" "useGPU: true" "gpuPrecision: $PREC" > "$WORK/B_$name.gpu.yaml"
    run "$NEW" "$WORK/B_$name.gpu.yaml" "$WORK/B_$name.gpu.log" || { FAIL=1; return; }
    grep -q "useGPU: refused" "$WORK/B_$name.gpu.log" && {
        echo "  $name: the GPU path REFUSED the run:"; grep "useGPU: refused" "$WORK/B_$name.gpu.log"; FAIL=1; return; }
    python3 "$HERE/cmp_gpu.py" --cpu "$oc" --gpu "$og" --P "$P" --tag "$name/$PREC" || FAIL=1
  }

  caseB q8_fp64     8 "$DATA/data/g5k" fp64
  caseB q8_fp32     8 "$DATA/data/g5k" fp32
  caseB qmiss_fp64  8 "$MISS"          fp64
  caseB qmiss_fp32  8 "$MISS"          fp32
  if [ $QUICK -eq 0 ]; then
    # 10^6 markers: the tolerance numbers that get reported come from here.
    for P in 1 8 32 128; do
      caseB "M1e6_P${P}_fp64" "$P" "$DATA/data/g1m" fp64
      caseB "M1e6_P${P}_fp32" "$P" "$DATA/data/g1m" fp32
    done
  fi
  echo
fi

if [ $FAIL -eq 0 ]; then echo "ACCEPTANCE: PASS"; else echo "ACCEPTANCE: FAIL"; fi
exit $FAIL
