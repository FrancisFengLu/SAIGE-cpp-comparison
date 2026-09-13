#!/bin/bash
# P=1 byte-identity gate for the multi-trait work (MULTITRAIT_DESIGN.md section 5).
#
# Every multi-trait commit must leave single-trait output BYTE-IDENTICAL. This
# runs the same config matrix through a baseline binary and the binary under
# test and cmp's every output file. "Relative error is tiny" is not a pass.
#
# usage: run_p1_regression.sh <new binary> <baseline binary> [workdir]
#
#   the baseline binary is normally built from the commit before the change:
#     git worktree add /tmp/base <sha> && (cd /tmp/base/.../step2_saige-step2 && make -j8)
#
# env overrides:
#   P1_PLINK_BIG   default /opt/saige/data/mid      (50k samples x 40k markers)
#   P1_PLINK_SMALL default /opt/saige/data/mid2k    (same samples, 2k markers)
#   P1_MODEL       default /opt/saige/logs/step2/null_y4_arma
#   P1_VR          default /opt/saige/logs/mp/single_y4.varianceRatio.txt
#   P1_GROUP       default /opt/saige/logs/step2/mid_group_400.txt
#
# Coverage, one config per code path that can reach the output writer:
#   R1 plain scalar-key run                      R6 same run via `models:` (P=1)
#   R2 isMoreOutput + Firth                      R7 per-model override beats top level
#   R3 blockSize=8 (scoreTestFast_block)         R8 LOCO, reads chr<N>/
#   R4 isnoadjCov (scoreTestFast_noadjCov)       R9 LOCO silent fallback
#   R5 region / group test
#
# R3/R4 use the small plink set on purpose: blockSize>1 prefetches every
# marker's N-vector, so 40k markers x 50k samples exhausts memory (a
# pre-existing issue, design section 7.1).

set -u

NEW="${1:?usage: run_p1_regression.sh <new binary> <baseline binary> [workdir]}"
BASE="${2:?usage: run_p1_regression.sh <new binary> <baseline binary> [workdir]}"
WORK="${3:-$(mktemp -d)}"

BIG="${P1_PLINK_BIG:-/opt/saige/data/mid}"
SMALL="${P1_PLINK_SMALL:-/opt/saige/data/mid2k}"
MODEL="${P1_MODEL:-/opt/saige/logs/step2/null_y4_arma}"
VR="${P1_VR:-/opt/saige/logs/mp/single_y4.varianceRatio.txt}"
GROUP="${P1_GROUP:-/opt/saige/logs/step2/mid_group_400.txt}"
MAKE_LOCO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/make_loco_model.py"

for f in "$MODEL/nullmodel.json" "$VR" "$BIG.bed" "$SMALL.bed"; do
    [ -e "$f" ] || { echo "SKIP: test data not found: $f"; exit 0; }
done
[ -x "$NEW" ]  || { echo "FAIL: binary not found: $NEW";  exit 1; }
[ -x "$BASE" ] || { echo "FAIL: baseline not found: $BASE"; exit 1; }

CFG="$WORK/cfg"; mkdir -p "$CFG" "$WORK/base" "$WORK/new"

# ---- synthetic LOCO models -------------------------------------------------
# loco_a has chr1/ (so chrom=1 really swaps files in); loco_b lists only chr2,
# so chrom=1 hits the loader's silent fallback to the genome-wide fit. The two
# must produce DIFFERENT numbers, otherwise the LOCO configs prove nothing.
mk_loco () {  # mk_loco <dst> <chrom> <factor> <loco_chroms json>
    python3 - "$MODEL" "$1" "$2" "$3" "$4" "$MAKE_LOCO" <<'PY'
import importlib.util, json, os, shutil, sys
src, dst, chrom, factor, chroms, helper = sys.argv[1:7]
spec = importlib.util.spec_from_file_location("mlm", helper)
mlm = importlib.util.module_from_spec(spec); spec.loader.exec_module(mlm)
if os.path.exists(dst): shutil.rmtree(dst)
shutil.copytree(src, dst)
cdir = os.path.join(dst, "chr" + chrom); os.makedirs(cdir, exist_ok=True)
for name in mlm.PER_CHROM:
    header, rows, cols, vals = mlm.read_arma(os.path.join(src, name + ".arma"))
    mlm.write_arma(os.path.join(cdir, name + ".arma"), header, rows, cols,
                   mlm.perturb(name, vals, float(factor)))
j = json.load(open(os.path.join(dst, "nullmodel.json")))
j["loco"] = True; j["loco_chroms"] = json.loads(chroms)
json.dump(j, open(os.path.join(dst, "nullmodel.json"), "w"))
PY
}
mk_loco "$WORK/loco_a" 1 0.90 '[1]' || { echo "FAIL: could not build LOCO model"; exit 1; }
mk_loco "$WORK/loco_b" 2 1.15 '[2]' || { echo "FAIL: could not build LOCO model"; exit 1; }

common () {  # common <plinkFile>
    cat <<EOF
genoType:   plink
plinkFile:  $1
minMAF: 0
minMAC: 1
maxMissRate: 0.15
AlleleOrder: alt-first
MACCutoffforER: 4
nThreads: 1
EOF
}

cat > "$CFG/R1_scalar.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile: __OUT__/R1_scalar.txt
isMoreOutput: false
isFirth: false
$(common "$BIG")
EOF
cat > "$CFG/R2_more_firth.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile: __OUT__/R2_more_firth.txt
isMoreOutput: true
isFirth: true
is_Firth_beta: true
pCutoffforFirth: 0.05
$(common "$BIG")
EOF
cat > "$CFG/R3_block8.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile: __OUT__/R3_block8.txt
isMoreOutput: false
isFirth: false
blockSize: 8
$(common "$SMALL")
EOF
cat > "$CFG/R4_noadjcov.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile: __OUT__/R4_noadjcov.txt
isnoadjCov: true
isMoreOutput: false
isFirth: false
$(common "$SMALL")
EOF
cat > "$CFG/R8_loco.yaml" <<EOF
modelFile:         $WORK/loco_a
varianceRatioFile: $VR
outputFile: __OUT__/R8_loco.txt
LOCO: true
chrom: "1"
isMoreOutput: true
isFirth: false
$(common "$SMALL")
EOF
cat > "$CFG/R9_loco_fallback.yaml" <<EOF
modelFile:         $WORK/loco_b
varianceRatioFile: $VR
outputFile: __OUT__/R9_loco_fallback.txt
LOCO: true
chrom: "1"
isMoreOutput: true
isFirth: false
$(common "$SMALL")
EOF
if [ -e "$GROUP" ]; then
cat > "$CFG/R5_region.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile: __OUT__/R5_region.txt
groupFile:  $GROUP
annotationList:
  - "lof"
maxMAFList:
  - 0.001
  - 0.01
r_corr: 0
MACCutoff_to_CollapseUltraRare: 10
markers_per_chunk_in_groupTest: 500
groups_per_chunk: 100
isFirth: false
$(common "$BIG")
EOF
fi

# `models:` form -- new binary only; compared against the LEGACY baselines,
# which is the whole point (design gates G0.2 / G0.3).
cat > "$CFG/R6_models_p1.yaml" <<EOF
models:
  - traitName: t0
    modelFile:         $MODEL
    varianceRatioFile: $VR
    outputFile: __OUT__/R6_models_p1.txt
isMoreOutput: false
isFirth: false
$(common "$BIG")
EOF
cat > "$CFG/R7_models_p1_ov.yaml" <<EOF
isnoadjCov: false
models:
  - traitName: t0
    modelFile:         $MODEL
    varianceRatioFile: $VR
    outputFile: __OUT__/R7_models_p1_ov.txt
    isnoadjCov: true
isMoreOutput: false
isFirth: false
$(common "$SMALL")
EOF

run_one () {  # run_one <binary> <outdir> <config name>
    local bin="$1" out="$2" n="$3"
    sed "s#__OUT__#$out#g" "$CFG/$n.yaml" > "$out/$n.yaml"
    "$bin" "$out/$n.yaml" > "$out/$n.log" 2>&1
    echo "    $n exit=$?"
}

LEGACY="R1_scalar R2_more_firth R3_block8 R4_noadjcov R8_loco R9_loco_fallback"
[ -e "$CFG/R5_region.yaml" ] && LEGACY="$LEGACY R5_region"
MODELSFORM="R6_models_p1 R7_models_p1_ov"

echo "=== baseline: $BASE ==="
for n in $LEGACY; do run_one "$BASE" "$WORK/base" "$n"; done
echo "=== under test: $NEW ==="
for n in $LEGACY $MODELSFORM; do run_one "$NEW" "$WORK/new" "$n"; done

pass=0; fail=0
chk () {  # chk <new file> <baseline file>
    if [ ! -s "$WORK/new/$1" ]; then echo "  FAIL $1: missing or empty"; fail=$((fail+1)); return; fi
    if cmp -s "$WORK/new/$1" "$WORK/base/$2"; then
        echo "  PASS $1 == base/$2  $(md5sum < "$WORK/new/$1" | cut -d' ' -f1)"; pass=$((pass+1))
    else
        echo "  FAIL $1 != base/$2"; fail=$((fail+1))
    fi
}

echo
echo "=== byte comparison ==="
for n in $LEGACY; do chk "$n.txt" "$n.txt"; done
[ -e "$WORK/base/R5_region.txt.singleAssoc.txt" ] && \
    chk R5_region.txt.singleAssoc.txt R5_region.txt.singleAssoc.txt
chk R6_models_p1.txt    R1_scalar.txt
chk R7_models_p1_ov.txt R4_noadjcov.txt

# The two LOCO configs must disagree, or they are both silently reading the
# top-level fit and prove nothing about the LOCO branch.
if cmp -s "$WORK/base/R8_loco.txt" "$WORK/base/R9_loco_fallback.txt"; then
    echo "  FAIL LOCO configs produced identical output -- the LOCO branch is not being exercised"
    fail=$((fail+1))
else
    echo "  PASS LOCO chr1 fit != LOCO genome-wide fallback"; pass=$((pass+1))
fi

echo
echo "=== $pass passed, $fail failed  (workdir $WORK) ==="
[ $fail -eq 0 ]
