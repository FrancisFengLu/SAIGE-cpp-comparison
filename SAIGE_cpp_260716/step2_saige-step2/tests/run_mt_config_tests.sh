#!/bin/bash
# Multi-trait config front end: acceptance tests (MULTITRAIT_DESIGN.md Phase 0).
#
# Two things are checked:
#   1. G0.2/G0.3 — a one-entry `models:` run is BYTE-IDENTICAL to the same run
#      written with the legacy scalar modelFile/varianceRatioFile/outputFile
#      keys, and a per-model override beats a conflicting top-level one.
#   2. Every malformed multi-trait config is rejected with a message that names
#      the offending model and key. None of these may run and produce output.
#
# usage: run_mt_config_tests.sh <saige-step2 binary> [workdir]
# env:   MT_PLINK  (default /opt/saige/data/mid2k)
#        MT_MODEL  (default /opt/saige/logs/step2/null_y4_arma)
#        MT_VR     (default /opt/saige/logs/mp/single_y4.varianceRatio.txt)

set -u

BIN="${1:?usage: run_mt_config_tests.sh <saige-step2 binary> [workdir]}"
WORK="${2:-$(mktemp -d)}"
PLINK="${MT_PLINK:-/opt/saige/data/mid2k}"
MODEL="${MT_MODEL:-/opt/saige/logs/step2/null_y4_arma}"
VR="${MT_VR:-/opt/saige/logs/mp/single_y4.varianceRatio.txt}"

mkdir -p "$WORK"
pass=0; fail=0

ok ()   { echo "  PASS  $1"; pass=$((pass+1)); }
bad ()  { echo "  FAIL  $1"; fail=$((fail+1)); }

# expect_error <name> <substring> <<< config on stdin
expect_error () {
    local name="$1" want="$2" cfg="$WORK/$1.yaml" log="$WORK/$1.log"
    cat > "$cfg"
    "$BIN" "$cfg" > "$log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        bad "$name: exited 0, expected a hard error"
    elif grep -qF "$want" "$log"; then
        ok "$name: rejected -- $(grep -m1 -oF "$want" "$log")"
    else
        bad "$name: exited $rc but the message did not mention '$want'"
        sed -n '$p' "$log" | sed 's/^/        /'
    fi
}

common () {
    cat <<EOF
genoType: plink
plinkFile: $PLINK
minMAF: 0
minMAC: 1
maxMissRate: 0.15
AlleleOrder: alt-first
isMoreOutput: false
isFirth: false
MACCutoffforER: 4
nThreads: 1
EOF
}

echo "=== multi-trait config front end ==="
echo "binary: $BIN"
echo "workdir: $WORK"
echo

# ---------------------------------------------------------------
echo "-- P=1 equivalence (byte-identical output) --"

cat > "$WORK/p1_scalar.yaml" <<EOF
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile:        $WORK/p1_scalar.txt
$(common)
EOF
cat > "$WORK/p1_models.yaml" <<EOF
models:
  - traitName: t0
    modelFile:         $MODEL
    varianceRatioFile: $VR
    outputFile:        $WORK/p1_models.txt
$(common)
EOF
# top level says false, the model says true -> the model must win
cat > "$WORK/p1_override.yaml" <<EOF
isnoadjCov: false
models:
  - traitName: t0
    modelFile:         $MODEL
    varianceRatioFile: $VR
    outputFile:        $WORK/p1_override.txt
    isnoadjCov: true
$(common)
EOF
cat > "$WORK/p1_toplevel.yaml" <<EOF
isnoadjCov: true
modelFile:         $MODEL
varianceRatioFile: $VR
outputFile:        $WORK/p1_toplevel.txt
$(common)
EOF

for n in p1_scalar p1_models p1_override p1_toplevel; do
    "$BIN" "$WORK/$n.yaml" > "$WORK/$n.log" 2>&1 || bad "$n: run failed (see $WORK/$n.log)"
done

if [ -s "$WORK/p1_scalar.txt" ] && cmp -s "$WORK/p1_scalar.txt" "$WORK/p1_models.txt"; then
    ok "models[1 entry] == scalar keys  (md5 $(md5sum < "$WORK/p1_scalar.txt" | cut -d' ' -f1))"
else
    bad "models[1 entry] != scalar keys"
fi
if [ -s "$WORK/p1_toplevel.txt" ] && cmp -s "$WORK/p1_toplevel.txt" "$WORK/p1_override.txt"; then
    ok "per-model isnoadjCov beats top-level  (md5 $(md5sum < "$WORK/p1_override.txt" | cut -d' ' -f1))"
else
    bad "per-model override did not reproduce the top-level-only run"
fi
if cmp -s "$WORK/p1_scalar.txt" "$WORK/p1_override.txt"; then
    bad "isnoadjCov override had no effect (outputs identical -- the override is being ignored)"
else
    ok "isnoadjCov override actually changed the numbers"
fi

echo
echo "-- rejection of malformed multi-trait configs --"

expect_error both_forms "mutually exclusive" <<EOF
modelFile: $MODEL
models:
  - modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/x.txt
$(common)
EOF

expect_error empty_models "is empty" <<EOF
models: []
$(common)
EOF

expect_error models_not_sequence "must be a sequence" <<EOF
models:
  modelFile: $MODEL
$(common)
EOF

expect_error missing_vr "models[0] missing required key: varianceRatioFile" <<EOF
models:
  - traitName: a
    modelFile: $MODEL
    outputFile: $WORK/x.txt
$(common)
EOF

expect_error missing_out "models[0] missing required key: outputFile" <<EOF
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
$(common)
EOF

expect_error duplicate_output "reuses outputFile" <<EOF
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/dup.txt
  - traitName: b
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/dup.txt
$(common)
EOF

# P > 1 runs now. What must still fail loudly is a model set that cannot
# legitimately share one genotype stream (design section 4.3) -- silently
# analysing a permuted sample vector is the failure mode this prevents -- and
# region testing with P > 1 (design section 10).
PERM="$WORK/perm_model"
IMPM="$WORK/imp_model"
python3 - "$MODEL" "$PERM" "$IMPM" <<'PY2'
import json, os, shutil, sys
src, perm, imp = sys.argv[1:4]
for dst in (perm, imp):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
j = json.load(open(os.path.join(perm, "nullmodel.json")))
ids = list(j["sampleIDs"])
ids[0], ids[1] = ids[1], ids[0]          # same set, different order
j["sampleIDs"] = ids
json.dump(j, open(os.path.join(perm, "nullmodel.json"), "w"))
j = json.load(open(os.path.join(imp, "nullmodel.json")))
j["impute_method"] = "bestguess"
json.dump(j, open(os.path.join(imp, "nullmodel.json"), "w"))
PY2

# Different sample lists are legal since design 4.7; mtRequireSameSamples: true
# turns them back into a hard error.
expect_error sample_order_differs "identical sample IDs in identical order" <<EOF
mtRequireSameSamples: true
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/so_a.txt
  - traitName: b
    modelFile: $PERM
    varianceRatioFile: $VR
    outputFile: $WORK/so_b.txt
$(common)
EOF

expect_error impute_method_differs "the models must agree on it" <<EOF
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/im_a.txt
  - traitName: b
    modelFile: $IMPM
    varianceRatioFile: $VR
    outputFile: $WORK/im_b.txt
$(common)
EOF

# ...and without the key the same pair of models runs, each trait on its own
# sample list, and each output is byte-identical to that model run alone.
cat > "$WORK/so_run.yaml" <<EOF
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/so_run_a.txt
  - traitName: b
    modelFile: $PERM
    varianceRatioFile: $VR
    outputFile: $WORK/so_run_b.txt
$(common)
EOF
cat > "$WORK/so_single_b.yaml" <<EOF
modelFile:         $PERM
varianceRatioFile: $VR
outputFile:        $WORK/so_single_b.txt
$(common)
EOF
"$BIN" "$WORK/so_run.yaml" > "$WORK/so_run.log" 2>&1
so_rc=$?
"$BIN" "$WORK/so_single_b.yaml" > "$WORK/so_single_b.log" 2>&1
if [ $so_rc -ne 0 ]; then
    bad "different sample order without mtRequireSameSamples: exited $so_rc ($(tail -1 "$WORK/so_run.log"))"
elif cmp -s "$WORK/so_run_a.txt" "$WORK/p1_scalar.txt" && cmp -s "$WORK/so_run_b.txt" "$WORK/so_single_b.txt"; then
    ok "different sample order runs; both traits == their single-trait runs"
else
    bad "different sample order: a multi-trait output differs from its single-trait run"
fi

MT_GROUP="${MT_GROUP:-/opt/saige/logs/step2/mid_group_400.txt}"
if [ -e "$MT_GROUP" ]; then
expect_error region_with_multitrait "multi-trait region testing is not supported" <<EOF
groupFile: $MT_GROUP
annotationList:
  - "lof"
maxMAFList:
  - 0.01
models:
  - traitName: a
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/rg_a.txt
  - traitName: b
    modelFile: $MODEL
    varianceRatioFile: $VR
    outputFile: $WORK/rg_b.txt
$(common)
EOF
fi

# Legacy form keeps its original per-key messages.
expect_error legacy_missing_output "Config missing required key: outputFile" <<EOF
modelFile: $MODEL
varianceRatioFile: $VR
$(common)
EOF

echo
echo "=== $pass passed, $fail failed ==="
[ $fail -eq 0 ]
