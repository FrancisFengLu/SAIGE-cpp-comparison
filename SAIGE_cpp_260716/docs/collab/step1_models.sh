#!/bin/bash
# step1_models.sh -- the null models the test matrix needs (binary traits only).
#
#   (0) sparse GRM: reused when SPARSE_GRM and SPARSE_GRM_IDS exist (e.g. from FastSparseGRM),
#       otherwise built from STEP1_BFILE on all .fam samples (step 1, make_sparse_grm_only)
#   (1) full-GRM null model, all BIN_TRAITS in one multi-trait run   -> $OUT_ROOT/step1/full/
#   (2) sparse-GRM null model, all BIN_TRAITS in one multi-trait run -> $OUT_ROOT/step1/sparse/
#   (3) $OUT_ROOT/step1/sparse_nofast/: the same sparse models with isFastTest = false in
#       nullmodel.json. Not a new fit: symlinks to (2) plus a copy of nullmodel.json with the
#       flag changed. Step 2 no longer reads that flag from the model (it has its own
#       isFastTest, R's default false; run_matrix.sh writes it per cell), so this view only
#       keeps the recorded copy consistent with what the cell runs.
# Step-2 cells with P = 1, 8, 32, ... use the first P traits of BIN_TRAITS from these models.
#
# Environment (required): OUT_ROOT, STEP1_BFILE (genome-wide PLINK prefix used for the GRM),
#   PHENO (phenotype/covariate file), BIN_TRAITS (column names, space separated, or @file with
#   one name per line), COVARS (covariate columns, space separated; may be empty)
# Optional: QCOVARS (categorical covariates among COVARS), IID_COL (IID), NTHREADS (all cores),
#   STEP1_GPU (1), SPARSE_GRM + SPARSE_GRM_IDS, RELATEDNESS_CUTOFF (0.05, only when building).
# Each fit is skipped when its .done file exists. Logs: $OUT_ROOT/step1/*/step1.log (with
# /usr/bin/time -v at the end).
set -euo pipefail
source "$(dirname "$0")/common.sh"
: "${OUT_ROOT:?}" "${STEP1_BFILE:?}" "${PHENO:?}" "${BIN_TRAITS:?}"
COVARS=${COVARS:-}; QCOVARS=${QCOVARS:-}; IID_COL=${IID_COL:-IID}
NTHREADS=${NTHREADS:-$NPROC}; STEP1_GPU=${STEP1_GPU:-1}
[[ $BIN_TRAITS == @* ]] && BIN_TRAITS=$(grep -v '^\s*$' "${BIN_TRAITS#@}" | tr '\n' ' ')
read -r -a TR <<< "$BIN_TRAITS"
yl() { local IFS=,; echo "[$*]"; }        # yaml inline list
S=$OUT_ROOT/step1; mkdir -p "$S"
GPU=false; [ "$STEP1_GPU" = 1 ] && GPU=true
COVLINE="covar_cols: $(yl $COVARS)"; [ -n "$QCOVARS" ] && COVLINE="$COVLINE
  q_covar_cols: $(yl $QCOVARS)"

run_fit() {   # dir
  local d=$1
  [ -e "$d/step1.done" ] && { echo "skip $d (done)"; return; }
  ( cd "$d" && /usr/bin/time -v "$NULLBIN" -c step1.yaml > step1.log 2>&1 ) || { echo "step 1 failed: $d/step1.log"; tail -20 "$d/step1.log"; exit 1; }
  grep -c "^Converged: yes" "$d/step1.log" | xargs echo "  converged traits:"
  touch "$d/step1.done"
}

# (0) sparse GRM
if [ -n "${SPARSE_GRM:-}" ] && [ -s "$SPARSE_GRM" ] && [ -s "${SPARSE_GRM_IDS:-}" ]; then
  echo "sparse GRM: using $SPARSE_GRM"
else
  SPARSE_GRM=$S/grm/sparse_grm.mtx; SPARSE_GRM_IDS=$S/grm/sparse_grm.ids
  if [ ! -s "$SPARSE_GRM" ]; then
    mkdir -p "$S/grm"
    # all .fam samples, with a placeholder phenotype (the GRM does not depend on it)
    awk 'BEGIN{print "IID\ty"} {print $2 "\t" (NR % 2)}' "$STEP1_BFILE.fam" > "$S/grm/all_samples.txt"
    cat > "$S/grm/step1.yaml" <<YAML
paths:
  plinkFile: $STEP1_BFILE
  out_prefix: $S/grm/grm_run
  sparse_grm: $SPARSE_GRM
  sparse_grm_ids: $SPARSE_GRM_IDS
design:
  csv: $S/grm/all_samples.txt
  iid_col: IID
  y_col: y
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  make_sparse_grm_only: true
  relatedness_cutoff: ${RELATEDNESS_CUTOFF:-0.05}
  min_maf_grm: 0.01
  nthreads: $NTHREADS
YAML
    echo "building the sparse GRM"
    ( cd "$S/grm" && /usr/bin/time -v "$NULLBIN" -c step1.yaml > step1.log 2>&1 ) || { tail -20 "$S/grm/step1.log"; exit 1; }
  fi
fi
printf '%s\n%s\n' "$SPARSE_GRM" "$SPARSE_GRM_IDS" > "$S/sparse_grm.paths"

# (1) full GRM
mkdir -p "$S/full"
cat > "$S/full/step1.yaml" <<YAML
paths:
  plinkFile: $STEP1_BFILE
  out_prefix: $S/full/models
  out_prefix_vr: $S/full/vr
  overwrite_varratio: true
design:
  csv: $PHENO
  iid_col: $IID_COL
  y_cols: $(yl "${TR[@]}")
  $COVLINE
fit:
  trait: binary
  loco: false
  nthreads: $NTHREADS
  use_gpu: $GPU
  firth_beta: true
  p_cutoff_for_firth: 0.01
  spa_cutoff: 2.0
  fast_test: true
YAML
echo "step 1, full GRM, ${#TR[@]} traits"; run_fit "$S/full"

# (2) sparse GRM
mkdir -p "$S/sparse"
cat > "$S/sparse/step1.yaml" <<YAML
paths:
  plinkFile: $STEP1_BFILE
  out_prefix: $S/sparse/models
  out_prefix_vr: $S/sparse/vr
  sparse_grm: $SPARSE_GRM
  sparse_grm_ids: $SPARSE_GRM_IDS
  overwrite_varratio: true
design:
  csv: $PHENO
  iid_col: $IID_COL
  y_cols: $(yl "${TR[@]}")
  $COVLINE
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  use_sparse_grm_for_vr: true
  firth_beta: true
  p_cutoff_for_firth: 0.01
  spa_cutoff: 2.0
  fast_test: true
YAML
echo "step 1, sparse GRM, ${#TR[@]} traits"; run_fit "$S/sparse"

# (3) fast test off: same fit, flag flipped
for t in "${TR[@]}"; do
  src=$S/sparse/models/$t; dst=$S/sparse_nofast/models/$t
  [ -s "$src/nullmodel.json" ] || { echo "missing $src/nullmodel.json"; exit 1; }
  mkdir -p "$dst"
  for f in "$src"/*; do [ "$(basename "$f")" = nullmodel.json ] || ln -sfn "$f" "$dst/"; done
  sed 's/"isFastTest": *true/"isFastTest": false/' "$src/nullmodel.json" > "$dst/nullmodel.json"
  grep -q '"isFastTest": *false' "$dst/nullmodel.json" || { echo "could not set isFastTest in $dst"; exit 1; }
  ln -sfn "$S/sparse/vr_$t.varianceRatio.txt" "$S/sparse_nofast/vr_$t.varianceRatio.txt"
done
echo "models ready: $S/{full,sparse,sparse_nofast}/models/<trait>/"
