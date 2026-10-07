#!/bin/bash
# 07_sparse_grm.sh -- sparse-GRM workflow:
#   (a) build a sparse GRM from the genotypes (step 1, make_sparse_grm_only)
#   (b) fit three binary traits on it (step 1, use_sparse_grm_to_fit)
#   (c) test them (step 2, GPU on)
# Writes $WORK/sparse/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
O=$WORK/sparse
mkdir -p "$O/out"

# (a) sparse GRM only. Writes grm.mtx (MatrixMarket) + grm.ids (one IID per line).
cat > "$O/make_grm.yaml" <<YAML
paths:
  plinkFile: $D/geno
  out_prefix: $O/grm_run
  sparse_grm: $O/grm.mtx
  sparse_grm_ids: $O/grm.ids
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_col: b1
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  make_sparse_grm_only: true
  relatedness_cutoff: 0.05           # entries below this are dropped
  min_maf_grm: 0.01
YAML
"$S1" -c "$O/make_grm.yaml" > "$O/make_grm.log" 2>&1
head -3 "$O/grm.mtx"; head -2 "$O/grm.ids"

# (b) null models on the sparse GRM (it is read, not rebuilt, because the files exist).
cat > "$O/step1.yaml" <<YAML
paths:
  plinkFile: $D/geno                 # still needed: variance-ratio markers come from here
  out_prefix: $O/models
  out_prefix_vr: $O/vr
  sparse_grm: $O/grm.mtx
  sparse_grm_ids: $O/grm.ids
  overwrite_varratio: true
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [b1, b2, b3]
  covar_cols: [x1, x2]
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  use_sparse_grm_for_vr: true
  fast_test: true
YAML
"$S1" -c "$O/step1.yaml" > "$O/step1.log" 2>&1
grep -E "^Converged|\[sparse\] GRM" "$O/step1.log"

# (c) step 2
{
cat <<YAML
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first
minMAC: 1
nThreads: 8
useGPU: true
models:
YAML
for t in b1 b2 b3; do
cat <<YAML
  - traitName: $t
    modelFile: $O/models/$t
    varianceRatioFile: $O/vr_$t.varianceRatio.txt
    outputFile: $O/out/$t.txt
YAML
done
} > "$O/step2.yaml"
"$S2" "$O/step2.yaml" > "$O/step2.log" 2>&1
grep -E "useGPU: refused|gpuSparse|GPU coverage" "$O/step2.log" || true
wc -l "$O"/out/*.txt
