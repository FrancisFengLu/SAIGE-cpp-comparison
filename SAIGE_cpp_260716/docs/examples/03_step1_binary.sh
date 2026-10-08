#!/bin/bash
# 03_step1_binary.sh -- step 1 for four binary traits in one run, full (dense) GRM
# from the PLINK genotypes, GPU on. Writes $WORK/step1_bin/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
O=$WORK/step1_bin
mkdir -p "$O"

cat > "$O/step1.yaml" <<YAML
paths:
  plinkFile: $D/geno                 # .bed/.bim/.fam prefix
  out_prefix: $O/models              # one model directory per trait: $O/models/<trait>/
  out_prefix_vr: $O/vr               # variance ratios: $O/vr_<trait>.varianceRatio.txt
  overwrite_varratio: true           # allow re-running into the same place
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [b1, b2, b3, b4]
  covar_cols: [x1, x2]
fit:
  trait: binary
  loco: false
  nthreads: 8
  use_gpu: true
  firth_beta: true                   # stored in the model; step 2 can override
  p_cutoff_for_firth: 0.01
  spa_cutoff: 2.0
YAML

"$S1" -c "$O/step1.yaml" > "$O/step1.log" 2>&1
grep -E "^Converged|^Model artifact|^Variance ratio|GPU tier|GPU unavailable" "$O/step1.log"
ls "$O" "$O/models/b1"
