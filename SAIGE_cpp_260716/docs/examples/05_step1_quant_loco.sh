#!/bin/bash
# 05_step1_quant_loco.sh -- step 1 for two quantitative traits with LOCO (leave one
# chromosome out), full GRM, GPU on. Writes $WORK/step1_qt/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
O=$WORK/step1_qt
mkdir -p "$O"

cat > "$O/step1.yaml" <<YAML
paths:
  plinkFile: $D/geno
  out_prefix: $O/models
  out_prefix_vr: $O/vr
  overwrite_varratio: true
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [q1, q2]
  covar_cols: [x1, x2]
fit:
  trait: quantitative
  loco: true                         # needs >= 2 autosomes in the .bim
  inv_normalize: false               # true: rank-normalise the phenotype first
  nthreads: 8
  use_gpu: true
YAML

"$S1" -c "$O/step1.yaml" > "$O/step1.log" 2>&1
grep -E "^Converged|^LOCO|GPU tier" "$O/step1.log"
ls "$O/models/q1"
