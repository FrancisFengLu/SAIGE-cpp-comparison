#!/bin/bash
# 05_step1_quant_loco.sh -- step 1 for two quantitative traits with LOCO (leave one
# chromosome out), full GRM, GPU on. Writes $WORK/step1_qt/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

$SAIGE step1 \
  --plinkFile data/geno \
  --phenoFile data/pheno.txt \
  --phenoCol q1,q2 \
  --covarColList x1,x2 \
  --traitType quantitative \
  --LOCO=TRUE \
  --invNormalize=FALSE \
  --nThreads 8 \
  --useGPU \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outDir step1_qt > step1_qt.log 2>&1
# --LOCO=TRUE needs >= 2 autosomes in the .bim; --invNormalize=TRUE rank-normalises first

grep -E "^Converged|^LOCO|GPU tier" step1_qt.log
ls step1_qt/models/q1
