#!/bin/bash
# 03_step1_binary.sh -- step 1 for four binary traits in one run, full (dense) GRM
# from the PLINK genotypes, GPU on. Writes $WORK/step1_bin/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

$SAIGE step1 \
  --plinkFile data/geno \
  --phenoFile data/pheno.txt \
  --phenoCol b1,b2,b3,b4 \
  --covarColList x1,x2 \
  --traitType binary \
  --LOCO=FALSE \
  --nThreads 8 \
  --useGPU \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outDir step1_bin > step1_bin.log 2>&1

grep -E "^Converged|^Model artifact|^Variance ratio|GPU tier|GPU unavailable" step1_bin.log
ls step1_bin step1_bin/models/b1
