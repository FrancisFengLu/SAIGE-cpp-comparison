#!/bin/bash
# 07_sparse_grm.sh -- sparse-GRM workflow:
#   (a) build a sparse GRM from the genotypes (step 1, --makeSparseGRMOnly)
#   (b) fit three binary traits on it (step 1, --useSparseGRMtoFitNULL)
#   (c) test them (step 2, GPU on)
# Writes $WORK/sparse/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"
mkdir -p sparse

# (a) sparse GRM only. Writes grm.mtx (MatrixMarket) + grm.ids (one IID per line).
$SAIGE step1 \
  --plinkFile data/geno \
  --phenoFile data/pheno.txt \
  --phenoCol b1 \
  --traitType binary \
  --makeSparseGRMOnly \
  --sparseGRMFile sparse/grm.mtx \
  --sparseGRMSampleIDFile sparse/grm.ids \
  --relatednessCutoff 0.05 \
  --minMAFforGRM 0.01 \
  --outDir sparse/grm_run > sparse/make_grm.log 2>&1
# entries below --relatednessCutoff are dropped
head -3 sparse/grm.mtx; head -2 sparse/grm.ids

# (b) null models on the sparse GRM (it is read, not rebuilt, because the files exist).
#     --plinkFile is still needed: the variance-ratio markers come from it.
$SAIGE step1 \
  --plinkFile data/geno \
  --phenoFile data/pheno.txt \
  --phenoCol b1,b2,b3 \
  --covarColList x1,x2 \
  --traitType binary \
  --sparseGRMFile sparse/grm.mtx \
  --sparseGRMSampleIDFile sparse/grm.ids \
  --useSparseGRMtoFitNULL=TRUE \
  --useSparseGRMforVarRatio=TRUE \
  --is_fastTest=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outDir sparse/step1 > sparse/step1.log 2>&1
grep -E "^Converged|\[sparse\] GRM" sparse/step1.log

# (c) step 2. Step 2 has its own test flags (R's defaults); the model's copy of
# --is_fastTest is not read, so the fast test is asked for here again.
$SAIGE step2 \
  --step1Dir sparse/step1 \
  --plinkFile data/geno \
  --minMAC 1 \
  --LOCO=FALSE \
  --is_fastTest=TRUE \
  --nThreads 8 \
  --useGPU \
  --outDir sparse/step2 > sparse/step2.log 2>&1
grep -E "useGPU: refused|gpuSparse|GPU coverage" sparse/step2.log || true
wc -l sparse/step2/*.txt
