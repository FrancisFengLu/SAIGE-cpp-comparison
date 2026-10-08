#!/bin/bash
# 04_step2_binary.sh -- step 2 for the four binary traits of 03_step1_binary.sh,
# PLINK genotypes, GPU on, text output. Writes $WORK/step2_bin/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

$SAIGE step2 \
  --step1Dir step1_bin \
  --plinkFile data/geno \
  --minMAF 0 \
  --minMAC 1 \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth 0.01 \
  --nThreads 8 \
  --useGPU \
  --outDir step2_bin > step2_bin.log 2>&1

grep -E "useGPU|GPU coverage|device SPA|device Firth|device ER" step2_bin.log | head -20
wc -l step2_bin/*.txt
head -3 step2_bin/b1.txt
