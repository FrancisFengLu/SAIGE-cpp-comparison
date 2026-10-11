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
  --LOCO=FALSE \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth 0.01 \
  --is_noadjCov=TRUE \
  --impute_method best_guess \
  --is_fastTest=FALSE \
  --nThreads 8 \
  --useGPU \
  --outDir step2_bin > step2_bin.log 2>&1
# The defaults are R SAIGE 1.5.2's (--LOCO=TRUE, --is_Firth_beta=FALSE, --is_noadjCov=TRUE,
# --impute_method best_guess, --is_fastTest=FALSE); the models of 03 were fitted without LOCO,
# so --LOCO=FALSE is needed, and the rest is written out so the run does not depend on them.

grep -E "useGPU|GPU coverage|device SPA|device Firth|device ER" step2_bin.log | head -20
wc -l step2_bin/*.txt
head -3 step2_bin/b1.txt
