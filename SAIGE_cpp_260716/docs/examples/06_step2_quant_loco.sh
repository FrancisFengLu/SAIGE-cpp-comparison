#!/bin/bash
# 06_step2_quant_loco.sh -- step 2 for the LOCO models of 05_step1_quant_loco.sh:
# one step-2 run per chromosome (--LOCO=TRUE --chrom N). Writes $WORK/step2_qt/chr<N>/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

for CHR in 1 2; do
  $SAIGE step2 \
    --step1Dir step1_qt \
    --plinkFile data/geno \
    --minMAF 0.01 \
    --minMAC 1 \
    --LOCO=TRUE \
    --chrom "$CHR" \
    --nThreads 8 \
    --useGPU \
    --outDir step2_qt/chr$CHR > step2_qt_chr$CHR.log 2>&1
  grep -E "LOCO: restricting|useGPU: refused|GPU coverage" step2_qt_chr$CHR.log || true
done
wc -l step2_qt/chr*/q*.txt
