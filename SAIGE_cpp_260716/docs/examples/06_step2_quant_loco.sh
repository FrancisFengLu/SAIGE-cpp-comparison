#!/bin/bash
# 06_step2_quant_loco.sh -- step 2 for the LOCO models of 05_step1_quant_loco.sh:
# one step-2 run per chromosome (LOCO: true + chrom). Writes $WORK/step2_qt/chr<N>/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
M=$WORK/step1_qt
O=$WORK/step2_qt

for CHR in 1 2; do
  mkdir -p "$O/chr$CHR"
  cat > "$O/chr$CHR/step2.yaml" <<YAML
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first
minMAF: 0.01
minMAC: 1
LOCO: true
chrom: "$CHR"                        # only this chromosome's markers are tested
nThreads: 8
useGPU: true
models:
  - traitName: q1
    modelFile: $M/models/q1
    varianceRatioFile: $M/vr_q1.varianceRatio.txt
    outputFile: $O/chr$CHR/q1.txt
  - traitName: q2
    modelFile: $M/models/q2
    varianceRatioFile: $M/vr_q2.varianceRatio.txt
    outputFile: $O/chr$CHR/q2.txt
YAML
  "$S2" "$O/chr$CHR/step2.yaml" > "$O/chr$CHR/step2.log" 2>&1
  grep -E "LOCO: restricting|useGPU: refused|GPU coverage" "$O/chr$CHR/step2.log" || true
done
wc -l "$O"/chr*/q*.txt
