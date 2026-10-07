#!/bin/bash
# 04_step2_binary.sh -- step 2 for the four binary traits of 03_step1_binary.sh,
# PLINK genotypes, GPU on, text output. Writes $WORK/step2_bin/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
M=$WORK/step1_bin
O=$WORK/step2_bin
mkdir -p "$O/out"                    # step 2 does not create the output directory

{
cat <<YAML
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first               # PLINK .bed: Allele2 in the output = the .bim A1 column
minMAF: 0
minMAC: 1
maxMissRate: 0.15
LOCO: false
isFirth: true                        # binary traits: Firth correction for p < pCutoffforFirth
pCutoffforFirth: 0.01
MACCutoffforER: 4                    # binary traits: exact test for MAC <= 4
nThreads: 8
useGPU: true                         # every GPU sub-switch defaults to on
outputFormat: text
models:
YAML
for t in b1 b2 b3 b4; do
cat <<YAML
  - traitName: $t
    modelFile: $M/models/$t
    varianceRatioFile: $M/vr_$t.varianceRatio.txt
    outputFile: $O/out/$t.txt
YAML
done
} > "$O/step2.yaml"

"$S2" "$O/step2.yaml" > "$O/step2.log" 2>&1
grep -E "useGPU|GPU coverage|device SPA|device Firth|device ER" "$O/step2.log" | head -20
wc -l "$O"/out/*.txt
head -3 "$O/out/b1.txt"
