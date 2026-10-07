#!/bin/bash
# 09_step2_formats.sh -- step 2 on the other genotype formats, using the binary
# models b1 b2 b3 of 03_step1_binary.sh (b4 has its own sample set, which only
# PLINK and hard-call PGEN input support). Writes $WORK/formats/<format>/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
M=$WORK/step1_bin
O=$WORK/formats

run() {   # run <name> <genotype keys (YAML lines)>
  local N=$1 G=$2
  mkdir -p "$O/$N"
  {
    echo "$G"
    cat <<YAML
minMAC: 1
nThreads: 8
useGPU: true
models:
YAML
    for t in b1 b2 b3; do
      cat <<YAML
  - traitName: $t
    modelFile: $M/models/$t
    varianceRatioFile: $M/vr_$t.varianceRatio.txt
    outputFile: $O/$N/$t.txt
YAML
    done
  } > "$O/$N/step2.yaml"
  "$S2" "$O/$N/step2.yaml" > "$O/$N/step2.log" 2>&1
  echo "== $N: $(grep -h -m1 -E 'useGPU: refused.*' "$O/$N/step2.log" || echo 'GPU path used')"
  wc -l < "$O/$N/b1.txt"
}

run pgen "genoType: pgen
pgenFile: $D/geno.pgen
pvarFile: $D/geno.pvar
psamFile: $D/geno.psam
AlleleOrder: ref-first"

run pgen_dosage "genoType: pgen
pgenFile: $D/dosage.pgen
pvarFile: $D/dosage.pvar
psamFile: $D/dosage.psam
AlleleOrder: ref-first"

run bgen "genoType: bgen
bgenFile: $D/geno.bgen
bgenSampleFile: $D/geno.sample
AlleleOrder: ref-first"

run vcf "genoType: vcf
vcfFile: $D/geno.vcf.gz
vcfField: DS"
