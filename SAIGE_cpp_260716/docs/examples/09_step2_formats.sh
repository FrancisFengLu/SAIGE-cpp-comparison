#!/bin/bash
# 09_step2_formats.sh -- step 2 on the other genotype formats, using the binary
# models b1 b2 b3 of 03_step1_binary.sh (b4 has its own sample set, which only
# PLINK and hard-call PGEN input support). Writes $WORK/formats/<format>/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

run() {   # run <name> <genotype flags...>
  local N=$1; shift
  $SAIGE step2 "$@" \
    --step1Dir step1_bin \
    --phenoCol b1,b2,b3 \
    --minMAC 1 \
    --LOCO=FALSE \
    --is_fastTest=FALSE \
    --nThreads 8 \
    --useGPU \
    --outDir formats/$N > formats_$N.log 2>&1
  echo "== $N: $(grep -h -m1 -E 'useGPU: refused.*' formats_$N.log || echo 'GPU path used')"
  wc -l < formats/$N/b1.txt
}

run pgen        --pgenPrefix data/geno                                   # .pgen/.pvar/.psam
run pgen_dosage --pgenPrefix data/dosage
run bgen        --bgenFile data/geno.bgen --sampleFile data/geno.sample
run vcf         --vcfFile data/geno.vcf.gz --vcfField DS
