#!/bin/bash
# 11_hpc_concat.sh -- after every chromosome job of 10_hpc_chrom_job.sh has
# finished: one genome-wide file per trait, header written once.
# Writes $WORK/hpc/all/<trait>.txt.
set -euo pipefail
source "$(dirname "$0")/env.sh"
H=$WORK/hpc
mkdir -p "$H/all"
for t in q1 q2; do
  first=1
  for CHR in 1 2; do                 # 1..22 for a real genome
    if [ $first = 1 ]; then cat "$H/chr$CHR/$t.txt"; first=0
    else tail -n +2 "$H/chr$CHR/$t.txt"; fi
  done > "$H/all/$t.txt"
done
wc -l "$H"/all/*.txt
