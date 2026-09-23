#!/bin/bash
set -u
G=/opt/saige/logs/tg2_step2/data/g1m
cd /opt/saige/logs/mtfold
for P in 8 32; do
  for F in 0 1; do
    for R in 1 2; do
      T=g1m_P${P}_f${F}_r${R}
      echo "===== $T $(date +%T) ====="
      bash scripts/run_one.sh $T $G $F q:$P
      grep -E "Elapsed \(wall" runs/$T/time.txt
      rm -rf runs/$T/out
      df -h /opt/saige | tail -1
    done
  done
done
echo "===== done $(date +%T) ====="
