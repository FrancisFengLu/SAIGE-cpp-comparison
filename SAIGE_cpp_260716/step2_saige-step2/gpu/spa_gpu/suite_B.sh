#!/bin/bash
# suite_B.sh -- one pause window: imbalanced-only traits with the extreme-z
# bands weighted up, and a biobank-N throughput run; both implementations.
# Run via: bash run_gate.sh LOG --script suite_B.sh
set -u
OUT=/opt/saige/logs/spa-gpu-lib
IMB="--n 50000 --markers 2500 --threads 4 --blocks 256 --seed 23 --traits 0.01,0.02,0.05,0.1 --zmix 0.3,0.3,0.2,0.2"
BIG="--n 400000 --markers 200 --threads 4 --blocks 128 --seed 5 --traits 0.5,0.01"
t() { echo; echo "=================== $(date +%T) $* ==================="; }
t "imbalanced gate, this library";               ./spa_gpu_test $IMB --bench --tsv $OUT/gate_imb_mine.tsv
t "imbalanced gate, integrator kernel";          ./spa_gpu_test_integ --impl integ $IMB --bench --tsv $OUT/gate_imb_integ.tsv
t "N=400k, this library";                        ./spa_gpu_test $BIG --bench --tsv $OUT/gate_big_mine.tsv
t "N=400k, integrator kernel";                   ./spa_gpu_test_integ --impl integ $BIG --bench --tsv $OUT/gate_big_integ.tsv
