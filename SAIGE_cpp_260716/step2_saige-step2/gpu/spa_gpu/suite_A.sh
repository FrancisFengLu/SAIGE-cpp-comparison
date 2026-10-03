#!/bin/bash
# suite_A.sh -- one pause window: the tail study on the device, the main
# accuracy gate for both implementations on identical pairs, and the same
# gate with the two other erfc modes (to show what they would cost).
# Run via: bash run_gate.sh LOG --script suite_A.sh
set -u
OUT=/opt/saige/logs/spa-gpu-lib
MAIN="--n 50000 --markers 4000 --threads 4 --blocks 256 --seed 11"
t() { echo; echo "=================== $(date +%T) $* ==================="; }
t "erfc sweep on the device";                    ./spa_gpu_test --erfc-sweep
t "main gate, this library";                     ./spa_gpu_test $MAIN --bench --tsv $OUT/gate_main_mine.tsv
t "main gate, integrator kernel 201c747b";       ./spa_gpu_test_integ --impl integ $MAIN --bench --tsv $OUT/gate_main_integ.tsv
t "main gate, this library, erfcMode 0 (CUDA erfc)";  ./spa_gpu_test $MAIN --erfc-mode 0 --tsv $OUT/gate_main_mine_erfc0.tsv
t "main gate, this library, erfcMode 2 (literal port)"; ./spa_gpu_test $MAIN --erfc-mode 2 --tsv $OUT/gate_main_mine_erfc2.tsv
