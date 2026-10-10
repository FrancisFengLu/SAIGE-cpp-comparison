#!/bin/bash
# sparse_gate.sh <binary>: sparse-GRM cases of bingpu_test (bt: 10 binary, same list; bm: 8 binary own lists;
# qm: 8 quantitative own lists; mix: bm + qm) -- the new binary's GPU path against its own CPU path
# (useGPU false, block inverse), text output, route dumps; both under the R defaults unless noted.
set -u
BIN=$1
D=/opt/saige/data/bingpu_test
S=/opt/saige/logs/gpuprep/sparse
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
log(){ echo "$(date '+%F %T') $*"; }
models(){  # models <set> <traits...>: the models: block
  local set=$1; shift
  local root
  case $set in qm_sparse_fast) root=$D/qm/models/qm_sparse_fast/out;; *) root=$D/models/$set/out;; esac
  for t in "$@"; do
    echo "  - traitName: $t"; echo "    modelFile: $root/m/$t"
    echo "    varianceRatioFile: $root/mvr_$t.varianceRatio.txt"; echo "    outputFile: OUTDIR/$t.txt"
  done
}
BT="c01_1 c01_2 c02_1 c05_1 c05_2 c10_1 c10_2 c25_1 c50_1 c50_2"
BM="bm1 bm2 bm3 bm4 bm5 bm6 bm7 bm8"
QM="q1 q2 q3 q4 q5 q6 q7 q8"
run(){  # run <name> <gpu 0|1> <keys...> ; models block on stdin
  local name=$1 gpu=$2; shift 2
  local O=$S/$name/$([ $gpu = 1 ] && echo ${GDIR:-gpu} || echo cpu); mkdir -p $O/routes
  { echo "genoType: plink"; echo "plinkFile: $D/geno/bt"; echo "AlleleOrder: alt-first"
    echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 8"
    echo "MACCutoffforER: 4"; echo "blockSparseSigma: true"; echo "blockSparseSigmaSolveMargin: 0.01"
    [ $gpu = 1 ] && echo "useGPU: true" || echo "useGPU: false"
    for k in "$@"; do echo "$k"; done
    echo "models:"; sed "s|OUTDIR|$O|"; } > $O/cfg.yaml
  local t0=$(date +%s.%N)
  SAIGE_STEP2_ROUTE_DUMP=$O/routes $BIN $O/cfg.yaml > $O/log.txt 2>&1; local rc=$?
  log "$name $([ $gpu = 1 ] && echo gpu || echo cpu) rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuSparse: the sparse\|gpuSparse: not active' $O/log.txt | head -2 | tr '\n' ';' | cut -c1-200)"
}
case1(){ local name=$1 set=$2 traits=$3; shift 3
  for g in ${GLIST:-1 0}; do models $set $traits | run $name $g "$@"; done; }
case1 bt_nofast_def   sparse_nofast "$BT"
case1 bt_nofast_defF  sparse_nofast "$BT" "is_Firth_beta: true" "pCutoffforFirth: 0.01"
case1 bt_nofast_adj   sparse_nofast "$BT" "isnoadjCov: false" "impute_method: mean"
case1 bt_fast_def     sparse_fast   "$BT" "isFastTest: true"
case1 bt_fast_adjF    sparse_fast   "$BT" "isFastTest: true" "isnoadjCov: false" "impute_method: mean" "is_Firth_beta: true" "pCutoffforFirth: 0.01"
case1 bm_nofast_def   bm_sparse_nofast "$BM"
case1 bm_nofast_adjF  bm_sparse_nofast "$BM" "isnoadjCov: false" "impute_method: mean" "is_Firth_beta: true" "pCutoffforFirth: 0.01"
case1 bm_fast_def     bm_sparse_fast "$BM" "isFastTest: true"
case1 qm_nofast_def   qm_sparse_fast "$QM"
case1 qm_nofast_adj   qm_sparse_fast "$QM" "isnoadjCov: false" "impute_method: mean"
case1 qm_fast_def     qm_sparse_fast "$QM" "isFastTest: true"
case1 qm_fast_adj     qm_sparse_fast "$QM" "isFastTest: true" "isnoadjCov: false" "impute_method: mean"
for g in ${GLIST:-1 0}; do { models bm_sparse_fast $BM; models qm_sparse_fast $QM; } | run mix_fast_defF $g "isFastTest: true" "is_Firth_beta: true" "pCutoffforFirth: 0.01"; done
for g in ${GLIST:-1 0}; do { models bm_sparse_nofast $BM; models qm_sparse_fast $QM; } | run mix_nofast_def $g; done
log "sparse runs done"
