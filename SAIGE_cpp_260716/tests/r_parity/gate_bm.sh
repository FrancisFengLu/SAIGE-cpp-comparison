#!/bin/bash
# gate_bm.sh <new-bin> <base-bin>: bingpu_test bm (8 binary traits, 4 missing-phenotype patterns, N = 50,000,
# 26,060 markers, 3,640 with missing calls) under the four rdefaults configurations:
#   new GPU multi-trait (device stats on own sample lists) vs the same binary's CPU single-trait scalar path,
#   and vs the base binary's GPU multi-trait (host tail). Route bytes of the two GPU runs are compared too.
set -u
NEW=$1; BASE=$2
B=${OWNB:-/opt/saige/logs/ownstats/bm}
D=/opt/saige/data/bingpu_test
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
log(){ echo "$(date '+%F %T') $*"; }
keys(){
  case $1 in
    def)  ;;
    defF) echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adj)  echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjF) echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
common(){
  echo "genoType: plink"; echo "plinkFile: $D/geno/bt"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 8"
}
run(){  # run <bin> <cfg> <path: single|gpu> <label> [trait]
  local BIN=$1 cfg=$2 path=$3 lab=$4 t=${5:-}
  local O=$B/$cfg/$lab; mkdir -p $O
  local C=$O/cfg${t:+_$t}.yaml L=$O/log${t:+_$t}.txt
  { common; keys $cfg
    if [ $path = single ]; then
      echo "modelFile: $D/models/bm_full/out/m/$t"; echo "varianceRatioFile: $D/models/bm_full/out/mvr_$t.varianceRatio.txt"; echo "outputFile: $O/$t.txt"
    else
      echo "useGPU: true"
      echo "models:"
      for k in 1 2 3 4 5 6 7 8; do
        echo "  - traitName: bm$k"; echo "    modelFile: $D/models/bm_full/out/m/bm$k"
        echo "    varianceRatioFile: $D/models/bm_full/out/mvr_bm$k.varianceRatio.txt"; echo "    outputFile: $O/bm$k.txt"
      done
    fi; } > $C
  local t0=$(date +%s.%N)
  if [ $path = gpu ]; then mkdir -p $O/routes; SAIGE_STEP2_ROUTE_DUMP=$O/routes $BIN $C > $L 2>&1; else $BIN $C > $L 2>&1; fi
  local rc=$?
  log "$cfg $lab${t:+ $t} rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuDeviceStats: o' $L | head -2 | tr '\n' ';' | cut -c1-220)"
}
for cfg in def defF adj adjF; do
  while [ -e /opt/saige/logs/TIMING_LOCK ]; do sleep 30; done
  run $NEW $cfg gpu gpu_new
  run $BASE $cfg gpu gpu_base
  for k in 1 2 3 4 5 6 7 8; do run $NEW $cfg single single bm$k; done
done
log "bm runs done"
