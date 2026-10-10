#!/bin/bash
# GPU-path runs of a given binary on the rdefaults datasets x configurations, into OUT/<ds>/<cfg>/gpu
#   cpp_gpu_runs.sh <binary> <OUT> [datasets...]
set -u
BIN=$1; OUT=$2; shift 2
DS=${@:-"audit bvs qt12"}
B=/opt/saige/logs/rdefaults
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
log(){ echo "$(date '+%F %T') $*"; }
plinkOf(){ case $1 in audit) echo /opt/saige/logs/audit/data/geno;; bvs) echo /opt/saige/logs/batch-vs-single/in/g;; qt12) echo $B/data/qt12/g;; esac; }
traitsOf(){ case $1 in audit|bvs) echo "b1 b2 b3 b4";; qt12) echo "q1 q2 b1 b2";; esac; }
vrOf(){ case $1 in audit) echo /opt/saige/logs/audit/r_step1/$2.varianceRatio.txt;; bvs) echo /opt/saige/logs/batch-vs-single/r_$2/m.varianceRatio.txt;; qt12) echo $B/R/qt12_s1/$2.varianceRatio.txt;; esac; }
cfgsOf(){ case $1 in qt12) echo "def adj defF adjF";; *) echo "def defF adj adjF";; esac; }
keys(){
  case $1 in
    def)  ;;
    defF) echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adj)  echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjF) echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
common(){ echo "genoType: plink"; echo "plinkFile: $(plinkOf $1)"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"; }
for ds in $DS; do for cfg in $(cfgsOf $ds); do
  O=$OUT/$ds/$cfg/gpu; mkdir -p $O/routes
  { common $ds; keys $cfg; echo "useGPU: true"; echo "models:"
    for tt in $(traitsOf $ds); do
      echo "  - traitName: $tt"; echo "    modelFile: $B/models/$ds/$tt"
      echo "    varianceRatioFile: $(vrOf $ds $tt)"; echo "    outputFile: $O/$tt.txt"
    done; } > $O/cfg.yaml
  t0=$(date +%s.%N)
  SAIGE_STEP2_ROUTE_DUMP=$O/routes $BIN $O/cfg.yaml > $O/log.txt 2>&1; rc=$?
  log "$ds $cfg gpu rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuDeviceStats: on --\|gpuSparse:' $O/log.txt | head -2 | tr '\n' ';' | cut -c1-200)"
done; done
log "gpu runs done"
