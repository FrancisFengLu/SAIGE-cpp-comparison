#!/bin/bash
# C++ step-2 runs for the rdefaults gate: every dataset x configuration on the
# three paths (CPU single-trait scalar, CPU multi-trait batch, GPU multi-trait),
# on the SAME null models as R (rda_to_arma.R of R's .rda).
#   cpp_runs.sh <saige-step2 binary> [datasets...]
set -u
BIN=$1; shift
DS=${@:-"audit bvs qt12"}
B=${RDEF:-/opt/saige/logs/rdefaults}
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
RDA2ARMA=$(cd "$(dirname "$0")/../../step2_saige-step2/tools" && pwd)/rda_to_arma.R
log(){ echo "$(date '+%F %T') $*"; }
plinkOf(){ case $1 in audit) echo /opt/saige/logs/audit/data/geno;; bvs) echo /opt/saige/logs/batch-vs-single/in/g;; qt12) echo $B/data/qt12/g;; esac; }
traitsOf(){ case $1 in audit|bvs) echo "b1 b2 b3 b4";; qt12) echo "q1 q2 b1 b2";; esac; }
rdaOf(){ case $1 in audit) echo /opt/saige/logs/audit/r_step1/$2.rda;; bvs) echo /opt/saige/logs/batch-vs-single/r_$2/m.rda;; qt12) echo $B/R/qt12_s1/$2.rda;; esac; }
vrOf(){ case $1 in audit) echo /opt/saige/logs/audit/r_step1/$2.varianceRatio.txt;; bvs) echo /opt/saige/logs/batch-vs-single/r_$2/m.varianceRatio.txt;; qt12) echo $B/R/qt12_s1/$2.varianceRatio.txt;; esac; }
cfgsOf(){ case $1 in qt12) echo "def adj defF adjF";; *) echo "def defF adj adjF";; esac; }
keys(){  # config keys of a configuration
  case $1 in
    def)  ;;
    defF) echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adj)  echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjF) echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
# ---- models: convert R's .rda once per dataset / trait ----
( source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
  for ds in $DS; do for t in $(traitsOf $ds); do
    M=$B/models/$ds/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $(rdaOf $ds $t) $M > $M/convert.log 2>&1 || { log "convert $ds $t FAILED"; cat $M/convert.log; }
  done; done )
common(){  # common <ds>
  echo "genoType: plink"; echo "plinkFile: $(plinkOf $1)"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
}
run(){  # run <ds> <cfg> <path: single|multi|gpu> [trait]
  local ds=$1 cfg=$2 path=$3 t=${4:-}
  local O=$B/cpp/$ds/$cfg/$path; mkdir -p $O
  local C=$O/cfg${t:+_$t}.yaml L=$O/log${t:+_$t}.txt
  { common $ds; keys $cfg
    if [ $path = single ]; then
      echo "modelFile: $B/models/$ds/$t"; echo "varianceRatioFile: $(vrOf $ds $t)"; echo "outputFile: $O/$t.txt"
    else
      [ $path = gpu ] && echo "useGPU: true"
      echo "models:"
      for tt in $(traitsOf $ds); do
        echo "  - traitName: $tt"; echo "    modelFile: $B/models/$ds/$tt"
        echo "    varianceRatioFile: $(vrOf $ds $tt)"; echo "    outputFile: $O/$tt.txt"
      done
    fi; } > $C
  local t0=$(date +%s.%N)
  $BIN $C > $L 2>&1; local rc=$?
  log "$ds $cfg $path${t:+ $t} rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuDeviceStats: o' $L | head -2 | tr '\n' ';' | cut -c1-200)"
}
for ds in $DS; do for cfg in $(cfgsOf $ds); do
  for t in $(traitsOf $ds); do run $ds $cfg single $t; done
  run $ds $cfg multi
  run $ds $cfg gpu
done; done
log "cpp runs done"
