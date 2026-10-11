#!/bin/bash
# C++ step-2 runs for the covariate-count gate (branch covlimit): covp (gen_covp.py) at
# p = 4 9 13 24 40, the same null models as R (rda_to_arma.R of R's .rda, r_runs_covp.sh).
#   cpp_runs_covp.sh <saige-step2 binary> <OUT> [variants...]
# Configurations as in r_runs.sh (def / defF / adj / adjF). Variants:
#   single    CPU single-trait scalar path, one run per trait (b1..b4)
#   same      GPU multi-trait, b1 + b2 (one sample list; the non-own kernels)
#   own       GPU multi-trait, b1..b4 (b3 / b4 have their own lists: gpuOwnSampleSets)
#   sameown   same, with gpuSpaImpl: own (gpu/gpu_spa.cu; it takes one list only)
#   same32    same, gpuPrecisionSPA / gpuPrecisionFirth: fp32
#   own32     own, gpuPrecisionSPA / gpuPrecisionFirth: fp32
set -u
BIN=$1; OUT=$2; shift 2
VARS=${@:-"single same own sameown same32 own32"}
B=${COVLIM:-/opt/saige/logs/covlimit}
PS=${PS:-"4 9 13 24 40"}
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
RDA2ARMA=$(cd "$(dirname "$0")/../../step2_saige-step2/tools" && pwd)/rda_to_arma.R
log(){ echo "$(date '+%F %T') $*"; }
keys(){
  case $1 in
    def)  ;;
    defF) echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adj)  echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjF) echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
# ---- models: convert R's .rda once per p / trait ----
( source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
  for p in $PS; do for t in b1 b2 b3 b4; do
    M=$B/models/covp/p$p/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $B/R/covp/s1/p$p/$t.rda $M > $M/convert.log 2>&1 || { log "convert p=$p $t FAILED"; cat $M/convert.log; }
  done; done )
common(){
  echo "genoType: plink"; echo "plinkFile: $B/data/covp/g"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
}
model(){ echo "  - traitName: $2"; echo "    modelFile: $B/models/covp/p$1/$2"
         echo "    varianceRatioFile: $B/R/covp/s1/p$1/$2.varianceRatio.txt"; echo "    outputFile: $3/$2.txt"; }
run(){  # run <p> <cfg> <variant> [trait]
  local p=$1 cfg=$2 v=$3 t=${4:-}
  local O=$OUT/p$p/$cfg/$v; mkdir -p $O/routes
  local C=$O/cfg${t:+_$t}.yaml L=$O/log${t:+_$t}.txt
  { common; keys $cfg
    case $v in
      single) echo "modelFile: $B/models/covp/p$p/$t"; echo "varianceRatioFile: $B/R/covp/s1/p$p/$t.varianceRatio.txt"; echo "outputFile: $O/$t.txt";;
      *) echo "useGPU: true"
         case $v in sameown) echo "gpuSpaImpl: own";; esac
         case $v in same32|own32) echo "gpuPrecisionSPA: fp32"; echo "gpuPrecisionFirth: fp32";; esac
         echo "models:"
         case $v in own|own32) for tt in b1 b2 b3 b4; do model $p $tt $O; done;;
                    *)         for tt in b1 b2; do model $p $tt $O; done;; esac;;
    esac; } > $C
  local t0=$(date +%s.%N)
  SAIGE_STEP2_ROUTE_DUMP=$O/routes $BIN $C > $L 2>&1; local rc=$?
  log "p=$p $cfg $v${t:+ $t} rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuSpa: device setup failed\|gpuFirth: device setup\|gpuFirth: ignored\|device SPA: \|device Firth: ' $L | head -3 | tr '\n' ';' | cut -c1-260)"
}
for p in $PS; do for cfg in def defF adj adjF; do for v in $VARS; do
  if [ $v = single ]; then for t in b1 b2 b3 b4; do run $p $cfg single $t; done; else run $p $cfg $v; fi
done; done; done
log "cpp runs done"
