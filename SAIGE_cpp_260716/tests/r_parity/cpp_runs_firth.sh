#!/bin/bash
# C++ step-2 runs for the Firth-status gate (firthstress: gen_firthstress.py, R references
# r_runs_firth.sh): the same null models as R (rda_to_arma.R of R's .rda), R's defaults +
# is_Firth_beta with pCutoffforFirth 0.05 and outputFirthStatus on, on
#   single   the CPU single-trait scalar path, one run per trait (b1..b8)
#   multi    the CPU multi-trait path, all 8
#   gpu      the GPU multi-trait path, fp64 (device SPA / Firth)
#   gpu32    the GPU path with gpuPrecisionFirth: fp32
#   gpusgs   gpu with outputFormat: sgs (converted back with sgs2txt)
#   cpu_def  single-trait runs without outputFirthStatus (the default output format)
#   cpp_runs_firth.sh <saige-step2 binary> [OUT]      (OUT default $FIRTHB/cpp)
set -u
BIN=$1
B=${FIRTHB:-/opt/saige/logs/integrate/firth}
OUT=${2:-$B/cpp}
SGS2TXT=$(dirname $BIN)/sgs2txt
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
RDA2ARMA=$(cd "$(dirname "$0")/../../step2_saige-step2/tools" && pwd)/rda_to_arma.R
TR="b1 b2 b3 b4 b5 b6 b7 b8"
log(){ echo "$(date '+%F %T') $*"; }
( source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
  for t in $TR; do
    M=$B/models/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $B/R/s1/$t.rda $M > $M/convert.log 2>&1 || { log "convert $t FAILED"; cat $M/convert.log; }
  done )
common(){
  echo "genoType: plink"; echo "plinkFile: $B/data/g"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
  echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.05"
}
model(){ echo "  - traitName: $1"; echo "    modelFile: $B/models/$1"
         echo "    varianceRatioFile: $B/R/s1/$1.varianceRatio.txt"; echo "    outputFile: $2/$1.txt"; }
run(){  # run <variant> [trait]
  local v=$1 t=${2:-}
  local O=$OUT/$v; mkdir -p $O
  local C=$O/cfg${t:+_$t}.yaml L=$O/log${t:+_$t}.txt
  { common
    [ $v != cpu_def ] && echo "outputFirthStatus: true"
    case $v in
      single|cpu_def) echo "modelFile: $B/models/$t"; echo "varianceRatioFile: $B/R/s1/$t.varianceRatio.txt"; echo "outputFile: $O/$t.txt";;
      *) case $v in gpu*) echo "useGPU: true";; esac
         case $v in gpu32) echo "gpuPrecisionFirth: fp32";; gpusgs) echo "outputFormat: sgs";; esac
         echo "models:"; for tt in $TR; do model $tt $O; done;;
    esac; } > $C
  local t0=$(date +%s.%N)
  $BIN $C > $L 2>&1; local rc=$?
  log "$v${t:+ $t} rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuFirth: device setup\|gpuFirth: ignored' $L | head -2 | tr '\n' ';' | cut -c1-200)"
  if [ $v = gpusgs ]; then $SGS2TXT $O/*.sgs > $O/sgs2txt.log 2>&1 || log "sgs2txt FAILED"; fi
}
for t in $TR; do run single $t; done
for t in $TR; do run cpu_def $t; done
run multi
run gpu
run gpu32
run gpusgs
log "firth runs done"
