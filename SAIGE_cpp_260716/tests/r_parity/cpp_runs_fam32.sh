#!/bin/bash
# C++ step-2 runs on fam32 (tests/r_parity/gen_fam32.py, R references r_runs_fam32.sh): the 32 full-GRM
# models (every trait its own sample list) and the 8 sparse-GRM models, converted from R's .rda, on the
# three paths (single-trait scalar, CPU multi-trait, GPU multi-trait), under R's defaults and
# is_noadjCov FALSE + mean, Firth on / off, and for the sparse models fast test off / on.
#   cpp_runs_fam32.sh <saige-step2 binary>
set -u
BIN=$1
B=${GPB:-/opt/saige/logs/gpuprep}
Q=$B/data/fam32
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
RDA2ARMA=$(cd "$(dirname "$0")/../../step2_saige-step2/tools" && pwd)/rda_to_arma.R
log(){ echo "$(date '+%F %T') $*"; }
BINT="b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16"
QNT="q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16"
SP="b1 b2 b3 b4 q1 q2 q3 q4"
MTX=$B/R/fam32_sparse/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
IDS=$MTX.sampleIDs.txt
# ---- models: convert R's .rda once ----
( source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
  for t in $BINT $QNT; do
    M=$B/models/fam32/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $B/R/fam32_s1/$t.rda $M > $M/convert.log 2>&1 || { log "convert fam32 $t FAILED"; cat $M/convert.log; }
  done
  for t in $SP; do
    M=$B/models/fam32sp/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $B/R/fam32_s1sp/$t.rda $M $MTX $IDS > $M/convert.log 2>&1 || { log "convert fam32sp $t FAILED"; cat $M/convert.log; }
  done )
keys(){
  case $1 in
    def)   ;;
    defF)  echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adj)   echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjF)  echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    defT)  echo "isFastTest: true";;
    adjT)  echo "isFastTest: true"; echo "isnoadjCov: false"; echo "impute_method: mean";;
    defTF) echo "isFastTest: true"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
common(){
  echo "genoType: plink"; echo "plinkFile: $Q/g"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
  echo "blockSparseSigma: true"
}
run(){  # run <ds> <cfg> <path: single|multi|gpu> <traits...> [single: one trait]
  local ds=$1 cfg=$2 path=$3; shift 3
  local O=$B/cpp/$ds/$cfg/$path; mkdir -p $O
  local t=""; [ $path = single ] && t=$1
  local C=$O/cfg${t:+_$t}.yaml L=$O/log${t:+_$t}.txt
  { common; keys $cfg
    if [ $path = single ]; then
      echo "modelFile: $B/models/$ds/$t"
      echo "varianceRatioFile: $B/R/$([ $ds = fam32sp ] && echo fam32_s1sp || echo fam32_s1)/$t.varianceRatio.txt"
      echo "outputFile: $O/$t.txt"
    else
      [ $path = gpu ] && echo "useGPU: true"
      echo "models:"
      for tt in "$@"; do
        echo "  - traitName: $tt"; echo "    modelFile: $B/models/$ds/$tt"
        echo "    varianceRatioFile: $B/R/$([ $ds = fam32sp ] && echo fam32_s1sp || echo fam32_s1)/$tt.varianceRatio.txt"; echo "    outputFile: $O/$tt.txt"
      done
    fi; } > $C
  local t0=$(date +%s.%N)
  $BIN $C > $L 2>&1; local rc=$?
  log "$ds $cfg $path${t:+ $t} rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|gpuDeviceStats: on --\|gpuSparse: the sparse' $L | head -2 | tr '\n' ';' | cut -c1-160)"
}
for cfg in def adj; do
  for t in $BINT $QNT; do run fam32 $cfg single $t; done
  run fam32 $cfg multi $BINT $QNT
  run fam32 $cfg gpu $BINT $QNT
done
for cfg in defF adjF; do
  for t in $BINT; do run fam32 $cfg single $t; done
  run fam32 $cfg multi $BINT
  run fam32 $cfg gpu $BINT
done
for cfg in def adj defT adjT; do
  for t in $SP; do run fam32sp $cfg single $t; done
  run fam32sp $cfg multi $SP
  run fam32sp $cfg gpu $SP
done
for cfg in defF defTF; do
  for t in b1 b2 b3 b4; do run fam32sp $cfg single $t; done
  run fam32sp $cfg multi b1 b2 b3 b4
  run fam32sp $cfg gpu b1 b2 b3 b4
done
log "cpp runs done"
