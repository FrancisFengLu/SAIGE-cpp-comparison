#!/bin/bash
# R SAIGE 1.5.2 reference runs for the Firth-status gate (firthstress, gen_firthstress.py).
#   r_runs_firth.sh [s1|s2|all]      (default all)
# Step 1: one binary null model per trait (b1..b8, covariates x1 x2 x3). Step 2: R's defaults
# + --is_Firth_beta=TRUE --pCutoffforFirth=0.05 (defF05), so hundreds of rare-variant pairs per
# trait go through the Firth fit, many of them to maxit. Paths: FIRTHB (default
# /opt/saige/logs/integrate/firth).
set -u
WHAT=${1:-all}
source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
R1=/opt/saige/SAIGE-upstream/extdata/step1_fitNULLGLMM.R
R2=/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R
B=${FIRTHB:-/opt/saige/logs/integrate/firth}
TR="b1 b2 b3 b4 b5 b6 b7 b8"
log(){ echo "$(date '+%F %T') $*"; }
s1(){  # s1 <trait>
  local t=$1
  local O=$B/R/s1
  mkdir -p $O
  [ -f $O/$t.rda ] && return
  log "step1 firthstress $t"
  Rscript $R1 --plinkFile=$B/data/g --phenoFile=$B/data/pheno.txt --phenoCol=$t --traitType=binary \
    --covarColList=x1,x2,x3 --sampleIDColinphenoFile=IID \
    --LOCO=FALSE --nThreads=2 --IsOverwriteVarianceRatioFile=TRUE --outputPrefix=$O/$t \
    > $O/$t.log 2>&1 || log "step1 firthstress $t FAILED"
}
s2(){  # s2 <trait> <cfg> <extra flags...>
  local t=$1 cfg=$2; shift 2
  local O=$B/R/s2/$cfg; mkdir -p $O
  [ -f $O/$t.txt ] && return
  log "step2 firthstress $cfg $t"
  Rscript $R2 --bedFile=$B/data/g.bed --bimFile=$B/data/g.bim --famFile=$B/data/g.fam \
    --GMMATmodelFile=$B/R/s1/$t.rda --varianceRatioFile=$B/R/s1/$t.varianceRatio.txt \
    --SAIGEOutputFile=$O/$t.txt --LOCO=FALSE --nThreads=1 "$@" > $O/$t.log 2>&1 || log "step2 firthstress $cfg $t FAILED"
}
par(){ while [ $(jobs -r | wc -l) -ge ${NPAR:-3} ]; do sleep 2; done; "$@" & }
if [ $WHAT = s1 ] || [ $WHAT = all ]; then
  for t in $TR; do par s1 $t; done; wait
  log "step1 done"
fi
if [ $WHAT = s2 ] || [ $WHAT = all ]; then
  for t in $TR; do par s2 $t defF05 --is_Firth_beta=TRUE --pCutoffforFirth=0.05; done; wait
  log "step2 done"
fi
