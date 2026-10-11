#!/bin/bash
# R SAIGE 1.5.2 reference runs for the covariate-count gate (branch covlimit).
#   r_runs_covp.sh [s1|s2|all]      (default all)
# Data: gen_covp.py (covp: N = M = 5000, 4 binary traits, 40 covariates; covt: the timing
# set, N = M = 20,000, b1 / b2 only). Step 1 once per p in 4 9 13 24 40 (the first p - 1
# covariates, so p counts the intercept); step 2 on covp for every p x configuration
# (def / defF / adj / adjF as in r_runs.sh) x trait.
set -u
WHAT=${1:-all}
source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
R1=/opt/saige/SAIGE-upstream/extdata/step1_fitNULLGLMM.R
R2=/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R
B=${COVLIM:-/opt/saige/logs/covlimit}
PS=${PS:-"4 9 13 24 40"}
log(){ echo "$(date '+%F %T') $*"; }
covs(){ local p=$1; seq -s, 1 $((p-1)) | sed 's/\([0-9]\+\)/x\1/g'; }
s1(){  # s1 <set> <p> <trait>
  local ds=$1 p=$2 t=$3
  local O=$B/R/$ds/s1/p$p
  mkdir -p $O
  [ -f $O/$t.rda ] && return
  log "step1 $ds p=$p $t"
  Rscript $R1 --plinkFile=$B/data/$ds/g --phenoFile=$B/data/$ds/pheno.txt --phenoCol=$t --traitType=binary \
    --covarColList=$(covs $p) --sampleIDColinphenoFile=IID \
    --LOCO=FALSE --nThreads=2 --IsOverwriteVarianceRatioFile=TRUE --outputPrefix=$O/$t \
    > $O/$t.log 2>&1 || log "step1 $ds p=$p $t FAILED"
}
s2(){  # s2 <p> <trait> <cfg> <extra flags...>
  local p=$1 t=$2 cfg=$3; shift 3
  local O=$B/R/covp/s2/p$p/$cfg; mkdir -p $O
  [ -f $O/$t.txt ] && return
  log "step2 covp p=$p $cfg $t"
  Rscript $R2 --bedFile=$B/data/covp/g.bed --bimFile=$B/data/covp/g.bim --famFile=$B/data/covp/g.fam \
    --GMMATmodelFile=$B/R/covp/s1/p$p/$t.rda --varianceRatioFile=$B/R/covp/s1/p$p/$t.varianceRatio.txt \
    --SAIGEOutputFile=$O/$t.txt --LOCO=FALSE --nThreads=1 "$@" > $O/$t.log 2>&1 || log "step2 covp p=$p $cfg $t FAILED"
}
par(){ while [ $(jobs -r | wc -l) -ge ${NPAR:-3} ]; do sleep 2; done; "$@" & }
if [ $WHAT = s1 ] || [ $WHAT = all ]; then
  for p in $PS; do for t in b1 b2 b3 b4; do par s1 covp $p $t; done; done; wait
  [ -d $B/data/covt ] && { for p in $PS; do for t in b1 b2; do NPAR=2 par s1 covt $p $t; done; done; wait; }
  log "step1 done"
fi
if [ $WHAT = s2 ] || [ $WHAT = all ]; then
  ADJ="--is_noadjCov=FALSE --impute_method=mean"
  FIR="--is_Firth_beta=TRUE --pCutoffforFirth=0.01"
  for p in $PS; do for t in b1 b2 b3 b4; do
    par s2 $p $t def; par s2 $p $t defF $FIR; par s2 $p $t adj $ADJ; par s2 $p $t adjF $ADJ $FIR
  done; done; wait
  log "step2 done"
fi
