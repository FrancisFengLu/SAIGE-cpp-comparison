#!/bin/bash
# R SAIGE 1.5.2 reference runs for the rdefaults gate.
#   step 1 (qt12 only; audit and batch-vs-single have theirs), then step 2 for
#   every dataset x configuration:
#     def   R defaults (is_noadjCov TRUE, impute best_guess, no Firth)
#     defF  + is_Firth_beta TRUE (binary only)
#     adj   is_noadjCov FALSE, impute mean
#     adjF  adj + Firth
set -u
source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
R1=/opt/saige/SAIGE-upstream/extdata/step1_fitNULLGLMM.R
R2=/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R
B=${RDEF:-/opt/saige/logs/rdefaults}
Q=$B/data/qt12
mkdir -p $B/R/qt12_s1
log(){ echo "$(date '+%F %T') $*"; }
s1(){  # s1 <trait> <type>
  local t=$1 ty=$2
  [ -f $B/R/qt12_s1/$t.rda ] && return
  log "step1 qt12 $t"
  Rscript $R1 --plinkFile=$Q/g --phenoFile=$Q/pheno.txt --phenoCol=$t --traitType=$ty \
    --covarColList=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12 --sampleIDColinphenoFile=IID \
    --LOCO=FALSE --nThreads=2 --IsOverwriteVarianceRatioFile=TRUE --outputPrefix=$B/R/qt12_s1/$t \
    > $B/R/qt12_s1/$t.log 2>&1 || log "step1 $t FAILED"
}
s1 q1 quantitative; s1 q2 quantitative; s1 b1 binary; s1 b2 binary
s2(){  # s2 <dataset> <plink prefix> <rda> <vr> <trait> <cfg> <extra flags...>
  local ds=$1 pf=$2 rda=$3 vr=$4 t=$5 cfg=$6; shift 6
  local O=$B/R/$ds/$cfg; mkdir -p $O
  [ -f $O/$t.txt ] && return
  log "step2 $ds $cfg $t"
  Rscript $R2 --bedFile=$pf.bed --bimFile=$pf.bim --famFile=$pf.fam \
    --GMMATmodelFile=$rda --varianceRatioFile=$vr --SAIGEOutputFile=$O/$t.txt \
    --LOCO=FALSE --nThreads=1 "$@" > $O/$t.log 2>&1 || log "step2 $ds $cfg $t FAILED"
}
ADJ="--is_noadjCov=FALSE --impute_method=mean"
FIR="--is_Firth_beta=TRUE --pCutoffforFirth=0.01"
# audit: 4 binary traits, p = 3
A=/opt/saige/logs/audit
for t in b1 b2 b3 b4; do
  s2 audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t def
  s2 audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t defF $FIR
  s2 audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t adj $ADJ
  s2 audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t adjF $ADJ $FIR
done
# batch-vs-single: 4 binary traits with missing phenotypes, p = 4, flips, 5% missing genotypes
V=/opt/saige/logs/batch-vs-single
for t in b1 b2 b3 b4; do
  s2 bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t def
  s2 bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t defF $FIR
  s2 bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t adj $ADJ
  s2 bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t adjF $ADJ $FIR
done
# qt12: 2 quantitative + 2 binary, p = 13
for t in q1 q2 b1 b2; do
  s2 qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t def
  s2 qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t adj $ADJ
  case $t in b*)
    s2 qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t defF $FIR
    s2 qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t adjF $ADJ $FIR;;
  esac
done
log "R runs done"
