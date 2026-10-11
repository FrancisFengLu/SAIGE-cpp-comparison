#!/bin/bash
# R SAIGE 1.5.2 references with is_fastTest=TRUE (full-GRM models: the fast test recomputes the pairs whose
# first-pass p is below pval_cutoff_for_fastTest with the covariate-adjusted score; with is_noadjCov FALSE and
# no sparse GRM R turns is_fastTest off itself, so adjT = adj):
#   defT  R defaults + --is_fastTest=TRUE          defTF  + Firth
#   adjT  is_noadjCov FALSE + mean + fast test     adjTF  + Firth
# audit / bvs / qt12 (rdefaults datasets) and fam32 (32 traits, own missing-phenotype patterns).
set -u
source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
R2=/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R
B=${RDEF:-/opt/saige/logs/rdefaults}
G=${GPB:-/opt/saige/logs/gpuprep}
log(){ echo "$(date '+%F %T') $*"; }
s2(){  # s2 <outroot> <dataset> <plink prefix> <rda> <vr> <trait> <cfg> <extra flags...>
  local root=$1 ds=$2 pf=$3 rda=$4 vr=$5 t=$6 cfg=$7; shift 7
  local O=$root/R/$ds/$cfg; mkdir -p $O
  [ -f $O/$t.txt ] && return
  log "step2 $ds $cfg $t"
  Rscript $R2 --bedFile=$pf.bed --bimFile=$pf.bim --famFile=$pf.fam \
    --GMMATmodelFile=$rda --varianceRatioFile=$vr --SAIGEOutputFile=$O/$t.txt \
    --LOCO=FALSE --nThreads=1 "$@" > $O/$t.log 2>&1 || log "step2 $ds $cfg $t FAILED"
}
FT="--is_fastTest=TRUE"
ADJ="--is_noadjCov=FALSE --impute_method=mean"
FIR="--is_Firth_beta=TRUE --pCutoffforFirth=0.01"
A=/opt/saige/logs/audit
for t in b1 b2 b3 b4; do
  s2 $B audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t defT $FT
  s2 $B audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t defTF $FT $FIR
  s2 $B audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t adjT $FT $ADJ
  s2 $B audit $A/data/geno $A/r_step1/$t.rda $A/r_step1/$t.varianceRatio.txt $t adjTF $FT $ADJ $FIR
done
V=/opt/saige/logs/batch-vs-single
for t in b1 b2 b3 b4; do
  s2 $B bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t defT $FT
  s2 $B bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t defTF $FT $FIR
  s2 $B bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t adjT $FT $ADJ
  s2 $B bvs $V/in/g $V/r_$t/m.rda $V/r_$t/m.varianceRatio.txt $t adjTF $FT $ADJ $FIR
done
Q=$B/data/qt12
for t in q1 q2 b1 b2; do
  s2 $B qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t defT $FT
  s2 $B qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t adjT $FT $ADJ
  case $t in b*)
    s2 $B qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t defTF $FT $FIR
    s2 $B qt12 $Q/g $B/R/qt12_s1/$t.rda $B/R/qt12_s1/$t.varianceRatio.txt $t adjTF $FT $ADJ $FIR;;
  esac
done
F=$G/data/fam32
BINT="b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16"
QNT="q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16"
for t in $BINT $QNT; do
  s2 $G fam32 $F/g $G/R/fam32_s1/$t.rda $G/R/fam32_s1/$t.varianceRatio.txt $t defT $FT
  s2 $G fam32 $F/g $G/R/fam32_s1/$t.rda $G/R/fam32_s1/$t.varianceRatio.txt $t adjT $FT $ADJ
done
for t in $BINT; do
  s2 $G fam32 $F/g $G/R/fam32_s1/$t.rda $G/R/fam32_s1/$t.varianceRatio.txt $t defTF $FT $FIR
  s2 $G fam32 $F/g $G/R/fam32_s1/$t.rda $G/R/fam32_s1/$t.varianceRatio.txt $t adjTF $FT $ADJ $FIR
done
log "R runs done"
