#!/bin/bash
# R SAIGE 1.5.2 reference runs for fam32 (tests/r_parity/gen_fam32.py): 32 traits, every one
# with its own missing-phenotype pattern, 1000 families of 4.
#   step 1 full GRM for the 32 traits; a sparse GRM (createSparseGRM.R) and step 1 on it for
#   b1..b4 / q1..q4 (categorical variance ratios, as R's sparse example); then step 2:
#     full:   def (R defaults), adj (is_noadjCov FALSE + mean), defF / adjF (+ Firth, binary)
#     sparse: def / adj with is_fastTest FALSE (every marker on Sigma^-1), defT / adjT with is_fastTest TRUE
set -u
source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
R1=/opt/saige/SAIGE-upstream/extdata/step1_fitNULLGLMM.R
R2=/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R
RS=/opt/saige/SAIGE-upstream/extdata/createSparseGRM.R
B=${GPB:-/opt/saige/logs/gpuprep}
Q=$B/data/fam32
mkdir -p $B/R/fam32_s1 $B/R/fam32_s1sp
log(){ echo "$(date '+%F %T') $*"; }
COV=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
BIN="b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16"
QNT="q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16"
s1(){  # s1 <trait> <type>
  local t=$1 ty=$2
  [ -f $B/R/fam32_s1/$t.rda ] && return
  log "step1 fam32 $t"
  Rscript $R1 --plinkFile=$Q/g --phenoFile=$Q/pheno.txt --phenoCol=$t --traitType=$ty \
    --covarColList=$COV --sampleIDColinphenoFile=IID \
    --LOCO=FALSE --nThreads=2 --IsOverwriteVarianceRatioFile=TRUE --outputPrefix=$B/R/fam32_s1/$t \
    > $B/R/fam32_s1/$t.log 2>&1 || log "step1 $t FAILED"
}
for t in $BIN; do s1 $t binary; done
for t in $QNT; do s1 $t quantitative; done
# sparse GRM
SG=$B/R/fam32_sparse
mkdir -p $SG
if [ ! -f $SG/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx ]; then
  log "createSparseGRM"
  Rscript $RS --plinkFile=$Q/g --nThreads=4 --outputPrefix=$SG/sparseGRM --numRandomMarkerforSparseKin=2000 \
    --relatednessCutoff=0.125 > $SG/create.log 2>&1 || log "createSparseGRM FAILED"
fi
MTX=$SG/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
IDS=$MTX.sampleIDs.txt
s1sp(){  # s1sp <trait> <type>
  local t=$1 ty=$2
  [ -f $B/R/fam32_s1sp/$t.rda ] && return
  log "step1 sparse fam32 $t"
  Rscript $R1 --plinkFile=$Q/g --phenoFile=$Q/pheno.txt --phenoCol=$t --traitType=$ty \
    --covarColList=$COV --sampleIDColinphenoFile=IID \
    --sparseGRMFile=$MTX --sparseGRMSampleIDFile=$IDS --useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE --isCateVarianceRatio=TRUE \
    --LOCO=FALSE --nThreads=2 --IsOverwriteVarianceRatioFile=TRUE --outputPrefix=$B/R/fam32_s1sp/$t \
    > $B/R/fam32_s1sp/$t.log 2>&1 || log "step1 sparse $t FAILED"
}
for t in b1 b2 b3 b4; do s1sp $t binary; done
for t in q1 q2 q3 q4; do s1sp $t quantitative; done
s2(){  # s2 <ds> <rda> <vr> <trait> <cfg> <extra flags...>
  local ds=$1 rda=$2 vr=$3 t=$4 cfg=$5; shift 5
  local O=$B/R/$ds/$cfg; mkdir -p $O
  [ -f $O/$t.txt ] && return
  log "step2 $ds $cfg $t"
  Rscript $R2 --bedFile=$Q/g.bed --bimFile=$Q/g.bim --famFile=$Q/g.fam \
    --GMMATmodelFile=$rda --varianceRatioFile=$vr --SAIGEOutputFile=$O/$t.txt \
    --LOCO=FALSE --nThreads=1 "$@" > $O/$t.log 2>&1 || log "step2 $ds $cfg $t FAILED"
}
ADJ="--is_noadjCov=FALSE --impute_method=mean"
FIR="--is_Firth_beta=TRUE --pCutoffforFirth=0.01"
for t in $BIN $QNT; do
  s2 fam32 $B/R/fam32_s1/$t.rda $B/R/fam32_s1/$t.varianceRatio.txt $t def
  s2 fam32 $B/R/fam32_s1/$t.rda $B/R/fam32_s1/$t.varianceRatio.txt $t adj $ADJ
done
for t in $BIN; do
  s2 fam32 $B/R/fam32_s1/$t.rda $B/R/fam32_s1/$t.varianceRatio.txt $t defF $FIR
  s2 fam32 $B/R/fam32_s1/$t.rda $B/R/fam32_s1/$t.varianceRatio.txt $t adjF $ADJ $FIR
done
SPF="--sparseGRMFile=$MTX --sparseGRMSampleIDFile=$IDS"
for t in b1 b2 b3 b4 q1 q2 q3 q4; do
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t def  $SPF
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t adj  $SPF $ADJ
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t defT $SPF --is_fastTest=TRUE
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t adjT $SPF --is_fastTest=TRUE $ADJ
done
for t in b1 b2 b3 b4; do
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t defF $SPF $FIR
  s2 fam32sp $B/R/fam32_s1sp/$t.rda $B/R/fam32_s1sp/$t.varianceRatio.txt $t defTF $SPF --is_fastTest=TRUE $FIR
done
log "R runs done"
