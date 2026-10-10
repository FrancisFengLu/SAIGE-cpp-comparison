#!/bin/bash
# fastrecomp gate driver: the new binary (fastrecomp) on
#   1. the fast-test configurations (single-trait scalar + GPU) -> fast/        (gate_table_fast.py)
#   2. the default configurations of the rdefaults datasets and fam32 (GPU) -> defaults/, byte-compared with
#      the gpuprep outputs (the default path is untouched)
#   3. the bingpu_test sparse cases, GPU only, GDIR=gpu_fastrecomp -> /opt/saige/logs/gpuprep/sparse/<case>/gpu_fastrecomp
#      (sparse_table.py GDIR=gpu_fastrecomp, against the CPU runs already there)
set -u
F=/opt/saige/logs/fastrecomp
T=/opt/saige/worktrees/fastrecomp/SAIGE_cpp_260716/tests/r_parity
NEW=/opt/saige/worktrees/fastrecomp/SAIGE_cpp_260716/step2_saige-step2/saige-step2
G=/opt/saige/logs/gpuprep
log(){ echo "$(date '+%F %T') $*"; }
# ---- 1. fast test on ----
$T/cpp_fast_runs.sh $NEW $F/fast 2>&1 | tee $F/fast_runs.log
# ---- 2. defaults unchanged ----
$T/cpp_gpu_runs.sh $NEW $F/defaults 2>&1 | tee $F/default_runs.log
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
BINT="b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16"
QNT="q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16"
keys(){ case $1 in def) ;; defF) echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  adj) echo "isnoadjCov: false"; echo "impute_method: mean";;
  adjF) echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";; esac; }
for cfg in def adj defF adjF; do
  case $cfg in *F) TR="$BINT";; *) TR="$BINT $QNT";; esac
  O=$F/defaults/fam32/$cfg/gpu; mkdir -p $O
  { echo "genoType: plink"; echo "plinkFile: $G/data/fam32/g"; echo "AlleleOrder: alt-first"
    echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
    echo "blockSparseSigma: true"; keys $cfg; echo "useGPU: true"; echo "models:"
    for tt in $TR; do echo "  - traitName: $tt"; echo "    modelFile: $G/models/fam32/$tt"
      echo "    varianceRatioFile: $G/R/fam32_s1/$tt.varianceRatio.txt"; echo "    outputFile: $O/$tt.txt"; done; } > $O/cfg.yaml
  $NEW $O/cfg.yaml > $O/log.txt 2>&1; log "fam32 $cfg gpu rc=$?"
done | tee -a $F/default_runs.log
{
echo "## defaults (fast test off): new binary's GPU outputs byte-compared with gpuprep's"
n=0; d=0
for f in $(cd $F/defaults && find . -name '*.txt' -path '*/gpu/*' ! -name log.txt | sort); do
  case $f in ./fam32/*) ref=$G/cpp/${f#./};; *) ref=$G/rparity/new2/${f#./};; esac
  n=$((n+1)); cmp -s $F/defaults/$f $ref || { d=$((d+1)); echo "DIFF $f"; }
done
echo "$n files compared, $d differ"
} | tee $F/defaults_cmp.md
# ---- 3. sparse cases, GPU only ----
GLIST=1 GDIR=gpu_fastrecomp $T/sparse_gate.sh $NEW 2>&1 | tee $F/sparse_runs.log
GDIR=gpu_fastrecomp python3 $T/sparse_table.py > $F/sparse_table_fastrecomp.md 2>&1
# ---- tables ----
python3 $T/gate_table_fast.py $F/fast > $F/gate_table_fast.md 2>&1
log "gates done"
