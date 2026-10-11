#!/bin/bash
# C++ step-2 runs with the fast test on (isFastTest: true), the single-trait scalar path (P = 1) and the GPU
# multi-trait path, into OUT/<ds>/<cfg>/{single,gpu}:
#   audit / bvs / qt12 (rdefaults datasets, R references r_runs_fast.sh) x defT / defTF / adjT / adjTF
#   fam32 (32 traits, own missing-phenotype patterns; R references r_runs_fast.sh) x the same
#   fam32sp (8 traits, sparse GRM; R references r_runs_fam32.sh) x defT / adjT / defTF
# defT = R defaults + fast test; adjT = is_noadjCov FALSE + mean + fast test (no sparse GRM: R turns the fast
# test off, this port skips the recompute whose context equals the first pass); F = Firth.
#   cpp_fast_runs.sh <saige-step2 binary> <OUT> [datasets...]
set -u
BIN=$1; OUT=$2; shift 2
DS=${@:-"audit bvs qt12 fam32 fam32sp"}
B=${RDEF:-/opt/saige/logs/rdefaults}
G=${GPB:-/opt/saige/logs/gpuprep}
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
log(){ echo "$(date '+%F %T') $*"; }
BINT="b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16"
QNT="q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16"
plinkOf(){ case $1 in audit) echo /opt/saige/logs/audit/data/geno;; bvs) echo /opt/saige/logs/batch-vs-single/in/g;; qt12) echo $B/data/qt12/g;; fam32|fam32sp) echo $G/data/fam32/g;; esac; }
traitsOf(){  # traitsOf <ds> <cfg>
  case $1 in
    audit|bvs) echo "b1 b2 b3 b4";;
    qt12) case $2 in *F) echo "b1 b2";; *) echo "q1 q2 b1 b2";; esac;;
    fam32) case $2 in *F) echo "$BINT";; *) echo "$BINT $QNT";; esac;;
    fam32sp) case $2 in *F) echo "b1 b2 b3 b4";; *) echo "b1 b2 b3 b4 q1 q2 q3 q4";; esac;;
  esac; }
modelOf(){ case $1 in fam32|fam32sp) echo $G/models/$1/$2;; *) echo $B/models/$1/$2;; esac; }
vrOf(){ case $1 in audit) echo /opt/saige/logs/audit/r_step1/$2.varianceRatio.txt;; bvs) echo /opt/saige/logs/batch-vs-single/r_$2/m.varianceRatio.txt;;
  qt12) echo $B/R/qt12_s1/$2.varianceRatio.txt;; fam32) echo $G/R/fam32_s1/$2.varianceRatio.txt;; fam32sp) echo $G/R/fam32_s1sp/$2.varianceRatio.txt;; esac; }
cfgsOf(){ case $1 in fam32sp) echo "defT adjT defTF";; *) echo "defT defTF adjT adjTF";; esac; }
keys(){
  case $1 in
    defT)  echo "isFastTest: true";;
    defTF) echo "isFastTest: true"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
    adjT)  echo "isFastTest: true"; echo "isnoadjCov: false"; echo "impute_method: mean";;
    adjTF) echo "isFastTest: true"; echo "isnoadjCov: false"; echo "impute_method: mean"; echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01";;
  esac
}
common(){ echo "genoType: plink"; echo "plinkFile: $(plinkOf $1)"; echo "AlleleOrder: alt-first"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
  [ $1 = fam32sp ] && echo "blockSparseSigma: true"; true; }
for ds in $DS; do for cfg in $(cfgsOf $ds); do
  O=$OUT/$ds/$cfg/single; mkdir -p $O
  for t in $(traitsOf $ds $cfg); do
    { common $ds; keys $cfg; echo "modelFile: $(modelOf $ds $t)"; echo "varianceRatioFile: $(vrOf $ds $t)"; echo "outputFile: $O/$t.txt"; } > $O/cfg_$t.yaml
    $BIN $O/cfg_$t.yaml > $O/log_$t.txt 2>&1 || log "$ds $cfg single $t FAILED"
  done
  log "$ds $cfg single done"
  O=$OUT/$ds/$cfg/gpu; mkdir -p $O/routes
  { common $ds; keys $cfg; echo "useGPU: true"; echo "models:"
    for tt in $(traitsOf $ds $cfg); do
      echo "  - traitName: $tt"; echo "    modelFile: $(modelOf $ds $tt)"
      echo "    varianceRatioFile: $(vrOf $ds $tt)"; echo "    outputFile: $O/$tt.txt"
    done; } > $O/cfg.yaml
  t0=$(date +%s.%N)
  SAIGE_STEP2_ROUTE_DUMP=$O/routes $BIN $O/cfg.yaml > $O/log.txt 2>&1; rc=$?
  log "$ds $cfg gpu rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s $(grep -h 'useGPU: refused\|fast-test recompute (dense' $O/log.txt | head -2 | tr '\n' ';' | cut -c1-220)"
done; done
log "fast-test runs done"
