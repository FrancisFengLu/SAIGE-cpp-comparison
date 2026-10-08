#!/usr/bin/env bash
# A/B bit-identity: base (main 610ec5e9, embedded R RNG) vs step1-rng (C++ RNG port).
# Usage: run_ab.sh case [case ...]   (cases defined below). nthreads 1 everywhere.
set -o pipefail
source /home/francisfenglu4/miniforge3/etc/profile.d/conda.sh
conda activate saige-build
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
RH="$CONDA_PREFIX/lib/R"
BASE=/opt/saige/worktrees/step1-rng-base/SAIGE_cpp_260716/step1_saige-null/saige-null
NEW=/opt/saige/worktrees/step1-rng/SAIGE_cpp_260716/step1_saige-null/saige-null
W=${W:-/opt/saige/logs/step1-rng/run}
S=/opt/saige/SAIGE-upstream/extdata/input
SP=$S/plinkforGRM_1000samples_10kMarkers
SPH=$S/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt
BT=/opt/saige/data/bingpu_test

cfg() { # case -> writes $W/$case/cfg.yaml ; echoes extra CLI args
  local c=$1 d=$W/$1; mkdir -p $d
  local plink csv ycols trait loco sparse gpu=0 nthr=${NTHR:-1} sg="" sgi=""
  case $c in
    s_b_cpu)      plink=$SP csv=$SPH ycols="[y_binary]" trait=binary loco=false sparse=false;;
    s_b_gpu)      plink=$SP csv=$SPH ycols="[y_binary]" trait=binary loco=false sparse=false gpu=1;;
    s_b_loco_cpu) plink=$SP csv=$SPH ycols="[y_binary]" trait=binary loco=true sparse=false;;
    s_b_loco_gpu) plink=$SP csv=$SPH ycols="[y_binary]" trait=binary loco=true sparse=false gpu=1;;
    s_q_cpu)      plink=$SP csv=$SPH ycols="[y_quantitative]" trait=quantitative loco=false sparse=false;;
    s_q_gpu)      plink=$SP csv=$SPH ycols="[y_quantitative]" trait=quantitative loco=false sparse=false gpu=1;;
    s_q_loco_cpu) plink=$SP csv=$SPH ycols="[y_quantitative]" trait=quantitative loco=true sparse=false;;
    s_q_loco_gpu) plink=$SP csv=$SPH ycols="[y_quantitative]" trait=quantitative loco=true sparse=false gpu=1;;
    m_b_sparse_mt) plink=/opt/saige/data/mid csv=$BT/pheno/bt.pheno.txt ycols="[c01_1, c05_1, c10_1, c25_1, c50_1]" trait=binary loco=false sparse=true;;
    m_q_sparse_mt) plink=/opt/saige/data/mid csv=$BT/qm/pheno/qm.pheno.txt ycols="[q1, q2, q3, q4]" trait=quantitative loco=false sparse=true;;
    m_b_full_gpu_mt) plink=/opt/saige/data/mid csv=$BT/pheno/bt.pheno.txt ycols="[c10_1, c50_1]" trait=binary loco=false sparse=false gpu=1;;
    m_q_full_gpu_mt) plink=/opt/saige/data/mid csv=$BT/qm/pheno/qm.pheno.txt ycols="[q1, q2]" trait=quantitative loco=false sparse=false gpu=1;;
    *) echo "unknown case $c" >&2; return 1;;
  esac
  if [[ $sparse == true ]]; then sg=/opt/saige/data/mid.fam10.sgrm.mtx sgi=/opt/saige/data/mid.fam10.sgrm.ids; fi
  cat > $d/cfg.yaml <<YAML
paths:
  plinkFile: $plink
  out_prefix: $d/out/m
  out_prefix_vr: $d/out/mvr
  sparse_grm: "$sg"
  sparse_grm_ids: "$sgi"
  overwrite_varratio: true
design:
  csv: $csv
  iid_col: IID
  covar_cols: [x1, x2]
  y_cols: $ycols
fit:
  trait: $trait
  loco: $loco
  nthreads: $nthr
  maxiter: 20
  tol: 0.02
  tolPCG: 1.0e-05
  maxiterPCG: 500
  nrun: 30
  num_markers_for_vr: 30
  min_maf_grm: 0.01
  use_sparse_grm_to_fit: $sparse
  use_sparse_grm_for_vr: $sparse
  use_pcg_with_sparse_grm: false
  multi_lockstep: false
  use_gpu: false
  fast_test: true
  overwrite_vr: true
YAML
  [[ $gpu == 1 ]] && echo "--gpu" || true
}

for c in "$@"; do
  d=$W/$c; rm -rf $d; extra=$(cfg $c) || continue
  for side in base new; do
    rm -rf $d/out; mkdir -p $d/out
    if [[ $side == base ]]; then
      ( cd $d && env R_HOME=$RH /usr/bin/time -f "%e s %M KB" $BASE -c cfg.yaml $extra > $d/$side.log 2> $d/$side.err ); rc=$?
    else
      ( cd $d && env -u R_HOME -u SAIGE_R_HOME /usr/bin/time -f "%e s %M KB" $NEW -c cfg.yaml $extra > $d/$side.log 2> $d/$side.err ); rc=$?
    fi
    echo $rc > $d/$side.rc
    mv $d/out $d/${side}_out
  done
  # compare every output file byte for byte
  (cd $d/base_out && find . -type f | sort) > $d/files.base
  (cd $d/new_out && find . -type f | sort) > $d/files.new
  nfiles=$(wc -l < $d/files.base); ndiff=0; difflist=""
  if ! cmp -s $d/files.base $d/files.new; then difflist="FILELIST "; ndiff=$((ndiff+1)); fi
  while read -r f; do
    cmp -s "$d/base_out/$f" "$d/new_out/$f" || { ndiff=$((ndiff+1)); difflist+="$f "; }
  done < $d/files.base
  # tau / iteration lines from the logs
  pat='[Tt]au|[Ii]ter|trace|Trace|ratio'
  grep -E "$pat" $d/base.log | grep -v TIMER > $d/base.tau; grep -E "$pat" $d/new.log | grep -v TIMER > $d/new.tau
  ntau=$(wc -l < $d/base.tau)
  if cmp -s $d/base.tau $d/new.tau; then taueq=same; else taueq=DIFF; fi
  # whole log minus wall-clock timer lines and the removed "[R] Embedded R runtime" line
  nlogdiff=$(diff <(grep -vE 'TIMER|Embedded R runtime|[0-9.]+ ?(s|sec|seconds|ms)\b' $d/base.log) <(grep -vE 'TIMER|[0-9.]+ ?(s|sec|seconds|ms)\b' $d/new.log) | grep -c '^[<>]')
  nrng=$(grep -c 'trace RNG seed' $d/new.log)
  ngpu=$(grep -ciE 'gpu.*(enabled|active|create|tier)' $d/new.log)
  echo "$c rc=$(cat $d/base.rc)/$(cat $d/new.rc) files=$nfiles differ=$ndiff [$difflist] tau_iter_lines=$ntau $taueq loglines_diff=$nlogdiff rng_seed_calls=$nrng gpu_lines=$ngpu time base=$(tail -1 $d/base.err) new=$(tail -1 $d/new.err)" | tee -a $W/SUMMARY.txt
done
