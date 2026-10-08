#!/usr/bin/env bash
# Run one saige-null binary on the A/B cases of tests/r_rng/run_ab_vs_main.sh
# and byte-compare every output file against a reference run directory
# (default: the main 610ec5e9 outputs <REF>/<case>/base_out from that script).
# No R in the environment: R_HOME / SAIGE_R_HOME are unset for the run.
#
# usage: [BIN=..] [REF=..] [W=..] [NTHR=1] [TAG=new] [CMPTAG=<tag of an earlier run in W to compare to instead of REF>]
#        run_vs_ref.sh case [case ...]
# extra cases beyond run_ab_vs_main.sh: m_b_full_cpu, m_q_full_cpu (N=50k dense GRM, CPU, one trait);
# s_b_sgrm_build (sparse-GRM fit with no GRM file: the GRM is built in place, which runs
# the related-pair search reducer). BIN may be an older binary that embeds R if R_HOME is
# passed through KEEP_R_HOME=1.
# Paths are this machine's.
set -o pipefail
C=/home/francisfenglu4/miniforge3/envs/saige-build
BIN=${BIN:-/opt/saige/worktrees/step1-rng/SAIGE_cpp_260716/step1_saige-null/saige-null}
REF=${REF:-/opt/saige/logs/step1-rng/run}
W=${W:-/opt/saige/logs/step1-rng/noR/run}
TAG=${TAG:-new}
NTHR=${NTHR:-1}
S=/opt/saige/SAIGE-upstream/extdata/input
SP=$S/plinkforGRM_1000samples_10kMarkers
SPH=$S/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt
BT=/opt/saige/data/bingpu_test

cfg() { # case dir -> writes dir/cfg.yaml ; echoes extra CLI args
  local c=$1 d=$2; mkdir -p $d
  local plink csv ycols trait loco sparse gpu=0 sg="" sgi=""
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
    m_b_full_cpu) plink=/opt/saige/data/mid csv=$BT/pheno/bt.pheno.txt ycols="[c10_1]" trait=binary loco=false sparse=false;;
    s_b_sgrm_build) plink=$SP csv=$SPH ycols="[y_binary]" trait=binary loco=false sparse=true sg=none;;
    m_q_full_cpu) plink=/opt/saige/data/mid csv=$BT/qm/pheno/qm.pheno.txt ycols="[q1]" trait=quantitative loco=false sparse=false;;
    *) echo "unknown case $c" >&2; return 1;;
  esac
  if [[ $sg == none ]]; then sg=""; elif [[ $sparse == true ]]; then sg=/opt/saige/data/mid.fam10.sgrm.mtx sgi=/opt/saige/data/mid.fam10.sgrm.ids; fi
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
  nthreads: $NTHR
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

mkdir -p $W
for c in "$@"; do
  d=$W/$c/$TAG; rm -rf $d; extra=$(cfg $c $d) || continue
  mkdir -p $d/out
  if [[ ${KEEP_R_HOME:-0} == 1 ]]; then unset_r=(); else unset_r=(-u R_HOME -u SAIGE_R_HOME LD_LIBRARY_PATH=); fi
  ( cd $d && env "${unset_r[@]}" /usr/bin/time -f "%e s %M KB" $BIN -c cfg.yaml $extra > $d/run.log 2> $d/run.err ); echo $? > $d/rc
  if [[ -n ${CMPTAG:-} ]]; then ref=$W/$c/$CMPTAG/out; else ref=$REF/$c/base_out; fi
  (cd $ref && find . -type f | sort) > $d/files.ref
  (cd $d/out && find . -type f | sort) > $d/files.out
  nfiles=$(wc -l < $d/files.ref); ndiff=0; difflist=""
  cmp -s $d/files.ref $d/files.out || { difflist="FILELIST "; ndiff=1; }
  while read -r f; do cmp -s "$ref/$f" "$d/out/$f" || { ndiff=$((ndiff+1)); difflist+="$f "; }; done < $d/files.ref
  echo "$c $TAG nthr=$NTHR rc=$(cat $d/rc) files=$nfiles differ=$ndiff [$difflist] ref=$ref time=$(tail -1 $d/run.err)" | tee -a $W/SUMMARY.txt
done
