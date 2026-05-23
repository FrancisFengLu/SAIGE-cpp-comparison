#!/usr/bin/env bash
# Whole-UKB Step-1 benchmark sweep for one trait on one node.
# Sweep: nThreads 64 → 32 → 16, R Docker first then cpp-cpu at each setting.
# Usage:  ./run_whole_ukb_sweep.sh {ldl|t2d}
set -euo pipefail
TRAIT="${1:-}"
if [ "$TRAIT" != "ldl" ] && [ "$TRAIT" != "t2d" ]; then
  echo "Usage: $0 ldl|t2d" >&2; exit 1
fi
DATE=2026-05-06
PROJECT=/media/leelabsg-storage0/seokho/seokho92-SAIGE-cpp-comparison/code_copy/cpp_standalone
OUTBASE="${PROJECT}/benchmark_results/ukb_${TRAIT}_whole_recap_${DATE}"
CONFIG="${PROJECT}/configs_whole/config_ukb_${TRAIT}_whole.yaml"
LOCAL_BENCH=/data/UKB_imp/ukb_step1_bench
PRUNED_STEM=ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned
WHOLE_PHENO_TSV=Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_whole.tsv
MAPPED_FAM=${PRUNED_STEM}.mapped.fam
SAIGE_ENV=/data/home/seokhojeong/.local/share/mamba/envs/saige-build
DOCKER_IMG=wzhou88/saige:1.5.1

case "$TRAIT" in
  ldl) PHENO_COL=f.30780.0.0; TRAIT_TYPE=quantitative; EXTRA_R_ARGS="--invNormalize=TRUE" ;;
  t2d) PHENO_COL=XX250.2;     TRAIT_TYPE=binary;       EXTRA_R_ARGS="" ;;
  *) echo "Unknown trait: $TRAIT"; exit 1 ;;
esac

run_one() {
  local impl="$1" nthr="$2"
  local logdir="${OUTBASE}/${impl}_n${nthr}"
  mkdir -p "$logdir"
  # Resume guard: skip if a previous successful run produced its completion marker
  # R Docker prints "Final  <tau0> <tau1> :"; cpp-cpu prints "Final tau: [..]".
  if [ -f "${logdir}/stdout.log" ] && grep -qE "(^Final[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+|Final tau:[[:space:]]*\[)" "${logdir}/stdout.log" 2>/dev/null; then
    printf '[%s] SKIP %s n%s (already complete)\n' "$(date +%FT%T)" "$impl" "$nthr"
    grep -oE "Final tau:?[[:space:]]*\[[^]]*\]|^Final[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+" "${logdir}/stdout.log" | tail -1 || true
    return 0
  fi
  printf '\n========== [%s] START %s n%s on %s ==========\n' "$(date +%FT%T)" "$impl" "$nthr" "$(hostname)"
  case "$impl" in
    rdocker)
      /usr/bin/time -v docker run --rm \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -v ${LOCAL_BENCH}:/data_local \
        -v /media/leelabsg-storage0:/media/leelabsg-storage0 \
        "$DOCKER_IMG" step1_fitNULLGLMM.R \
        --bedFile=/data_local/${PRUNED_STEM}.bed \
        --bimFile=/data_local/${PRUNED_STEM}.bim \
        --famFile=/data_local/${MAPPED_FAM} \
        --phenoFile=/data_local/${WHOLE_PHENO_TSV} \
        --phenoCol=${PHENO_COL} \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=${TRAIT_TYPE} \
        ${EXTRA_R_ARGS} \
        --outputPrefix="${logdir}/saige_out" \
        --nThreads=${nthr} --LOCO=FALSE \
        --skipVarianceRatioEstimation=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE \
        > "${logdir}/stdout.log" 2> "${logdir}/time_stderr.log" || true
      ;;
    cpp_cpu)
      # Generate a per-run YAML by substituting out_prefix in the template (the binary's
      # multi-`-o paths.*` override is buggy: it wipes the entire paths map).
      local run_cfg="${logdir}/run_config.yaml"
      sed -e "s|^  out_prefix:.*|  out_prefix: ${logdir}/saige_out|" \
          -e "s|^  out_prefix_vr:.*|  out_prefix_vr: ${logdir}/saige_out_vr|" \
          "$CONFIG" > "$run_cfg"
      export PATH="${SAIGE_ENV}/bin:${PATH}"
      export LD_LIBRARY_PATH="${SAIGE_ENV}/lib:${SAIGE_ENV}/lib/R/lib:${LD_LIBRARY_PATH:-}"
      export LD_PRELOAD="${SAIGE_ENV}/lib/libtbb.so.12"
      /usr/bin/time -v "${PROJECT}/saige-null" -c "$run_cfg" -t "${nthr}" \
        > "${logdir}/stdout.log" 2> "${logdir}/time_stderr.log" || true
      unset LD_PRELOAD
      ;;
  esac
  printf '[%s] DONE %s n%s\n' "$(date +%FT%T)" "$impl" "$nthr"
  grep -E "wall clock|Maximum resident" "${logdir}/time_stderr.log" | tee "${logdir}/time_summary.txt"
  grep -oE "Final tau:?\s*\[[^]]*\]|tau\s*=\s*[^ ]+" "${logdir}/stdout.log" | tail -1 || true
}

mkdir -p "${OUTBASE}"
echo "[$(date +%FT%T)] sweep start: trait=${TRAIT} host=$(hostname) → ${OUTBASE}"
for thr in 64 32 16; do
  run_one rdocker "$thr"
  run_one cpp_cpu "$thr"
done
echo "[$(date +%FT%T)] sweep complete on $(hostname) for trait=${TRAIT}"
