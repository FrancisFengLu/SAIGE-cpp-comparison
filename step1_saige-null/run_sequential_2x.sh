#!/usr/bin/env bash
set -euo pipefail

CPP_BIN="/media/leelabsg-storage0/seokho/seokho92-SAIGE-cpp-comparison/code_copy/cpp_standalone/saige-null"
CPP_CFG="/media/leelabsg-storage0/seokho/seokho92-SAIGE-cpp-comparison/code_copy/cpp_standalone/config_ukb_t2d_10k.yaml"

export LD_LIBRARY_PATH="/data/home/seokhojeong/.local/share/mamba/envs/saige-build/lib:/data/home/seokhojeong/.local/share/mamba/envs/saige-build/lib/R/lib"
export LD_PRELOAD="/data/home/seokhojeong/.local/share/mamba/envs/saige-build/lib/libtbb.so.12"
export R_HOME="/data/home/seokhojeong/.local/share/mamba/envs/saige-build/lib/R"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

OUTDIR="/data/ukb_cpp_test/output/sequential_2x"
mkdir -p "$OUTDIR"
DOCKER_IMG="wzhou88/saige:1.5.0.2"

for i in 1 2; do
  echo "=== C++ Run $i === $(date)"
  mkdir -p "${OUTDIR}/cpp_run${i}"
  /usr/bin/time -v "$CPP_BIN" -c "$CPP_CFG" \
    > "${OUTDIR}/cpp_run${i}/stdout.log" \
    2> "${OUTDIR}/cpp_run${i}/time.log" || true
  grep "TIMER" "${OUTDIR}/cpp_run${i}/stdout.log"
  grep -E "Final tau|Iterations" "${OUTDIR}/cpp_run${i}/stdout.log" | tail -2
  grep -E "wall clock|Maximum resident|User time" "${OUTDIR}/cpp_run${i}/time.log"
  echo ""

  echo "=== R Run $i === $(date)"
  mkdir -p "${OUTDIR}/r_run${i}"
  /usr/bin/time -v docker run --rm \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    -v /data:/data \
    "$DOCKER_IMG" step1_fitNULLGLMM.R \
    --bedFile=/data/ukb_cpp_test/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.bed \
    --bimFile=/data/ukb_cpp_test/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.bim \
    --famFile=/data/ukb_cpp_test/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam \
    --phenoFile=/data/ukb_cpp_test/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random10000.tsv \
    --phenoCol=X250.2 \
    --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
    --qCovarCol=Sex \
    --sampleIDColinphenoFile=eid \
    --traitType=binary \
    --outputPrefix="${OUTDIR}/r_run${i}/saige_out" \
    --nThreads=1 --LOCO=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE \
    > "${OUTDIR}/r_run${i}/stdout.log" \
    2> "${OUTDIR}/r_run${i}/time.log" || true
  grep -E "wall clock|Maximum resident|User time" "${OUTDIR}/r_run${i}/time.log"
  echo ""
done

echo "=== SUMMARY ==="
for i in 1 2; do
  echo "--- Run $i ---"
  grep -E "wall clock|User time|Maximum" "${OUTDIR}/cpp_run${i}/time.log" | sed 's/^/  C++: /'
  grep -E "wall clock|User time|Maximum" "${OUTDIR}/r_run${i}/time.log" | sed 's/^/  R:   /'
done
