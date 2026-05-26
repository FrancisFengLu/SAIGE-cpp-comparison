#!/usr/bin/env bash
# ===========================================================================
# benchmark.sh — Run matching R (SAIGE docker) and C++ configs with
#                /usr/bin/time for wall-clock + peak-memory comparison
#
# Usage:
#   bash benchmark.sh [all|sparse_x1|sparse_x1x2|dense_x1|dense_x1x2|ukb_ldl|small]
#   bash benchmark.sh small          # sparse_x1 + sparse_x1x2 only (fast)
#   bash benchmark.sh all            # run every config pair
#   bash benchmark.sh sparse_x1      # run one config only
#
# Prerequisites:
#   - saige-build mamba env at /data/home/seokhojeong/.local/share/mamba/envs/saige-build
#   - Docker image: wzhou88/saige:1.5.0.2
# ===========================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CPP_BIN="${SCRIPT_DIR}/saige-null"
RESULTS_DIR="${SCRIPT_DIR}/benchmark_results"

# Build environment
SAIGE_BUILD_ENV="/data/home/seokhojeong/.local/share/mamba/envs/saige-build"
export PATH="${SAIGE_BUILD_ENV}/bin:$PATH"
export CONDA_PREFIX="${SAIGE_BUILD_ENV}"
export LD_LIBRARY_PATH="${SAIGE_BUILD_ENV}/lib:${SAIGE_BUILD_ENV}/lib/R/lib:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="${SAIGE_BUILD_ENV}/lib/libtbb.so.12"

# Single-threaded BLAS for fair comparison
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Shared data paths
PLINK_PREFIX="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"
PHENO_QUANT="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
SPARSE_GRM="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
SPARSE_GRM_IDS="${SPARSE_GRM}.sampleIDs.txt"
DOCKER_IMG="wzhou88/saige:1.5.0.2"
DOCKER_MOUNT="-v /media/leelabsg-storage0:/media/leelabsg-storage0"

# UKB paths
UKB_PLINK="/media/leelabsg-storage0/DATA/UKBB/cal/pruned/output/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned"
UKB_FAM="/media/leelabsg-storage0/seokho/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam"
UKB_PHENO="/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random1000.tsv"

# ---------------------------------------------------------------------------
run_timed() {
    local label="$1"; shift
    local logdir="$1"; shift
    mkdir -p "$logdir"
    echo ""
    echo "================================================================"
    echo "[${label}] Starting: $*"
    echo "================================================================"
    /usr/bin/time -v "$@" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] Done."
}

run_r_timed() {
    local label="$1"; shift
    local logdir="$1"; shift
    mkdir -p "$logdir"
    echo ""
    echo "================================================================"
    echo "[${label}] R SAIGE (docker)"
    echo "================================================================"
    /usr/bin/time -v \
        docker run --rm \
        -e OMP_NUM_THREADS=1 \
        -e OPENBLAS_NUM_THREADS=1 \
        ${DOCKER_MOUNT} ${DOCKER_IMG} "$@" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] Done."
}

extract_tau() {
    local logdir="$1"
    local label="$2"
    if [ -f "${logdir}/stdout.log" ]; then
        echo -n "  [${label}] "
        grep -oE "Final tau: \[.*\]" "${logdir}/stdout.log" 2>/dev/null || \
        grep -oP "tau\s*=.*" "${logdir}/stdout.log" 2>/dev/null | tail -1 || \
        echo "(no tau found)"
    fi
}

# ===========================================================================
run_sparse_x1() {
    local base="${RESULTS_DIR}/sparse_x1"
    run_timed "sparse_x1/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_sparse_x1.yaml"
    run_r_timed "sparse_x1/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE --minMAFforGRM=0.01 --skipModelFitting=FALSE \
        --tol=0.02 --tolPCG=1e-5 --maxiterPCG=500 --maxiter=20 \
        --traceCVcutoff=0.0025 --isCovariateOffset=FALSE --isDiagofKinSetAsOne=TRUE \
        --useSparseGRMtoFitNULL=TRUE --usePCGwithSparseGRM=FALSE \
        --sparseGRMFile="${SPARSE_GRM}" --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}" \
        --numRandomMarkerforVarianceRatio=30 --IsOverwriteVarianceRatioFile=TRUE
}

run_sparse_x1x2() {
    local base="${RESULTS_DIR}/sparse_x1x2"
    run_timed "sparse_x1x2/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_sparse_x1x2.yaml"
    run_r_timed "sparse_x1x2/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE --minMAFforGRM=0.01 --skipModelFitting=FALSE \
        --tol=0.02 --tolPCG=1e-5 --maxiterPCG=500 --maxiter=20 \
        --traceCVcutoff=0.0025 --isCovariateOffset=FALSE --isDiagofKinSetAsOne=TRUE \
        --useSparseGRMtoFitNULL=TRUE --usePCGwithSparseGRM=FALSE \
        --sparseGRMFile="${SPARSE_GRM}" --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}" \
        --numRandomMarkerforVarianceRatio=30 --IsOverwriteVarianceRatioFile=TRUE
}

run_dense_x1() {
    local base="${RESULTS_DIR}/dense_x1"
    run_timed "dense_x1/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_dense_x1.yaml"
    run_r_timed "dense_x1/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE --minMAFforGRM=0.01 --skipModelFitting=FALSE \
        --tol=0.02 --tolPCG=1e-5 --maxiterPCG=500 --maxiter=20 \
        --traceCVcutoff=0.0025 --isCovariateOffset=FALSE --isDiagofKinSetAsOne=TRUE \
        --numRandomMarkerforVarianceRatio=30 --IsOverwriteVarianceRatioFile=TRUE
}

run_dense_x1x2() {
    local base="${RESULTS_DIR}/dense_x1x2"
    run_timed "dense_x1x2/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_dense_x1x2.yaml"
    run_r_timed "dense_x1x2/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE --minMAFforGRM=0.01 --skipModelFitting=FALSE \
        --tol=0.02 --tolPCG=1e-5 --maxiterPCG=500 --maxiter=20 \
        --traceCVcutoff=0.0025 --isCovariateOffset=FALSE --isDiagofKinSetAsOne=TRUE \
        --numRandomMarkerforVarianceRatio=30 --IsOverwriteVarianceRatioFile=TRUE
}

run_ukb_ldl() {
    local base="${RESULTS_DIR}/ukb_ldl"
    run_timed "ukb_ldl/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_ukb_ldl.yaml"
    run_r_timed "ukb_ldl/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --bedFile="${UKB_PLINK}.bed" \
        --bimFile="${UKB_PLINK}.bim" \
        --famFile="${UKB_FAM}" \
        --phenoFile="${UKB_PHENO}" \
        --phenoCol=f.30780.0.0 \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE \
        --skipVarianceRatioEstimation=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# CONFIG 6: ukb_t2d (binary trait, dense GRM, UKB data, PheCode X250.2)
# ===========================================================================
run_ukb_t2d() {
    local base="${RESULTS_DIR}/ukb_t2d"
    run_timed "ukb_t2d/cpp" "${base}/cpp" \
        "${CPP_BIN}" -c "${SCRIPT_DIR}/config_ukb_t2d.yaml"
    run_r_timed "ukb_t2d/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --bedFile="${UKB_PLINK}.bed" \
        --bimFile="${UKB_PLINK}.bim" \
        --famFile="${UKB_FAM}" \
        --phenoFile="${UKB_PHENO}" \
        --phenoCol=X250.2 \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=binary \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 --LOCO=FALSE \
        --skipVarianceRatioEstimation=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# Dispatch
# ===========================================================================
TARGET="${1:-all}"

echo "=== Building C++ in release mode ==="
(cd "${SCRIPT_DIR}" && make clean && make -j4)
echo ""

mkdir -p "${REPO_DIR}/output"

case "$TARGET" in
    sparse_x1)    run_sparse_x1 ;;
    sparse_x1x2)  run_sparse_x1x2 ;;
    dense_x1)     run_dense_x1 ;;
    dense_x1x2)   run_dense_x1x2 ;;
    ukb_ldl)      run_ukb_ldl ;;
    ukb_t2d)      run_ukb_t2d ;;
    small)
        run_sparse_x1
        run_sparse_x1x2
        ;;
    all)
        run_sparse_x1
        run_sparse_x1x2
        run_dense_x1
        run_dense_x1x2
        run_ukb_ldl
        run_ukb_t2d
        ;;
    *)
        echo "Unknown config: $TARGET"
        echo "Usage: $0 [all|sparse_x1|sparse_x1x2|dense_x1|dense_x1x2|ukb_ldl|ukb_t2d|small]"
        exit 1
        ;;
esac

# ===========================================================================
# Summary table
# ===========================================================================
echo ""
echo "============================================"
echo "=== BENCHMARK SUMMARY ==="
echo "============================================"
printf "%-15s %-6s %15s %15s\n" "Config" "Lang" "Wall" "Peak RSS(KB)"
echo "-------------------------------------------------------"
for cfg_dir in "${RESULTS_DIR}"/*/; do
    cfg_name="$(basename "$cfg_dir")"
    for lang_dir in "${cfg_dir}"{cpp,r}; do
        [ -d "$lang_dir" ] || continue
        lang="$(basename "$lang_dir")"
        summary="${lang_dir}/time_summary.txt"
        if [ -f "$summary" ]; then
            wall=$(grep "wall clock" "$summary" | sed 's/.*: //')
            rss=$(grep "Maximum resident" "$summary" | awk '{print $NF}')
            printf "%-15s %-6s %15s %15s\n" "$cfg_name" "$lang" "$wall" "$rss"
        fi
    done
done
echo "============================================"
echo ""
echo "--- Tau Comparison ---"
for cfg_dir in "${RESULTS_DIR}"/*/; do
    cfg_name="$(basename "$cfg_dir")"
    echo "[$cfg_name]"
    extract_tau "${cfg_dir}/cpp" "C++"
    extract_tau "${cfg_dir}/r" "R"
done
echo "============================================"
