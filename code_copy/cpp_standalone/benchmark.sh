#!/usr/bin/env bash
# ===========================================================================
# benchmark.sh — Run matching R (SAIGE docker) and C++ configs under
#                valgrind/callgrind + /usr/bin/time for memory & CPU comparison
#
# Usage:
#   bash benchmark.sh [all|sparse_x1|sparse_x1x2|dense_x1|dense_x1x2|ukb_ldl]
#   bash benchmark.sh all          # run every config pair
#   bash benchmark.sh sparse_x1    # run one config only
#
# Prerequisites:
#   - C++ binary built: make profile   (in this directory)
#   - Docker image:     wzhou88/saige:1.5.0.2
#   - valgrind installed
#
# Output goes to benchmark_results/<config_name>/{r,cpp}/
# ===========================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_BIN="${SCRIPT_DIR}/saige-null"
RESULTS_DIR="${SCRIPT_DIR}/benchmark_results"

# Shared data paths
PLINK_PREFIX="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"
PHENO_QUANT="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
PHENO_BINARY="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/input/pheno_1000samples.txt"
SPARSE_GRM="/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
SPARSE_GRM_IDS="${SPARSE_GRM}.sampleIDs.txt"
DOCKER_IMG="wzhou88/saige:1.5.0.2"

# Docker mount: map the entire storage root so all absolute paths work
DOCKER_MOUNT="-v /media/leelabsg-storage0:/media/leelabsg-storage0"

# UKB paths
UKB_PLINK="/media/leelabsg-storage0/DATA/UKBB/cal/pruned/output/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned"
UKB_FAM="/media/leelabsg-storage0/seokho/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam"
UKB_PHENO="/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random1000.tsv"

# ---------------------------------------------------------------------------
# Helper: run command with /usr/bin/time and capture peak RSS + wall time
# ---------------------------------------------------------------------------
run_timed() {
    local label="$1"; shift
    local logdir="$1"; shift
    mkdir -p "$logdir"
    echo "[${label}] Starting: $*"
    /usr/bin/time -v "$@" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    # Extract key metrics from /usr/bin/time -v output
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] Done. Logs in ${logdir}/"
    echo ""
}

# ---------------------------------------------------------------------------
# Helper: run C++ under callgrind
# ---------------------------------------------------------------------------
run_cpp_callgrind() {
    local label="$1"; shift
    local logdir="$1"; shift
    local config="$1"; shift
    mkdir -p "$logdir"
    echo "[${label}] callgrind: ./saige-null -c ${config}"
    /usr/bin/time -v \
        valgrind --tool=callgrind \
                 --callgrind-out-file="${logdir}/callgrind.out" \
                 --collect-jumps=yes \
                 --cache-sim=yes \
                 --branch-sim=yes \
                 "${CPP_BIN}" -c "${config}" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] Done. callgrind output: ${logdir}/callgrind.out"
    echo ""
}

# ---------------------------------------------------------------------------
# Helper: run R SAIGE under callgrind (via docker)
#   Callgrind attaches to the R process inside docker.
#   We use --tool=callgrind on the docker entrypoint.
# ---------------------------------------------------------------------------
run_r_timed() {
    local label="$1"; shift
    local logdir="$1"; shift
    mkdir -p "$logdir"
    echo "[${label}] R SAIGE (docker, timed): $*"
    /usr/bin/time -v \
        docker run --rm ${DOCKER_MOUNT} ${DOCKER_IMG} "$@" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] Done. Logs in ${logdir}/"
    echo ""
}

# ===========================================================================
# CONFIG 1: sparse_x1  (quantitative, sparse GRM, covariate: x1)
# ===========================================================================
run_sparse_x1() {
    local base="${RESULTS_DIR}/sparse_x1"

    # --- C++ ---
    run_cpp_callgrind "sparse_x1/cpp" "${base}/cpp" \
        "${SCRIPT_DIR}/config_sparse_x1.yaml"

    # --- R ---
    # Matches config_sparse_x1.yaml:
    #   trait=quantitative, covar=x1, sparse GRM, LOCO=FALSE, usePCGwithSparseGRM=FALSE
    run_r_timed "sparse_x1/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 \
        --LOCO=FALSE \
        --minMAFforGRM=0.01 \
        --skipModelFitting=FALSE \
        --tol=0.02 \
        --tolPCG=1e-5 \
        --maxiterPCG=500 \
        --maxiter=20 \
        --traceCVcutoff=0.0025 \
        --isCovariateOffset=FALSE \
        --isDiagofKinSetAsOne=TRUE \
        --useSparseGRMtoFitNULL=TRUE \
        --usePCGwithSparseGRM=FALSE \
        --sparseGRMFile="${SPARSE_GRM}" \
        --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}" \
        --numMarkersForVarRatio=30 \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# CONFIG 2: sparse_x1x2  (quantitative, sparse GRM, covariates: x1, x2)
# ===========================================================================
run_sparse_x1x2() {
    local base="${RESULTS_DIR}/sparse_x1x2"

    # --- C++ ---
    run_cpp_callgrind "sparse_x1x2/cpp" "${base}/cpp" \
        "${SCRIPT_DIR}/config_sparse_x1x2.yaml"

    # --- R ---
    run_r_timed "sparse_x1x2/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 \
        --LOCO=FALSE \
        --minMAFforGRM=0.01 \
        --skipModelFitting=FALSE \
        --tol=0.02 \
        --tolPCG=1e-5 \
        --maxiterPCG=500 \
        --maxiter=20 \
        --traceCVcutoff=0.0025 \
        --isCovariateOffset=FALSE \
        --isDiagofKinSetAsOne=TRUE \
        --useSparseGRMtoFitNULL=TRUE \
        --usePCGwithSparseGRM=FALSE \
        --sparseGRMFile="${SPARSE_GRM}" \
        --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}" \
        --numMarkersForVarRatio=30 \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# CONFIG 3: dense_x1  (quantitative, dense GRM, covariate: x1)
# ===========================================================================
run_dense_x1() {
    local base="${RESULTS_DIR}/dense_x1"

    # --- C++ ---
    run_cpp_callgrind "dense_x1/cpp" "${base}/cpp" \
        "${SCRIPT_DIR}/config_dense_x1.yaml"

    # --- R ---
    # Dense GRM = no sparse flags, no sparseGRMFile
    run_r_timed "dense_x1/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 \
        --LOCO=FALSE \
        --minMAFforGRM=0.01 \
        --skipModelFitting=FALSE \
        --tol=0.02 \
        --tolPCG=1e-5 \
        --maxiterPCG=500 \
        --maxiter=20 \
        --traceCVcutoff=0.0025 \
        --isCovariateOffset=FALSE \
        --isDiagofKinSetAsOne=TRUE \
        --numMarkersForVarRatio=30 \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# CONFIG 4: dense_x1x2  (quantitative, dense GRM, covariates: x1, x2)
# ===========================================================================
run_dense_x1x2() {
    local base="${RESULTS_DIR}/dense_x1x2"

    # --- C++ ---
    run_cpp_callgrind "dense_x1x2/cpp" "${base}/cpp" \
        "${SCRIPT_DIR}/config_dense_x1x2.yaml"

    # --- R ---
    run_r_timed "dense_x1x2/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${PLINK_PREFIX}" \
        --phenoFile="${PHENO_QUANT}" \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 \
        --LOCO=FALSE \
        --minMAFforGRM=0.01 \
        --skipModelFitting=FALSE \
        --tol=0.02 \
        --tolPCG=1e-5 \
        --maxiterPCG=500 \
        --maxiter=20 \
        --traceCVcutoff=0.0025 \
        --isCovariateOffset=FALSE \
        --isDiagofKinSetAsOne=TRUE \
        --numMarkersForVarRatio=30 \
        --IsOverwriteVarianceRatioFile=TRUE
}

# ===========================================================================
# CONFIG 5: ukb_ldl  (quantitative, dense GRM, UKB data, inv_normalize,
#                      covariates: Sex, Age, Batch, PC1-PC4)
# ===========================================================================
run_ukb_ldl() {
    local base="${RESULTS_DIR}/ukb_ldl"

    # --- C++ ---
    run_cpp_callgrind "ukb_ldl/cpp" "${base}/cpp" \
        "${SCRIPT_DIR}/config_ukb_ldl.yaml"

    # --- R ---
    # Note: R SAIGE reads the phenotype TSV directly.
    # Sex is categorical (qCovarCol), others are numeric.
    run_r_timed "ukb_ldl/R" "${base}/r" \
        step1_fitNULLGLMM.R \
        --plinkFile="${UKB_PLINK}" \
        --famFile="${UKB_FAM}" \
        --phenoFile="${UKB_PHENO}" \
        --phenoCol=f.30780.0.0 \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix="${base}/r/saige_out" \
        --nThreads=1 \
        --LOCO=FALSE \
        --minMAFforGRM=0.01 \
        --skipModelFitting=FALSE \
        --tol=0.02 \
        --tolPCG=1e-5 \
        --maxiterPCG=500 \
        --maxiter=20 \
        --traceCVcutoff=0.0025 \
        --isCovariateOffset=FALSE \
        --isDiagofKinSetAsOne=TRUE \
        --numMarkersForVarRatio=30 \
        --IsOverwriteVarianceRatioFile=TRUE
}


# ===========================================================================
# Dispatch
# ===========================================================================
TARGET="${1:-all}"

# Rebuild C++ in profile mode (debug symbols + -O1 for callgrind)
echo "=== Building C++ in profile mode ==="
(cd "${SCRIPT_DIR}" && make profile)
echo ""

case "$TARGET" in
    sparse_x1)    run_sparse_x1 ;;
    sparse_x1x2)  run_sparse_x1x2 ;;
    dense_x1)     run_dense_x1 ;;
    dense_x1x2)   run_dense_x1x2 ;;
    ukb_ldl)      run_ukb_ldl ;;
    all)
        run_sparse_x1
        run_sparse_x1x2
        run_dense_x1
        run_dense_x1x2
        run_ukb_ldl
        ;;
    *)
        echo "Unknown config: $TARGET"
        echo "Usage: $0 [all|sparse_x1|sparse_x1x2|dense_x1|dense_x1x2|ukb_ldl]"
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
printf "%-15s %-6s %12s %15s\n" "Config" "Lang" "Wall(s)" "Peak RSS(KB)"
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
            printf "%-15s %-6s %12s %15s\n" "$cfg_name" "$lang" "$wall" "$rss"
        fi
    done
done
echo "============================================"
echo ""
echo "Callgrind outputs in: ${RESULTS_DIR}/<config>/cpp/callgrind.out"
echo "View with: callgrind_annotate <file>  or  kcachegrind <file>"
