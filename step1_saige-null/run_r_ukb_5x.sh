#!/usr/bin/env bash
# Run R SAIGE 5 times each for ukb_ldl (quantitative) and ukb_t2d (binary)
# on leelabsg02 with local /data/ukb_cpp_test/ data
set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DATA="/data/ukb_cpp_test"
BED="${DATA}/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.bed"
BIM="${DATA}/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.bim"
FAM="${DATA}/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam"
PHENO="${DATA}/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random1000.tsv"
DOCKER_IMG="wzhou88/saige:1.5.0.2"
RESULTS="${DATA}/benchmark_results"

run_r() {
    local label="$1"; shift
    local logdir="$1"; shift
    mkdir -p "$logdir"
    echo "[${label}] $(date +%H:%M:%S) Starting..."
    /usr/bin/time -v \
        docker run --rm \
        -e OMP_NUM_THREADS=1 \
        -e OPENBLAS_NUM_THREADS=1 \
        -v /data:/data \
        ${DOCKER_IMG} "$@" \
        > "${logdir}/stdout.log" \
        2> "${logdir}/time_stderr.log" \
        || true
    grep -E "(wall clock|Maximum resident)" "${logdir}/time_stderr.log" \
        | tee "${logdir}/time_summary.txt"
    echo "[${label}] $(date +%H:%M:%S) Done."
    echo ""
}

# ===== UKB LDL (quantitative) x5 =====
echo "=========================================="
echo "=== UKB LDL (quantitative) — 5 runs ==="
echo "=========================================="
for i in 1 2 3 4 5; do
    run_r "ukb_ldl/R/run${i}" "${RESULTS}/ukb_ldl/r_run${i}" \
        step1_fitNULLGLMM.R \
        --bedFile="${BED}" \
        --bimFile="${BIM}" \
        --famFile="${FAM}" \
        --phenoFile="${PHENO}" \
        --phenoCol=f.30780.0.0 \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix="${RESULTS}/ukb_ldl/r_run${i}/saige_out" \
        --nThreads=1 --LOCO=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE
done

# ===== UKB T2D (binary) x5 =====
echo "=========================================="
echo "=== UKB T2D (binary) — 5 runs ==="
echo "=========================================="
for i in 1 2 3 4 5; do
    run_r "ukb_t2d/R/run${i}" "${RESULTS}/ukb_t2d/r_run${i}" \
        step1_fitNULLGLMM.R \
        --bedFile="${BED}" \
        --bimFile="${BIM}" \
        --famFile="${FAM}" \
        --phenoFile="${PHENO}" \
        --phenoCol=X250.2 \
        --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 \
        --qCovarCol=Sex \
        --sampleIDColinphenoFile=eid \
        --traitType=binary \
        --outputPrefix="${RESULTS}/ukb_t2d/r_run${i}/saige_out" \
        --nThreads=1 --LOCO=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE
done

# ===== Summary =====
echo ""
echo "============================================"
echo "=== RESULTS SUMMARY ==="
echo "============================================"
for trait in ukb_ldl ukb_t2d; do
    echo "--- ${trait} ---"
    for i in 1 2 3 4 5; do
        logdir="${RESULTS}/${trait}/r_run${i}"
        wall=$(grep "wall clock" "${logdir}/time_summary.txt" 2>/dev/null | sed 's/.*: //')
        rss=$(grep "Maximum resident" "${logdir}/time_summary.txt" 2>/dev/null | awk '{print $NF}')
        tau=$(grep -A1 'Variance component' "${logdir}/stdout.log" 2>/dev/null | grep '^\[1\]' | tail -1 | sed 's/\[1\] //')
        printf "  run%d: wall=%-12s RSS=%-10s tau=%s\n" "$i" "$wall" "${rss}KB" "$tau"
    done
    echo ""
done
