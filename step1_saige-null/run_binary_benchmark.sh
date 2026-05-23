#!/usr/bin/env bash
# Binary UKB T2D (X250.2) benchmark — matches LDL grid.
# Usage:  bash run_binary_benchmark.sh <size>   # size ∈ {1k, 10k, 100k, whole}
# Runs three paths in sequence and captures wall + peak-RSS for each.

set -u

SZ=${1:?usage: run_binary_benchmark.sh <1k|10k|100k|whole>}
HOST=$(hostname -s)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_BIN="${SCRIPT_DIR}/saige-null"
RESULTS="${SCRIPT_DIR}/benchmark_results/ukb_t2d_${SZ}"
CONFIG="${SCRIPT_DIR}/configs_binary/config_ukb_t2d_${SZ}.yaml"

mkdir -p "${RESULTS}/cpp_cpu_${HOST}" "${RESULTS}/cpp_gpu_${HOST}" "${RESULTS}/rdocker_${HOST}"

# --- env needed by saige-null (R + mamba libs) ---
ENV_BUILD=/data/home/seokhojeong/.local/share/mamba/envs/saige-build
export PATH="${ENV_BUILD}/bin:${PATH}"
export LD_LIBRARY_PATH="${ENV_BUILD}/lib:${ENV_BUILD}/lib/R/lib:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="${ENV_BUILD}/lib/libtbb.so.12"

echo "=== [$(date -Iseconds)] binary benchmark  size=${SZ}  host=${HOST} ==="
echo ""

# ----------- 1. cpp-cpu (float-fix) ------------------------------------------
OUT="${RESULTS}/cpp_cpu_${HOST}"
echo "---- [1/3] cpp-cpu ----  out=${OUT}"
/usr/bin/time -v "${CPP_BIN}" -c "${CONFIG}" --threads 16 \
    > "${OUT}/stdout.log" 2> "${OUT}/time_stderr.log" || echo "cpp-cpu exit=$?"
grep -E "(wall clock|Maximum resident|Exit)" "${OUT}/time_stderr.log"
grep -E "(Final tau|CONVERGED|Iterations)" "${OUT}/stdout.log" | head -3
echo ""

# ----------- 2. cpp-gpu (tier-3) ---------------------------------------------
# Skip GPU on nodes without NVIDIA card
if nvidia-smi --query-gpu=name --format=csv,noheader > /dev/null 2>&1; then
  OUT="${RESULTS}/cpp_gpu_${HOST}"
  echo "---- [2/3] cpp-gpu ----  out=${OUT}"
  SAIGE_GPU_TIER=3 /usr/bin/time -v "${CPP_BIN}" -c "${CONFIG}" --threads 16 --gpu \
      > "${OUT}/stdout.log" 2> "${OUT}/time_stderr.log" || echo "cpp-gpu exit=$?"
  grep -E "(wall clock|Maximum resident|Exit)" "${OUT}/time_stderr.log"
  grep -E "(Final tau|CONVERGED|tier=)" "${OUT}/stdout.log" | head -3
else
  echo "---- [2/3] cpp-gpu SKIPPED (no NVIDIA GPU on ${HOST}) ----"
fi
echo ""

# ----------- 3. R Docker with peak-RSS polling -------------------------------
OUT="${RESULTS}/rdocker_${HOST}"
echo "---- [3/3] R Docker ----  out=${OUT}"

UKB_PLINK=/media/leelabsg-storage0/DATA/UKBB/cal/pruned/output/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned
UKB_FAM=/media/leelabsg-storage0/seokho/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam
case $SZ in
  1k)    PHENO=/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random1000.tsv;    PHENO_COL=X250.2 ;;
  10k)   PHENO=/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random10000.tsv;   PHENO_COL=X250.2 ;;
  100k)  PHENO=/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_random100000.tsv;  PHENO_COL=XX250.2 ;;
  whole) PHENO=/media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/ukb_exp/Basic_trait_with_PheCode_table_ICD10_081123_whites_batches_famFiltered_whole.tsv;         PHENO_COL=XX250.2 ;;
esac

CID_FILE="${OUT}/container.cid"
PEAK_FILE="${OUT}/peak_rss.txt"
rm -f "${CID_FILE}" "${PEAK_FILE}"

# Launch docker detached so we get the CID for monitoring
docker run -d --cidfile "${CID_FILE}" \
    -e OMP_NUM_THREADS=16 -e OPENBLAS_NUM_THREADS=16 -e RCPP_PARALLEL_NUM_THREADS=16 \
    -v /media/leelabsg-storage0:/media/leelabsg-storage0 \
    wzhou88/saige:1.5.0.2 step1_fitNULLGLMM.R \
    --bedFile=${UKB_PLINK}.bed --bimFile=${UKB_PLINK}.bim --famFile=${UKB_FAM} \
    --phenoFile=${PHENO} --phenoCol=${PHENO_COL} \
    --covarColList=Sex,Age,Batch,PC1,PC2,PC3,PC4 --qCovarCol=Sex \
    --sampleIDColinphenoFile=eid --traitType=binary \
    --outputPrefix=${OUT}/saige_out --nThreads=16 --LOCO=FALSE \
    --skipVarianceRatioEstimation=TRUE --IsOverwriteVarianceRatioFile=TRUE \
    > "${OUT}/docker_run.log" 2>&1
CID=$(cat "${CID_FILE}")
echo "  container=${CID:0:12}  starting peak-RSS poller..."

# Poll peak RSS every 5 s until container stops.  `docker stats` prints
# MemUsage like "12.3GiB / 251GiB"; we extract the used side in bytes.
(
  peak=0
  while docker ps --no-trunc --format '{{.ID}}' 2>/dev/null | grep -q "${CID}"; do
      line=$(docker stats --no-stream --format '{{.MemUsage}}' "${CID}" 2>/dev/null || true)
      # Example: "4.57GiB / 251GiB"
      val=${line%% /*}
      # Convert to MiB for comparison
      num=${val%[KMGT]*}
      unit=${val##*[0-9]}
      case "$unit" in
        KiB|KB) mib=$(awk "BEGIN{print $num/1024}");;
        MiB|MB) mib=$num;;
        GiB|GB) mib=$(awk "BEGIN{print $num*1024}");;
        TiB|TB) mib=$(awk "BEGIN{print $num*1024*1024}");;
        *)      mib=0;;
      esac
      # awk max compare (peak is MiB, possibly float)
      if [ -n "$mib" ] && awk "BEGIN{exit !($mib > $peak)}"; then
        peak=$mib
      fi
      sleep 5
  done
  echo "$peak" > "${PEAK_FILE}"
) &
POLL_PID=$!

# Wait for container to exit, then wall+exit via docker wait
RDSTART=$(date +%s)
EXITCODE=$(docker wait "${CID}" 2>/dev/null)
RDEND=$(date +%s)
WALL=$((RDEND - RDSTART))
wait "${POLL_PID}" 2>/dev/null

# Capture R stdout/stderr via docker logs
docker logs "${CID}" > "${OUT}/stdout.log" 2> "${OUT}/stderr.log"
docker rm "${CID}" > /dev/null 2>&1

PEAK_MIB=$(cat "${PEAK_FILE}" 2>/dev/null)
PEAK_GB=$(awk "BEGIN{printf \"%.2f\", ${PEAK_MIB:-0}/1024}")
printf "  wall=%ds  (peak-RSS=%s GiB)  exit=%s\n" "$WALL" "$PEAK_GB" "$EXITCODE"
grep -A1 "^Tau:" "${OUT}/stdout.log" 2>/dev/null | tail -4
echo ""

echo "=== [$(date -Iseconds)] DONE ==="
