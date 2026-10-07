#!/bin/bash
# env_record.sh -- record the machine and software used for the tests into $OUT_ROOT/env/.
#   environment.txt   collected automatically (hardware, OS, GPU, CUDA, filesystems, tools)
#   manual.txt        a few lines for you to fill in (only if it does not exist yet)
# Optional: GENO (the step-2 genotype prefix, to record its filesystem), RSCRIPT.
set -uo pipefail
source "$(dirname "$0")/common.sh"
: "${OUT_ROOT:?}"
E=$OUT_ROOT/env; mkdir -p "$E"
sec() { echo; echo "## $1"; }
{
  echo "# environment record, $(date -Is)"
  sec OS;      grep -E '^(PRETTY_NAME|VERSION_ID)=' /etc/os-release; uname -r
  sec CPU;     lscpu | grep -E '^(Model name|Socket\(s\)|Core\(s\) per socket|Thread\(s\) per core|CPU\(s\)|NUMA node\(s\)|L3 cache)'
               echo "nproc (usable by this process): $(nproc)"
               echo "flags: $(grep -m1 -o -w -E 'avx2|avx512f|avx512vl|fma|bmi2' /proc/cpuinfo | sort -u | tr '\n' ' ')"
               [ -n "${SLURM_JOB_ID:-}" ] && echo "SLURM: cpus_on_node=${SLURM_CPUS_ON_NODE:-} cpus_per_task=${SLURM_CPUS_PER_TASK:-} mem=${SLURM_MEM_PER_NODE:-}"
  sec Memory;  free -g | head -2
  sec GPU
  if command -v nvidia-smi > /dev/null; then
    nvidia-smi --query-gpu=index,name,memory.total,driver_version,compute_cap,persistence_mode,pcie.link.gen.max,pcie.link.width.max --format=csv
    nvidia-smi | grep -m1 -o 'CUDA Version: [0-9.]*'
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
  else echo "nvidia-smi not found"; fi
  sec CUDA;    "$CUDA_HOME/bin/nvcc" --version 2>/dev/null | tail -2 || echo "nvcc not found under $CUDA_HOME"
  sec Filesystems
  for p in "$OUT_ROOT" ${GENO:+"$(dirname "$GENO")"}; do
    echo "$p: type=$(stat -f -c %T "$p") $(df -hP "$p" | awk 'NR==2{print "size=" $2 " free=" $4}')"
  done
  lsblk -d -o NAME,ROTA,TYPE,SIZE,MODEL 2>/dev/null | head -20
  sec Tools
  echo "GNU time: $(/usr/bin/time --version 2>&1 | head -1)"
  echo "python: $($PYTHON --version 2>&1)"
  echo "plink2: $("$PLINK2" --version 2>/dev/null | head -1 || echo not found)"
  echo "vmtouch: $(command -v vmtouch || echo no)   fincore: $(command -v fincore || echo no)"
  echo "passwordless sudo (for drop_caches): $(sudo -n true 2>/dev/null && echo yes || echo no)"
  echo "R: $(${RSCRIPT:-Rscript} --version 2>&1 | head -1)"
  echo "R SAIGE: $(${RSCRIPT:-Rscript} -e 'cat(as.character(packageVersion("SAIGE")))' 2>/dev/null | tail -1 || echo 'not found')"
  sec "SAIGE C++ build"
  echo "source commit: $(git -C "$SAIGE_HOME" rev-parse HEAD 2>/dev/null)"
  cat "$BIN_DIR/BUILD_INFO.txt" 2>/dev/null || echo "no $BIN_DIR/BUILD_INFO.txt (run build_bins.sh)"
} > "$E/environment.txt" 2>&1
[ -e "$E/manual.txt" ] || cat > "$E/manual.txt" <<'TXT'
# Fill in, then keep this file in $OUT_ROOT/env/ (collect_results.sh includes it).
site / cluster name:
machine type (cloud instance type or node model):
exclusive node during the runs? (yes/no; other jobs on the node?):
disk holding the genotype files (local NVMe / SSD / HDD / network: Lustre, GPFS, NFS, ...):
cohort size N in the step-2 genotype file:
chromosome tested, and number of markers in it:
number of binary traits, and their case fractions (range is enough):
how the sparse GRM was made (step1_models.sh / FastSparseGRM / other):
anything unusual:
TXT
cat "$E/environment.txt"
echo; echo "now fill in $E/manual.txt"
