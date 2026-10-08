# env.sh -- sourced by every example script. Set these for your machine (or export them first).
#
#   SAIGE_HOME  the SAIGE_cpp_260716 directory of this repository
#   CONDA_ENV   the conda environment with the build dependencies (docs/index.md, Installation)
#   CUDA_HOME   CUDA toolkit root (only needed for the GPU build)
#   WORK        where the examples write their data and results
#   CONDA_SH    conda.sh to source, when CONDA_ENV is not under <conda base>/envs/
SAIGE_HOME=${SAIGE_HOME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
CONDA_ENV=${CONDA_ENV:-$HOME/miniforge3/envs/saige-build}
CUDA_HOME=${CUDA_HOME:-/usr/local/cuda}
WORK=${WORK:-$PWD/saige_example}

# conda activate without needing `conda init` in this shell (conda's activate
# scripts read unset variables, so `set -u` is lifted around it)
case $- in *u*) _u=1; set +u;; *) _u=0;; esac
source "${CONDA_SH:-$(dirname "$(dirname "$CONDA_ENV")")/etc/profile.d/conda.sh}"
conda activate "$CONDA_ENV"
[ "$_u" = 1 ] && set -u
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

SAIGE=${SAIGE:-$SAIGE_HOME/bin/saige-gpu-cpp}     # the command-line program (step1, step2, sgs2txt)
S1=${S1:-$SAIGE_HOME/bin/saige-null}              # step-1 engine (YAML config), run by $SAIGE step1
S2=${S2:-$SAIGE_HOME/bin/saige-step2}             # step-2 engine (YAML config), run by $SAIGE step2
SGS2TXT=${SGS2TXT:-$SAIGE_HOME/bin/sgs2txt}       # .sgs -> text converter ($SAIGE sgs2txt runs it)
PLINK2=${PLINK2:-plink2}                          # only the examples need plink2
