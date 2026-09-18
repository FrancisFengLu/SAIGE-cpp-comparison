#!/usr/bin/env bash
# fit.fused_variance_ratio A/B gate: one binary, flag off vs on, every output
# file compared. See ab_gate.py for what each case proves.
#
#   run_ab_gate.sh <saige-null> [--cases a,b] [--outdir DIR] [--workdir DIR]
#
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -lt 1 ]]; then sed -n '2,6p' "$0"; exit 2; fi
if [[ -z "${R_HOME:-}" || -z "${CONDA_PREFIX:-}" ]]; then
  set +u
  source /home/francisfenglu4/miniforge3/etc/profile.d/conda.sh
  conda activate saige-build
  set -u
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  export R_HOME="$CONDA_PREFIX/lib/R"
fi
exec /usr/bin/python3 "$HERE/ab_gate.py" "$@"
