#!/usr/bin/env bash
# Multi-trait tolerance gate for saige-null (step 1).
#
#   run_gate.sh <saige-null> [--gpu|--cpu] [--small] [--cases a,b,c] [--workdir DIR]
#               [--lockstep] [--nthreads N] [--solo-vs-solo] [--no-cache] [--label NAME]
#               [--fit KEY=VALUE ...]
#
# --fit sets a fit.* key on the MULTI run only (the solo runs stay the plain
# P=1 reference), e.g. --fit mask_missing=true --fit mask_min_coverage=0.0.
#
# For every case (cases/*.yaml): each trait alone (P=1; cached under
# WORKDIR/ref_cache keyed by binary md5 + rendered config + input md5s), then all
# traits in one run; coverage checks; compare.py for every trait; a case x trait
# table at the end (also WORKDIR/runs/<label>/summary.txt). Exit 0 iff all pass.
# Default workdir /opt/saige/logs/mt_gate, default device CPU, default nthreads 8.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -lt 1 ]]; then sed -n '2,16p' "$0"; exit 2; fi
if [[ -z "${R_HOME:-}" || -z "${CONDA_PREFIX:-}" ]]; then
  # saige-null embeds R (RNG streams) and links the conda env's libraries
  set +u
  source /home/francisfenglu4/miniforge3/etc/profile.d/conda.sh
  conda activate saige-build
  set -u
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  export R_HOME="$CONDA_PREFIX/lib/R"
fi
exec /usr/bin/python3 "$HERE/gate.py" "$@"
