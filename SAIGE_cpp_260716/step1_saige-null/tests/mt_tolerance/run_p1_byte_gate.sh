#!/usr/bin/env bash
# P=1 byte-identity regression gate: every single-phenotype output file must be
# byte-identical between <new> and <baseline> saige-null (binary/quant x CPU/GPU,
# LOCO on 5 chromosomes, sparse direct and +PCG, nrun=2, covariate_offset off,
# covariate_qr off, no covariates, no VR, categorical VR, P=1 with missing
# phenotypes, mid-scale GPU). GPU cases fail without "[gpu_matvec] tier=4".
#
#   run_p1_byte_gate.sh <new saige-null> <baseline saige-null> [--cases a,b] [--workdir DIR] [--label L] [--no-cache]
#
# Outputs under WORKDIR/p1/<label>/ (default workdir /opt/saige/logs/mt_gate);
# baseline outputs cached under WORKDIR/p1_cache by baseline md5 + config.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -lt 2 ]]; then sed -n '2,12p' "$0"; exit 2; fi
if [[ -z "${R_HOME:-}" || -z "${CONDA_PREFIX:-}" ]]; then
  set +u
  source /home/francisfenglu4/miniforge3/etc/profile.d/conda.sh
  conda activate saige-build
  set -u
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  export R_HOME="$CONDA_PREFIX/lib/R"
fi
exec /usr/bin/python3 "$HERE/p1_byte_gate.py" "$@"
