#!/bin/bash
# run_firth_precision_gate.sh BIN_DIR WORKDIR [MODE]
#
# The Firth gate of gpuPrecisionFirth=MODE (default fp32) against all-fp64:
#   1. run_precision_gate.sh on its three sets (bt_full, bt_sparse, g200k;
#      Firth on, pCutoffforFirth 0.01);
#   2. g200k at pCutoffforFirth 0.05 (both runs), for more Firth pairs;
# with SAIGE_FIRTH_STATS=1 so each log carries the fit outcome counts, then
# firth_precision_report.py on every set. Simulated data only.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(cd "${1:?usage: run_firth_precision_gate.sh BIN_DIR WORKDIR [MODE]}" && pwd)
W=$(mkdir -p "${2:?WORKDIR}" && cd "$2" && pwd)
MODE=${3:-fp32}
export SAIGE_FIRTH_STATS=1

"$HERE/run_precision_gate.sh" "$BIN" "$W" gpuPrecisionFirth=$MODE

# g200k, cutoff 0.05 on both sides
SRC=/opt/saige/logs/s2-pgen/runs/gate_v1/g200k_8_f1/gpu_bed/cfg.yaml
S=$W/g200k_c05
rm -rf "$S"
for side in fp64 test; do
  python3 - "$SRC" "$S/$side" "$side" "$MODE" <<'EOF'
import sys, os, yaml
src, outdir, side, mode = sys.argv[1:5]
c = yaml.safe_load(open(src))
c['models'] = c['models'][:8]
os.makedirs(os.path.join(outdir, 'out'), exist_ok=True)
for m in c['models']:
    m['outputFile'] = os.path.join(outdir, 'out', m['traitName'] + '.txt')
c['useGPU'] = True
c['outputFormat'] = 'text'
c['pCutoffforFirth'] = 0.05
if side == 'test':
    c['gpuPrecisionFirth'] = mode
yaml.safe_dump(c, open(os.path.join(outdir, 'cfg.yaml'), 'w'), sort_keys=False)
EOF
  "$BIN/saige-step2" "$S/$side/cfg.yaml" > "$S/$side/log.txt" 2>&1 || echo "g200k_c05 $side failed"
done
echo "== g200k_c05 (pCutoffforFirth 0.05, gpuPrecisionFirth=$MODE)"
python3 "$HERE/../tools/precision_compare.py" "$S/fp64" "$S/test" --quiet --json "$W/g200k_c05.json" | sed 's/^/   /'

echo
python3 "$HERE/firth_precision_report.py" "$W"/bt_full "$W"/bt_sparse "$W"/g200k "$S"
