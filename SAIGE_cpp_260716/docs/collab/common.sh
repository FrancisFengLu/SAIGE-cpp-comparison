# common.sh -- sourced by the collaborator test scripts (docs/collaborator_tests.md).
# Loads the build environment of the user guide (docs/examples/env.sh) and sets:
#   COLLAB    this directory
#   BIN_DIR   where build_bins.sh puts the binaries (default $SAIGE_HOME/collab_bin)
#   NULLBIN   step 1 (GPU build)          S2BIN   step 2 (GPU build; also used for the CPU path)
#   S2PHASE   step 2, CPU build with PHASE_TIMING=1 (per-stage timers, block 3b only)
COLLAB=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$COLLAB/../examples/env.sh"
BIN_DIR=${BIN_DIR:-$SAIGE_HOME/collab_bin}
NULLBIN=$BIN_DIR/saige-null
S2BIN=$BIN_DIR/saige-step2
S2PHASE=$BIN_DIR/saige-step2.phase
PYTHON=${PYTHON:-python3}
NPROC=$(nproc)
