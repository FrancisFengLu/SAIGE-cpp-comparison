#!/bin/bash
# 08_step2_sgs.sh -- the step 2 run of 04_step2_binary.sh with the binary output
# format, then conversion back to text (as you would on a separate CPU machine).
# Writes $WORK/step2_sgs/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$WORK"

# same flags as 04, plus --outputFormat sgs
$SAIGE step2 \
  --step1Dir step1_bin \
  --plinkFile data/geno \
  --minMAF 0 \
  --minMAC 1 \
  --LOCO=FALSE \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth 0.01 \
  --is_noadjCov=TRUE \
  --impute_method best_guess \
  --is_fastTest=FALSE \
  --nThreads 8 \
  --useGPU \
  --outputFormat sgs \
  --outDir step2_sgs > step2_sgs.log 2>&1
ls -l step2_sgs        # b1.txt.sgs ... b4.txt.sgs + the shared b1.txt.markers.sgs

# Convert. The .sgs files and the shared markers file can be copied anywhere;
# -m names the markers file, -o the text file to write.
mkdir -p step2_sgs/moved step2_sgs/text
cp step2_sgs/*.sgs step2_sgs/moved/
for t in b1 b2 b3 b4; do
  $SAIGE sgs2txt -m step2_sgs/moved/b1.txt.markers.sgs -o step2_sgs/text/$t.txt step2_sgs/moved/$t.txt.sgs
done
# The text is identical to the text writer's output from 04_step2_binary.sh
for t in b1 b2 b3 b4; do cmp step2_sgs/text/$t.txt step2_bin/$t.txt && echo "$t identical"; done
