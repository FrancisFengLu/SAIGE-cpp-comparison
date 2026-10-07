#!/bin/bash
# 08_step2_sgs.sh -- the step 2 run of 04_step2_binary.sh with the binary output
# format, then conversion back to text (as you would on a separate CPU machine).
# Writes $WORK/step2_sgs/.
set -euo pipefail
source "$(dirname "$0")/env.sh"
O=$WORK/step2_sgs
mkdir -p "$O/out"
# same config as 04, different output paths and format
sed -e "s#$WORK/step2_bin/out/#$O/out/#" -e "s#^outputFormat: text#outputFormat: sgs#" \
    "$WORK/step2_bin/step2.yaml" > "$O/step2.yaml"
"$S2" "$O/step2.yaml" > "$O/step2.log" 2>&1
ls -l "$O/out"

# Convert. The .sgs files and the shared markers file can be copied anywhere;
# -m names the markers file, -o the text file to write.
mkdir -p "$O/moved" "$O/text"
cp "$O"/out/*.sgs "$O/moved/"
for t in b1 b2 b3 b4; do
  "$SGS2TXT" -m "$O/moved/b1.txt.markers.sgs" -o "$O/text/$t.txt" "$O/moved/$t.txt.sgs"
done
# The text is identical to the text writer's output from 04_step2_binary.sh
for t in b1 b2 b3 b4; do cmp "$O/text/$t.txt" "$WORK/step2_bin/out/$t.txt" && echo "$t identical"; done
