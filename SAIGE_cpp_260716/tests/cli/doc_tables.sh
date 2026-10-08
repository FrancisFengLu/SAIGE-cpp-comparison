#!/bin/bash
# doc_tables.sh [--write] -- the flag tables of docs/step1.md and docs/step2.md
# against `saige-gpu-cpp stepN --help-markdown`. Without --write: report whether
# they match (exit 1 if not). With --write: replace the text between the
# "<!-- flag tables ... -->" and "<!-- end flag tables -->" markers.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
H=$(cd "$HERE/../.." && pwd)
SAIGE=${SAIGE:-$H/bin/saige-gpu-cpp}
rc=0
for st in step1 step2; do
  doc=$H/docs/$st.md
  gen=$(mktemp); "$SAIGE" $st --help-markdown > "$gen"
  cur=$(mktemp); awk '/<!-- flag tables/{f=1;next} /<!-- end flag tables -->/{f=0} f' "$doc" > "$cur"
  if [ "${1:-}" = --write ]; then
    awk -v g="$gen" '/<!-- flag tables/{print; while ((getline l < g) > 0) print l; skip=1; next}
                     /<!-- end flag tables -->/{skip=0} !skip' "$doc" > "$doc.tmp" && mv "$doc.tmp" "$doc"
    echo "$st.md: tables written"
  elif cmp -s "$gen" "$cur"; then echo "$st.md: tables match --help-markdown"
  else echo "$st.md: tables differ from --help-markdown"; diff "$cur" "$gen" | head -20; rc=1; fi
  rm -f "$gen" "$cur"
done
exit $rc
