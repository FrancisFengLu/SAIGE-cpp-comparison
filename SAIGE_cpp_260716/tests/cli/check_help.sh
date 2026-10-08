#!/bin/bash
# check_help.sh [R_EXTDATA_DIR] -- --help completeness.
#  1. every flag line of `step1 --help` / `step2 --help` shows a [default: ...]
#  2. every make_option() flag of R SAIGE's step1_fitNULLGLMM.R / step2_SPAtests.R
#     (R_EXTDATA_DIR, the extdata/ of the SAIGE R package source) is either a
#     flag of the same name or in the refused list -- and a refused one stops the
#     run with "is an R SAIGE flag that ... does not support".
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SAIGE=${SAIGE:-$(cd "$HERE/../.." && pwd)/bin/saige-gpu-cpp}
R=${1:-}
nfail=0
for st in step1 step2; do
  H=$("$SAIGE" $st --help)
  # the flag list ends where the refused R flags start; a long flag name puts its
  # help (and default) on the next line
  F=$(sed '/^R SAIGE flags that are refused/,$d' <<<"$H")
  nflag=$(grep -cE '^  --' <<<"$F")
  nodef2=$(awk '/^  --/ { if ($0 ~ /\[default:/) next; if ((getline n) > 0 && n ~ /^      +.*\[default:/) next; print }' <<<"$F")
  if [ -n "$nodef2" ]; then echo "$st: flags without a default:"; echo "$nodef2"; nfail=$((nfail+1)); fi
  echo "$st: $nflag flags listed, each with [default: ...]"
  # defaults come from the engine (--print-defaults); none may be missing
  if grep -q 'engine default (' <<<"$F"; then
    echo "$st: defaults the engine did not report:"; grep 'engine default (' <<<"$F"; nfail=$((nfail+1))
  else echo "$st: every engine-key default was read from the engine (--print-defaults)"; fi
  [ -z "$R" ] && continue
  rs=$([ $st = step1 ] && echo step1_fitNULLGLMM.R || echo step2_SPAtests.R)
  missing=0; nr=0; nref=0
  for f in $(grep -oE 'make_option\("--[A-Za-z_.0-9]+"' "$R/$rs" | sed -E 's/.*"--([^"]+)"/\1/'); do
    nr=$((nr+1))
    if grep -qE -- "^  --$f(\[| |$)" <<<"$H"; then continue; fi
    if grep -qE -- "^  --$f: " <<<"$H"; then
      nref=$((nref+1))
      out=$("$SAIGE" $st --$f=__x__ 2>&1)
      grep -q "is an R SAIGE flag that" <<<"$out" || { echo "  --$f refused in --help but not by the parser: $out"; missing=$((missing+1)); }
      continue
    fi
    echo "  $st: R flag --$f is neither supported nor refused"; missing=$((missing+1))
  done
  echo "$st: $nr R flags: $((nr-nref-missing)) supported, $nref refused, $missing unaccounted"
  [ "$missing" = 0 ] || nfail=$((nfail+1))
done
[ "$nfail" = 0 ] && echo "OK" || { echo "FAILED"; exit 1; }
