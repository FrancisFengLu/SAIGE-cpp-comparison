#!/usr/bin/env bash
# Byte-compare the r_rng port against R for several seeds:
#   first N draws of runif() (= unif_rand) and of rbinom(N, 1, 0.5), plus the
#   continuous-stream case the trace estimator uses (one set.seed, then
#   rbinom(n,1,0.5) called repeatedly).
# Needs Rscript on PATH (e.g. conda activate saige-build) and a C++17 compiler.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP="${1:-$(mktemp -d)}"; mkdir -p "$TMP"
N="${N:-1000000}"
SEEDS="${SEEDS:-10 200 1 0 -1 123456789 2147483647 -2147483647}"
g++ -std=c++17 -O2 -o "$TMP/r_rng_dump" "$HERE/r_rng_dump.cpp" "$HERE/../../r_rng.cpp"
fail=0
for s in $SEEDS; do
  for kind in unif rbinom; do
    "$TMP/r_rng_dump" "$s" "$N" "$kind" > "$TMP/cpp_${s}_${kind}.bin"
    Rscript --vanilla -e "set.seed($s); x <- if ('$kind'=='unif') runif($N) else as.double(rbinom($N,1,0.5)); writeBin(x, '$TMP/r_${s}_${kind}.bin', size=8, endian='little')"
    if cmp -s "$TMP/cpp_${s}_${kind}.bin" "$TMP/r_${s}_${kind}.bin"; then
      echo "PASS seed=$s $kind n=$N"
    else
      echo "FAIL seed=$s $kind n=$N"; fail=1
    fi
  done
done
# Continuous stream across calls: R draws 30 probe vectors of length 1000 after
# one set.seed; the port must continue the same stream across calls too.
Rscript --vanilla -e "set.seed(10); x <- unlist(lapply(1:30, function(i) as.double(rbinom(1000,1,0.5)))); writeBin(x, '$TMP/r_stream.bin', size=8, endian='little')"
"$TMP/r_rng_dump" 10 30000 rbinom > "$TMP/cpp_stream.bin"
if cmp -s "$TMP/cpp_stream.bin" "$TMP/r_stream.bin"; then echo "PASS stream seed=10 30x1000"; else echo "FAIL stream"; fail=1; fi
exit $fail
