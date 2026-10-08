#!/bin/bash
# Build step 1 with no R reachable: env -i, PATH = a fake Rscript/R that logs the call
# and fails + the conda env's bin (compiler) + /usr/bin:/bin. After the build,
# <logfile>.rcalls must not exist. Then check: readelf -d / ldd / nm -D of saige-null.
# usage: build_without_r.sh <srcdir> <logfile> [make args, e.g. NVCC=/nonexistent for the CPU build]
src=$1; log=$(realpath -m "$2"); shift 2
C=/home/francisfenglu4/miniforge3/envs/saige-build
fb=$(mktemp -d); rm -f "$log.rcalls"
for n in Rscript R; do printf '#!/bin/sh\necho "%s $*" >> %s\nexit 1\n' $n "$log.rcalls" > $fb/$n; chmod +x $fb/$n; done
cd "$src" && env -i HOME=$HOME PATH=$fb:$C/bin:/usr/bin:/bin \
  CONDA_PREFIX=$C CXX=$C/bin/x86_64-conda-linux-gnu-c++ \
  make -j8 "$@" > "$log" 2>&1
echo "make rc=$?" >> "$log"; rm -rf $fb
