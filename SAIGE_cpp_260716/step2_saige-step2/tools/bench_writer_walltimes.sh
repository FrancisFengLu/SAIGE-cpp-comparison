#!/bin/bash
# Pull the real wall times out of the /usr/bin/time -v files (both bench2.sh and
# bench3.sh printed min_wall=0s: their awk split the "Elapsed (wall clock) time
# (h:mm:ss or m:ss): M:SS" line on "clock) " and then on ":", which lands on the
# wrong field). Everything else in those runs is unaffected.
set -u
sec() { awk -F': ' '/Elapsed .wall/{print $NF}' "$1" \
        | awk -F: '{if(NF==2)print $1*60+$2; else print $1*3600+$2*60+$3}'; }
best() { local b=999999 s; for f in "$@"; do [ -f "$f" ] || continue
           s=$(sec "$f"); b=$(echo "$s $b" | awk '{print ($1<$2)?$1:$2}'); done
         echo "$b"; }
row() { # row <label> <dir> <tag>
  local files=( "$2/$3".r*.time )
  printf "  %-12s %8.1f s   (n=%d)\n" "$1" "$(best "${files[@]}")" "${#files[@]}"
}
echo "=== before: bench2 (the run in the task's table) ==="
for P in 1 8 32 128; do
  printf "P=%-4s" $P; row "gpu fp64" /opt/saige/logs/gpu_step2/bench2 g64_P$P
done
echo
echo "=== after: bench3 ==="
for P in 1 8 32 128; do
  printf "P=%-4s\n" $P
  row "gpu text" /opt/saige/logs/gpu_step2/writer/bench3 gtxt_P$P
  row "gpu sgs"  /opt/saige/logs/gpu_step2/writer/bench3 gsgs_P$P
done
