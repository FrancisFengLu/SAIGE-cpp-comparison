#!/usr/bin/env python3
"""Separate setup/I-O from compute for the C++ runs, using saige-step2's own
[TIMING] instrumentation, and for the R runs using the elapsed time that
SPAGMMATtest/fitNULLGLMM report versus the process wall clock.

  bench_io_split.py <out_dir>

C++ [TIMING] phases:
  10_yaml_parsed / 20_null_model_loaded / 40_geno_reader_ready  -> setup + I/O
  60_after_main_loop                                            -> compute
R:
  process wall clock - the BENCH_*_ELAPSED_S the driver script prints
  -> R interpreter start-up + library(SAIGE) load, i.e. fixed overhead
"""
import glob
import os
import re
import sys


def wall_of(timefile):
    try:
        txt = open(timefile, errors="replace").read()
    except OSError:
        return None
    m = re.search(r"Elapsed \(wall clock\).*?:\s*([0-9:.]+)", txt)
    if not m:
        return None
    p = [float(x) for x in m.group(1).split(":")]
    while len(p) < 3:
        p.insert(0, 0.0)
    return p[0] * 3600 + p[1] * 60 + p[2]


def main():
    out = sys.argv[1]
    print("%-34s %10s %10s %10s %8s" % ("run", "wall_s", "setup_s", "compute_s", "setup%"))
    for log in sorted(glob.glob(os.path.join(out, "raw", "*.log"))):
        tag = os.path.basename(log)[:-4]
        if tag.startswith("prep"):
            continue
        w = wall_of(os.path.join(out, "raw", tag + ".time"))
        if w is None:
            continue
        txt = open(log, errors="replace").read()
        # saige-step2 prints its [TIMING] breakdown on stderr, which the harness
        # captures into the .time file alongside /usr/bin/time -v's own output.
        txt += open(os.path.join(out, "raw", tag + ".time"), errors="replace").read()
        ph = dict(re.findall(r"\[TIMING\] (\S+)\s+\+\S+\s+total=([0-9.e+-]+)s", txt))
        if "60_after_main_loop" in ph and "50_before_main_loop" in ph:
            setup = float(ph["50_before_main_loop"])
            compute = float(ph["60_after_main_loop"]) - setup
        else:
            m = re.search(r"BENCH_R_STEP[12]_ELAPSED_S\s*([0-9.]+)", txt)
            if not m:
                continue
            compute = float(m.group(1))
            setup = w - compute
        print("%-34s %10.2f %10.2f %10.2f %7.1f%%"
              % (tag, w, setup, compute, 100.0 * setup / w if w else 0))


if __name__ == "__main__":
    main()
