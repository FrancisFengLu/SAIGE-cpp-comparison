#!/usr/bin/env python3
"""Compare two runs' .sgs trait files for the mtVecQuantStats A/B.

usage: mtvec_cmp_sgs.py DIR_A DIR_B [PREFIX] [--dump-stat FILE]

mtVecQuantStats only changes how the chi-square(1) upper tail is evaluated, so
the gate has two halves:

  1. BETA / SE / Tstat / var must be BIT identical -- they are written by the
     same expressions on either path, so anything but 0 here is a bug, not a
     rounding difference.
  2. the p-value may move. .sgs stores the double the printed "%.6E" parses
     back to (sgs_format.hpp), i.e. 7 significant digits, so what this reports
     for p is the PRINTED resolution: |d(-log10 p)| below 4.34e-7 shows up as
     0. The true fp64 difference is measured by tools/mtvec_tail_accuracy on
     the stat values this script can dump with --dump-stat.
"""
import sys, os, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mtfold_sgsread import read_trait

args = [a for a in sys.argv[1:] if not a.startswith("--")]
dump = None
for i, a in enumerate(sys.argv):
    if a == "--dump-stat":
        dump = sys.argv[i + 1]
A, B = args[0], args[1]
pref = args[2] if len(args) > 2 else ""

files = sorted(os.path.basename(p) for p in glob.glob(A + "/*.sgs"))
files = [f for f in files if f.startswith(pref) and ".markers." not in f]

agg = {k: [0.0, 0.0, 0] for k in ("BETA", "SE", "Tstat", "var")}   # maxabs, maxrel, ndiff
npair = 0
nbyte = 0
pdiff = 0
maxdlp = 0.0
where = ""
ncross = 0
stats = []

def as_float(col):
    """p.value column is an object array: floats plus verbatim strings."""
    out = np.empty(len(col), dtype=float)
    for i, v in enumerate(col):
        if isinstance(v, str):
            try:
                out[i] = float(v)
            except ValueError:
                out[i] = np.nan
        else:
            out[i] = v
    return out

for fn in files:
    ba = open(A + "/" + fn, 'rb').read()
    bb = open(B + "/" + fn, 'rb').read()
    if ba == bb:
        nbyte += 1
    da, _ = read_trait(A + "/" + fn)
    db, _ = read_trait(B + "/" + fn)
    n = len(da["BETA"]); npair += n
    for k in agg:
        x, y = da[k], db[k]
        d = np.abs(x - y)
        den = np.maximum(np.abs(x), np.abs(y))
        rel = np.where(den > 0, d / np.maximum(den, 1e-300), 0.0)
        agg[k][0] = max(agg[k][0], float(np.nanmax(d)))
        agg[k][1] = max(agg[k][1], float(np.nanmax(rel)))
        agg[k][2] += int(np.sum(x != y))
    pa, pb = as_float(da["p.value"]), as_float(db["p.value"])
    ne = np.where(pa != pb)[0]
    pdiff += int(ne.size)
    if ne.size:
        with np.errstate(divide='ignore'):
            l1 = -np.log10(pa[ne]); l2 = -np.log10(pb[ne])
        d = np.abs(l1 - l2)
        j = int(np.nanargmax(d))
        if d[j] > maxdlp:
            maxdlp = float(d[j])
            where = "%s row %d  p %.7E -> %.7E" % (fn, ne[j], pa[ne[j]], pb[ne[j]])
    ncross += int(np.sum((pa <= 5e-8) != (pb <= 5e-8)))
    if dump:
        stats.append(da["Tstat"] ** 2 / da["var"])

print("files %d, byte-identical %d/%d, pairs %d" % (len(files), nbyte, len(files), npair))
for k in ("BETA", "SE", "Tstat", "var"):
    print("  %-6s differing %d   max|d| %.3e   max rel %.3e"
          % (k, agg[k][2], agg[k][0], agg[k][1]))
print("  p.value (printed, 7 digits) differing %d (%.6f%%)"
      % (pdiff, 100.0 * pdiff / max(npair, 1)))
print("  max |d(-log10 p)| at printed resolution = %.3e   %s" % (maxdlp, where))
print("  pairs crossing p=5e-8: %d" % ncross)
if dump:
    np.concatenate(stats).astype("<f8").tofile(dump)
    print("  stat values written to %s (%d doubles)" % (dump, npair))
