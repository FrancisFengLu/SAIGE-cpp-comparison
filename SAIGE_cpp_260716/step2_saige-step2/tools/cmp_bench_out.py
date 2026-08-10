#!/usr/bin/env python3
"""Compare a C++ saige-step2 output against the R SAIGE output for the same run.

  cmp_bench_out.py single <cpp.txt> <r.txt>
  cmp_bench_out.py region <cpp.txt> <r.txt>

Prints row counts on both sides, the number of matched keys, and the max
absolute / relative deviation of the p-values (and of -log10 p, which is the
scale people actually read).  Exit status is always 0 -- this is a report, not
an assertion.
"""
import math
import sys


def load(path, keycols, pcol="Pvalue"):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    if len(header) == 1:
        with open(path) as f:
            header = f.readline().split()
        sep = None
    else:
        sep = "\t"
    idx = {c: i for i, c in enumerate(header)}
    kc = [idx[c] for c in keycols if c in idx]
    if not kc:
        raise SystemExit("no key columns %s in %s (header=%s)" % (keycols, path, header))
    pc = None
    for cand in (pcol, "p.value", "Pvalue", "Pvalue_Burden", "p_value"):
        if cand in idx:
            pc = idx[cand]
            pcolname = cand
            break
    if pc is None:
        raise SystemExit("no p-value column in %s (header=%s)" % (path, header))
    rows = {}
    with open(path) as f:
        f.readline()
        for line in f:
            t = line.rstrip("\n").split(sep) if sep else line.split()
            if len(t) <= pc:
                continue
            try:
                v = float(t[pc])
            except ValueError:
                continue
            rows[tuple(canon(t[c]) for c in kc)] = v
    return rows, pcolname


def canon(s):
    """Normalise a key field. R prints max_MAF as '1e-04' where the C++ prints
    '0.0001'; that is pure formatting and must not split the join."""
    try:
        return repr(float(s))
    except ValueError:
        return s


def main():
    mode, cpp_path, r_path = sys.argv[1], sys.argv[2], sys.argv[3]
    if mode == "single":
        keys = ["CHR", "POS", "MarkerID"]
    else:
        keys = ["Region", "Group", "max_MAF"]
    a, pa = load(cpp_path, keys)
    b, pb = load(r_path, keys)
    common = set(a) & set(b)
    print("  cpp rows : %d   (p column %s)" % (len(a), pa))
    print("  R   rows : %d   (p column %s)" % (len(b), pb))
    print("  matched  : %d   cpp-only %d   R-only %d"
          % (len(common), len(a) - len(common), len(b) - len(common)))
    if not common:
        print("  NO OVERLAP")
        return
    if mode == "region":
        # Cauchy rows are omnibus combinations of the per-mask rows, so they
        # amplify whatever deviation the individual masks carry.  Report them
        # separately rather than letting them define the headline number.
        for lab, sel in (("per-mask", lambda k: k[1] != "Cauchy"),
                         ("Cauchy  ", lambda k: k[1] == "Cauchy")):
            ks = [k for k in common if sel(k)]
            if not ks:
                continue
            m = max(abs(a[k] - b[k]) for k in ks)
            r = max(abs(a[k] - b[k]) / max(abs(a[k]), abs(b[k]), 1e-300) for k in ks)
            print("  %s rows %5d  max|dp| %.3e  max rel %.3e" % (lab, len(ks), m, r))
    max_abs = max_rel = max_log = 0.0
    worst = None
    for k in common:
        x, y = a[k], b[k]
        d = abs(x - y)
        if d > max_abs:
            max_abs, worst = d, (k, x, y)
        if max(abs(x), abs(y)) > 0:
            max_rel = max(max_rel, d / max(abs(x), abs(y)))
        if x > 0 and y > 0:
            max_log = max(max_log, abs(math.log10(x) - math.log10(y)))
    print("  max |dp|          : %.3e" % max_abs)
    print("  max relative dp   : %.3e" % max_rel)
    print("  max |dlog10(p)|   : %.3e" % max_log)
    if worst:
        print("  worst key         : %s  cpp=%.6e  R=%.6e" % (worst[0], worst[1], worst[2]))


if __name__ == "__main__":
    main()
