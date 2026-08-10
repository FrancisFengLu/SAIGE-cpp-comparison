#!/usr/bin/env python3
"""Assertions over saige-step2 single-variant output files.

  cmp_pval.py differ  A B   -> exit 0 if the markers common to A and B have
                              MEANINGFULLY different p-values
  cmp_pval.py same    A B   -> exit 0 if the common markers have identical
                              p-values (bit-for-bit on the printed text)
  cmp_pval.py onlychr C F   -> exit 0 if every row of F is on chromosome C
"""
import sys


def read(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    key_cols = [idx[c] for c in ("CHR", "POS", "MarkerID") if c in idx]
    pcol = idx.get("p.value")
    if pcol is None:
        raise SystemExit("no p.value column in " + path)
    rows = {}
    with open(path) as f:
        f.readline()
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) <= pcol:
                continue
            rows[tuple(t[c] for c in key_cols)] = (t[pcol], t[idx["CHR"]] if "CHR" in idx else "")
    return rows


def main():
    mode = sys.argv[1]
    if mode == "onlychr":
        chrom = sys.argv[2].lstrip("chrCHR")
        rows = read(sys.argv[3])
        if not rows:
            print("  (onlychr: output file has no rows)")
            return 1
        bad = [k for k, v in rows.items() if v[1].lstrip("chrCHR") != chrom]
        if bad:
            print("  %d of %d rows are not on chromosome %s (e.g. %s)"
                  % (len(bad), len(rows), chrom, bad[0]))
            return 1
        print("  all %d rows on chromosome %s" % (len(rows), chrom))
        return 0

    a, b = read(sys.argv[2]), read(sys.argv[3])
    common = set(a) & set(b)
    if not common:
        print("  no markers in common between the two outputs")
        return 1
    ndiff = sum(1 for k in common if a[k][0] != b[k][0])
    if mode == "differ":
        print("  %d of %d common markers have different p-values" % (ndiff, len(common)))
        # Require essentially all of them to move, not just rounding noise.
        return 0 if ndiff > 0.9 * len(common) else 1
    if mode == "same":
        print("  %d of %d common markers differ" % (ndiff, len(common)))
        return 0 if ndiff == 0 else 1
    raise SystemExit("unknown mode " + mode)


if __name__ == "__main__":
    sys.exit(main())
