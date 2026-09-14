#!/usr/bin/env python3
"""Compare two saige-step2 single-variant output files field by field.

The multi-trait batch kernel is algebraically equal to the scalar score test
but sums over all N samples instead of over the carriers only, so the two are
allowed to differ in the last bits. `cmp` is still the primary gate -- this
tool exists for the case where cmp fails, to answer the only question that
matters then: is every difference a last-digit rounding artefact of the printed
6-significant-figure format, or did something actually change?

  cmp_assoc.py <A> <B> [--ptol 0] [--rtol 1e-10] [--show 5]

Exit 0 when: the two files have the same rows in the same order, the number of
differing p-value STRINGS is <= --ptol, and every numeric column agrees to
within --rtol relative. Prints a per-column summary either way.
"""

import sys

NUMERIC = ("AC_Allele2", "AF_Allele2", "MissingRate", "imputationInfo",
           "BETA", "SE", "Tstat", "var", "AF_case", "AF_ctrl",
           "N_case", "N_ctrl", "N", "N_case_hom", "N_case_het",
           "N_ctrl_hom", "N_ctrl_het", "BETA_c", "SE_c", "Tstat_c", "var_c")
STRINGY = ("p.value", "p.value.NA", "p.value_c", "p.value.NA_c", "Is.SPA")
KEY = ("CHR", "POS", "MarkerID", "Allele1", "Allele2")


def load(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]
    return header, rows


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    args = sys.argv[3:]

    def opt(name, default):
        return float(args[args.index(name) + 1]) if name in args else default

    ptol = int(opt("--ptol", 0))
    rtol = opt("--rtol", 1e-10)
    show = int(opt("--show", 5))

    ha, ra = load(a_path)
    hb, rb = load(b_path)
    bad = 0

    if ha != hb:
        print(f"  HEADER differs:\n    A {ha}\n    B {hb}")
        return 1
    if len(ra) != len(rb):
        print(f"  ROW COUNT differs: A {len(ra)}, B {len(rb)}")
        return 1

    idx = {name: i for i, name in enumerate(ha)}
    keycols = [idx[c] for c in KEY if c in idx]
    for n, (x, y) in enumerate(zip(ra, rb)):
        if [x[c] for c in keycols] != [y[c] for c in keycols]:
            print(f"  ROW {n} identifies a different marker:")
            print(f"    A {[x[c] for c in keycols]}")
            print(f"    B {[y[c] for c in keycols]}")
            return 1

    for name in ha:
        c = idx[name]
        if name in KEY:
            continue
        if name in STRINGY:
            diffs = [(n, x[c], y[c]) for n, (x, y) in enumerate(zip(ra, rb)) if x[c] != y[c]]
            if diffs:
                limit = ptol if name.startswith("p.value") else 0
                verdict = "ok" if len(diffs) <= limit else "FAIL"
                if verdict == "FAIL":
                    bad += 1
                print(f"  {name:14s} {len(diffs)} of {len(ra)} strings differ  [{verdict}]")
                for n, xv, yv in diffs[:show]:
                    print(f"      row {n}: {xv} vs {yv}")
            continue
        if name not in NUMERIC:
            continue
        worst, worst_at = 0.0, None
        ndiff = 0
        for n, (x, y) in enumerate(zip(ra, rb)):
            if x[c] == y[c]:
                continue
            ndiff += 1
            try:
                fx, fy = float(x[c]), float(y[c])
            except ValueError:
                worst, worst_at = float("inf"), (n, x[c], y[c])
                break
            d = abs(fx - fy) / max(abs(fx), abs(fy), 1e-300)
            if d > worst:
                worst, worst_at = d, (n, x[c], y[c])
        if ndiff:
            verdict = "ok" if worst <= rtol else "FAIL"
            if verdict == "FAIL":
                bad += 1
            print(f"  {name:14s} {ndiff} of {len(ra)} differ, max rel {worst:.3e}  [{verdict}]")
            if worst_at:
                print(f"      worst row {worst_at[0]}: {worst_at[1]} vs {worst_at[2]}")

    if bad == 0:
        print(f"  all columns within tolerance (ptol={ptol}, rtol={rtol:g})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
