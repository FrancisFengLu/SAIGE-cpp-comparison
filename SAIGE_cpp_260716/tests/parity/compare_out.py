#!/usr/bin/env python3
"""Column-wise comparison of two SAIGE Step-2 output tables (R vs C++).

Not a byte diff: the two writers disagree on float formatting, so every
numeric column is parsed and compared as a number, and p-values are compared
on the log10 scale as well (a relative difference on a p-value of 1e-300 is
meaningless; the exponent is what matters).

    compare_out.py <R output> <C++ output> [--rtol 1e-6] [--json out.json]

Exit status: 0 if every column agrees within tolerance and both sides carry
the same rows, 1 otherwise.
"""
import argparse
import json
import math
import sys

import numpy as np

# Candidate key columns, in priority order.  The first set fully contained in
# the header is used to join the two tables.
KEY_SETS = [
    ["Region", "Group", "max_MAF"],                       # region test
    ["CHR", "POS", "MarkerID", "Allele1", "Allele2"],      # single variant
    ["MarkerID"],
    ["Region", "Group"],
]

# Columns whose "value" is a p-value: compared on log10 as well as relatively.
PVAL_COLS_SUFFIX = ("p.value", "p.value.NA", "Pvalue", "Pvalue_Burden",
                    "Pvalue_SKAT", "Pvalue_cond", "Pvalue_Burden_cond",
                    "Pvalue_SKAT_cond")


def read_table(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    ncol = len(header)
    rows = [r for r in rows if len(r) == ncol]
    cols = {h: [r[i] for r in rows] for i, h in enumerate(header)}
    return header, cols, len(rows)


def pick_key(header):
    for ks in KEY_SETS:
        if all(k in header for k in ks):
            return ks
    return None


def as_float(vals):
    """Parse a column as float. Returns (array, ok) -- ok=False if not numeric."""
    out = np.empty(len(vals), dtype=float)
    for i, v in enumerate(vals):
        s = v.strip()
        if s in ("", "NA", "NaN", "nan", "None", "."):
            out[i] = math.nan
            continue
        try:
            out[i] = float(s)
        except ValueError:
            return None, False
    return out, True


def compare(r_path, c_path, rtol=1e-6, atol=0.0):
    rh, rc, rn = read_table(r_path)
    ch, cc, cn = read_table(c_path)

    rep = {"r_file": r_path, "cpp_file": c_path,
           "r_rows": rn, "cpp_rows": cn,
           "rtol": rtol, "atol": atol,
           "header_only_in_r": [h for h in rh if h not in ch],
           "header_only_in_cpp": [h for h in ch if h not in rh],
           "columns": [], "verdict": "IDENTICAL"}

    shared = [h for h in rh if h in ch]
    key = pick_key(shared)
    if key is None:
        rep["verdict"] = "NO_KEY"
        rep["error"] = ("no usable key column set in the shared header: %s"
                        % ",".join(shared))
        return rep
    rep["key"] = key

    def mk(cols, n):
        return ["\t".join(cols[k][i] for k in key) for i in range(n)]

    rkeys, ckeys = mk(rc, rn), mk(cc, cn)
    rindex = {k: i for i, k in enumerate(rkeys)}
    cindex = {k: i for i, k in enumerate(ckeys)}
    common = [k for k in rkeys if k in cindex]
    rep["rows_common"] = len(common)
    rep["rows_only_in_r"] = [k for k in rkeys if k not in cindex]
    rep["rows_only_in_cpp"] = [k for k in ckeys if k not in rindex]
    rep["n_only_in_r"] = len(rep["rows_only_in_r"])
    rep["n_only_in_cpp"] = len(rep["rows_only_in_cpp"])
    # keep the report small
    rep["rows_only_in_r"] = rep["rows_only_in_r"][:20]
    rep["rows_only_in_cpp"] = rep["rows_only_in_cpp"][:20]
    if len(set(rkeys)) != rn or len(set(ckeys)) != cn:
        rep["warn_duplicate_keys"] = True

    ri = np.array([rindex[k] for k in common], dtype=int)
    ci = np.array([cindex[k] for k in common], dtype=int)

    worst = 0.0
    for h in shared:
        if h in key:
            continue
        rv = [rc[h][i] for i in ri]
        cv = [cc[h][i] for i in ci]
        ra, rok = as_float(rv)
        ca, cok = as_float(cv)
        entry = {"column": h}
        if rok and cok:
            entry["type"] = "numeric"
            both_nan = np.isnan(ra) & np.isnan(ca)
            one_nan = np.isnan(ra) ^ np.isnan(ca)
            m = ~(np.isnan(ra) | np.isnan(ca))
            entry["n_nan_mismatch"] = int(one_nan.sum())
            entry["n_both_nan"] = int(both_nan.sum())
            if m.sum():
                d = np.abs(ra[m] - ca[m])
                den = np.maximum(np.abs(ra[m]), np.abs(ca[m]))
                rel = np.where(den > 0, d / np.where(den > 0, den, 1.0), 0.0)
                entry["max_abs"] = float(d.max())
                entry["max_rel"] = float(rel.max())
                bad = (d > atol) & (rel > rtol)
                entry["n_differing"] = int(bad.sum())
                if bad.any():
                    j = np.argmax(rel)
                    entry["worst_row"] = common[int(np.flatnonzero(m)[j])]
                    entry["worst_r"] = float(ra[m][j])
                    entry["worst_cpp"] = float(ca[m][j])
                if h in PVAL_COLS_SUFFIX:
                    pm = m & (ra > 0) & (ca > 0)
                    if pm.sum():
                        ld = np.abs(np.log10(ra[pm]) - np.log10(ca[pm]))
                        entry["max_abs_log10_p"] = float(ld.max())
            else:
                entry["max_abs"] = 0.0
                entry["max_rel"] = 0.0
                entry["n_differing"] = 0
            if entry["n_differing"] or entry["n_nan_mismatch"]:
                worst = max(worst, entry.get("max_rel", 0.0))
        else:
            entry["type"] = "string"
            neq = sum(1 for a, b in zip(rv, cv) if a.strip() != b.strip())
            entry["n_differing"] = neq
            if neq:
                for a, b, k in zip(rv, cv, common):
                    if a.strip() != b.strip():
                        entry["worst_row"] = k
                        entry["worst_r"] = a
                        entry["worst_cpp"] = b
                        break
        rep["columns"].append(entry)

    diff_cols = [c for c in rep["columns"]
                 if c.get("n_differing") or c.get("n_nan_mismatch")]
    rep["n_columns_identical"] = len(rep["columns"]) - len(diff_cols)
    rep["n_columns_differing"] = len(diff_cols)
    rep["max_rel_overall"] = worst
    if rep["header_only_in_r"] or rep["header_only_in_cpp"]:
        rep["verdict"] = "HEADER_MISMATCH"
    elif rep["n_only_in_r"] or rep["n_only_in_cpp"]:
        rep["verdict"] = "ROW_SET_MISMATCH"
    elif diff_cols:
        rep["verdict"] = "VALUES_DIFFER"
    return rep


def print_report(rep):
    print("R   : %s  (%d rows)" % (rep["r_file"], rep["r_rows"]))
    print("cpp : %s  (%d rows)" % (rep["cpp_file"], rep["cpp_rows"]))
    if "error" in rep:
        print("ERROR: " + rep["error"])
        return
    print("key : %s   common rows: %d   only-in-R: %d   only-in-cpp: %d"
          % ("+".join(rep["key"]), rep["rows_common"],
             rep["n_only_in_r"], rep["n_only_in_cpp"]))
    if rep["header_only_in_r"]:
        print("columns only in R  : %s" % ", ".join(rep["header_only_in_r"]))
    if rep["header_only_in_cpp"]:
        print("columns only in cpp: %s" % ", ".join(rep["header_only_in_cpp"]))
    if rep["n_only_in_r"]:
        print("  e.g. only in R  : %s" % "; ".join(rep["rows_only_in_r"][:3]))
    if rep["n_only_in_cpp"]:
        print("  e.g. only in cpp: %s" % "; ".join(rep["rows_only_in_cpp"][:3]))
    print()
    print("%-22s %-8s %10s %12s %12s %9s" %
          ("column", "type", "n_diff", "max_abs", "max_rel", "dlog10p"))
    print("-" * 78)
    for c in rep["columns"]:
        star = " " if not (c.get("n_differing") or c.get("n_nan_mismatch")) else "*"
        print("%s%-21s %-8s %10s %12s %12s %9s" % (
            star, c["column"], c["type"],
            c.get("n_differing", 0),
            ("%.3e" % c["max_abs"]) if "max_abs" in c else "-",
            ("%.3e" % c["max_rel"]) if "max_rel" in c else "-",
            ("%.2e" % c["max_abs_log10_p"]) if "max_abs_log10_p" in c else "-"))
        if star == "*" and "worst_row" in c:
            print("     worst: %s  R=%s  cpp=%s"
                  % (c["worst_row"], c["worst_r"], c["worst_cpp"]))
    print("-" * 78)
    print("identical columns: %d    differing columns: %d    max rel: %.3e"
          % (rep["n_columns_identical"], rep["n_columns_differing"],
             rep["max_rel_overall"]))
    print("VERDICT: %s" % rep["verdict"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r_file")
    ap.add_argument("cpp_file")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--json")
    a = ap.parse_args()
    rep = compare(a.r_file, a.cpp_file, a.rtol, a.atol)
    print_report(rep)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=2)
    sys.exit(0 if rep["verdict"] == "IDENTICAL" else 1)


if __name__ == "__main__":
    main()
