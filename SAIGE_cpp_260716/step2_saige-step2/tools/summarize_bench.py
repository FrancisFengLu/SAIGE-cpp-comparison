#!/usr/bin/env python3
"""Turn bench_cpp_vs_r.sh's timings.csv into the markdown tables for the report.

  summarize_bench.py <timings.csv>

Warm-up rows (warmup=1) are printed in their own table and excluded from every
median / speedup figure.  Spread is reported as min-max of the measured reps,
not a standard deviation -- with 3 reps an SD is not meaningful.
"""
import csv
import sys
from collections import defaultdict


def med(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def fmt(x):
    return "%.1f" % x if x >= 10 else "%.2f" % x


def main():
    rows = list(csv.DictReader(open(sys.argv[1])))
    for r in rows:
        for k in ("wall_s", "user_s", "sys_s"):
            try:
                r[k] = float(r[k])
            except ValueError:
                r[k] = float("nan")
        r["max_rss_kb"] = int(r["max_rss_kb"]) if r["max_rss_kb"].isdigit() else 0

    cells = defaultdict(list)
    warm = []
    for r in rows:
        key = (r["phase"], r["blas"], r["stage"], r["engine"], int(r["threads"]))
        (warm if r["warmup"] == "1" else cells[key]).append(r)

    agg = {}
    for k, v in cells.items():
        agg[k] = dict(
            n=len(v),
            wall=med([x["wall_s"] for x in v]),
            wall_lo=min(x["wall_s"] for x in v),
            wall_hi=max(x["wall_s"] for x in v),
            user=med([x["user_s"] for x in v]),
            sys=med([x["sys_s"] for x in v]),
            rss=max(x["max_rss_kb"] for x in v) / 1024.0,
            units=v[0]["n_units"],
            exits=",".join(sorted({x["exit"] for x in v})),
        )

    for phase, blas in sorted({(k[0], k[1]) for k in agg}):
        title = ("OPENBLAS_NUM_THREADS=1" if blas == "1"
                 else "BLAS threading at machine default (%s)" % blas)
        print("\n### Phase %s -- %s\n" % (phase, title))
        print("| stage | engine | nThreads | reps | wall median (s) | wall min-max (s) "
              "| user (s) | sys (s) | user/wall | max RSS (MB) | rows out |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        for k in sorted(agg, key=lambda k: (k[2], k[3], k[4])):
            if (k[0], k[1]) != (phase, blas):
                continue
            a = agg[k]
            uw = a["user"] / a["wall"] if a["wall"] else 0
            print("| %s | %s | %d | %d | %s | %s-%s | %s | %s | %.2f | %.0f | %s |"
                  % (k[2], k[3], k[4], a["n"], fmt(a["wall"]), fmt(a["wall_lo"]),
                     fmt(a["wall_hi"]), fmt(a["user"]), fmt(a["sys"]), uw,
                     a["rss"], a["units"]))

        # ---- C++ vs R ----
        print("\n**C++ vs R (same phase, R has no thread knob in step 2)**\n")
        print("| stage | nThreads (C++) | C++ wall (s) | R wall (s) | speedup (R/C++) |")
        print("|---|---|---|---|---|")
        for stage in ("step1", "region", "single"):
            rk = [k for k in agg if k[:3] == (phase, blas, stage) and k[3] == "r"]
            if not rk:
                continue
            for ck in sorted([k for k in agg if k[:4] == (phase, blas, stage, "cpp")],
                             key=lambda k: k[4]):
                # for step 1, compare like-for-like on nThreads; otherwise use
                # the single R number
                match = [k for k in rk if k[4] == ck[4]] or rk
                rw = agg[match[0]]["wall"]
                cw = agg[ck]["wall"]
                print("| %s | %d | %s | %s | %.2fx |"
                      % (stage, ck[4], fmt(cw), fmt(rw), rw / cw if cw else 0))

        # ---- C++ thread scaling ----
        print("\n**C++ thread scaling (relative to its own nThreads=1)**\n")
        print("| stage | nThreads | wall (s) | speedup vs T=1 | parallel efficiency |")
        print("|---|---|---|---|---|")
        for stage in ("step1", "region", "single"):
            base = [k for k in agg if k[:4] == (phase, blas, stage, "cpp") and k[4] == 1]
            if not base:
                continue
            b = agg[base[0]]["wall"]
            for ck in sorted([k for k in agg if k[:4] == (phase, blas, stage, "cpp")],
                             key=lambda k: k[4]):
                w = agg[ck]["wall"]
                sp = b / w if w else 0
                print("| %s | %d | %s | %.2fx | %.0f%% |"
                      % (stage, ck[4], fmt(w), sp, 100.0 * sp / ck[4]))

    if warm:
        print("\n### Warm-up runs (EXCLUDED from every number above)\n")
        print("| phase | blas | stage | engine | nThreads | wall (s) |")
        print("|---|---|---|---|---|---|")
        for r in warm:
            print("| %s | %s | %s | %s | %s | %s |"
                  % (r["phase"], r["blas"], r["stage"], r["engine"], r["threads"],
                     fmt(r["wall_s"])))

    bad = [r for r in rows if r["exit"] != "0"]
    if bad:
        print("\n**NON-ZERO EXITS:** %d run(s)" % len(bad))
        for r in bad:
            print("  - %s %s %s T%s rep%s exit=%s"
                  % (r["phase"], r["stage"], r["engine"], r["threads"], r["rep"], r["exit"]))


if __name__ == "__main__":
    main()
