#!/usr/bin/env python3
"""make_summary.py OUT_ROOT -- one row per run (summary.csv) and one row per comparison (compare.csv)
from $OUT_ROOT/cells/*/cell.json and $OUT_ROOT/compare/*.json. Column meanings: docs/collaborator_tests.md."""
import csv, glob, json, os, sys

R = sys.argv[1]
RUN_COLS = ["block", "cell", "path", "trait_type", "P", "grm", "fast_test", "firth", "geno", "binary", "extra",
            "nthreads", "rc", "wall_s", "user_s", "sys_s", "cpu_pct", "peak_rss_kb",
            "r_processes", "r_peak_rss_sum_kb", "r_proc_wall_max_s",
            "gpu_peak_mem_mib", "gpu_mem_base_mib", "gpu_buffers_mib_log", "gpu_name", "gpu_coverage", "gpu_refused",
            "startup_s", "main_loop_s", "n_markers_tested_max", "n_result_files",
            "n_pairs", "n_need_spa", "n_spa", "n_lowmac", "n_er", "n_need_firth", "n_firth", "n_firth_conv",
            "n_need_fast", "n_firth_log", "n_spa_is_spa_col", "counts_source",
            "cache_method", "cache_cold", "cache_resident_after_MB", "cache_input_MB", "fstypes", "fs_inputs_MB",
            "log_errors", "cpu_vs_gpu_identical", "compare_label", "md5_all", "commit", "start"]
CMP_COLS = ["label", "kind", "identical", "n_paired", "n_md5_identical", "files_only_a", "files_only_b",
            "rows_a", "rows_b", "only_a", "only_b", "rows_text_differ",
            "p_n_differ", "p_max_abs", "p_max_rel", "beta_n_differ", "beta_max_abs", "beta_max_rel",
            "se_n_differ", "se_max_abs", "se_max_rel", "na_mismatch", "max_abs_dlog10p",
            "n_p_lt_5e-8_a", "n_p_lt_5e-8_b", "cross_5e-8_a_only", "cross_5e-8_b_only",
            "n_p_lt_1e-5_a", "n_p_lt_1e-5_b", "cross_1e-5_a_only", "cross_1e-5_b_only",
            "is_spa_true_a", "is_spa_true_b", "is_spa_disagree"]

cmps = {}
rows = []
for f in sorted(glob.glob(os.path.join(R, "compare", "*.json"))):
    c = json.load(open(f))
    lab = c.get("label") or os.path.basename(f)[:-5]
    cmps[lab] = c
    t = c.get("total", {})
    row = {"label": lab, "kind": lab.split("__")[0], "identical": c.get("identical"),
           "n_paired": c.get("n_paired"), "n_md5_identical": c.get("n_md5_identical"),
           "files_only_a": c.get("files_only_a"), "files_only_b": c.get("files_only_b")}
    for k in ("rows_a", "rows_b", "only_a", "only_b", "rows_text_differ", "max_abs_dlog10p",
              "n_p_lt_5e-8_a", "n_p_lt_5e-8_b", "cross_5e-8_a_only", "cross_5e-8_b_only",
              "n_p_lt_1e-5_a", "n_p_lt_1e-5_b", "cross_1e-5_a_only", "cross_1e-5_b_only",
              "is_spa_true_a", "is_spa_true_b", "is_spa_disagree"):
        row[k] = t.get(k)
    for c_, p_ in (("p.value", "p"), ("BETA", "beta"), ("SE", "se")):
        st = t.get(c_, {})
        row[p_ + "_n_differ"] = st.get("n_differ"); row[p_ + "_max_abs"] = st.get("max_abs"); row[p_ + "_max_rel"] = st.get("max_rel")
    row["na_mismatch"] = sum(t.get(c_, {}).get("na_mismatch", 0) for c_ in ("p.value", "BETA", "SE")) if t else None
    rows.append(row)
with open(os.path.join(R, "compare.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, CMP_COLS); w.writeheader(); w.writerows(rows)

runs = []
for f in sorted(glob.glob(os.path.join(R, "cells", "*", "cell.json"))):
    c = json.load(open(f))
    n = c["cell"]
    lab = ""
    if n.startswith(("main_", "rare_")):
        lab = "cpu_vs_gpu__" + n.replace("_cpu_", "_").replace("_gpu_", "_")
    elif n.startswith("vsr_"):
        lab = "cpp_vs_R__" + n.split("_", 2)[2]
    elif n.startswith("pgen_"):
        lab = "pgen_vs_bed__" + n.split("_", 2)[2]
    c["compare_label"] = lab if lab in cmps else ""
    if n.startswith(("main_", "rare_")) and lab in cmps:
        c["cpu_vs_gpu_identical"] = cmps[lab].get("identical")
    runs.append(c)
order = {"main": 0, "rare": 1, "stage": 2, "pgen": 3, "vsr": 4}
runs.sort(key=lambda c: (order.get(c["block"], 9), c["P"], c["grm"], str(c["fast_test"]), c["firth"], c["path"]))
with open(os.path.join(R, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, RUN_COLS, extrasaction="ignore"); w.writeheader(); w.writerows(runs)
print("%d runs, %d comparisons" % (len(runs), len(rows)))
