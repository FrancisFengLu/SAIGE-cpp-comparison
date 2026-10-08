#!/usr/bin/env python3
"""cell_post.py CELL_DIR  -- summarise one finished run into CELL_DIR/cell.json.

Reads CELL_DIR/meta.json (written by run_matrix.sh), log.txt (step 2 stdout+stderr with
/usr/bin/time -v appended, or for R cells the per-trait logs), gpu_mem.txt (nvidia-smi samples,
MiB), cache.json (cache_evict.py), routes/*.route (SAIGE_STEP2_ROUTE_DUMP) and out/*.txt.
Writes md5.txt (one line per result file) and cell.json. Nothing per-variant is written.
Route file: 9 bytes per output row {u8 route, f64 pre-SPA p}; bit0 needSPA, bit1 needFirth,
bit2 needFast, bit3 MAC <= MACCutoffforER, bit4 SPA/ER converged, bit5 Firth fitted,
bit6 Firth converged; 0xFF = no output row.
"""
import glob, hashlib, json, os, re, sys

D = sys.argv[1]
meta = json.load(open(os.path.join(D, "meta.json")))
r = dict(meta)


def md5(fn):
    h = hashlib.md5()
    with open(fn, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tv(t, key, conv=float):
    m = re.search(re.escape(key) + r":\s*(\S+)", t)
    if not m:
        return None
    v = m.group(1)
    if key.startswith("Elapsed"):
        parts = [float(x) for x in v.split(":")]
        return round(sum(p * 60 ** i for i, p in enumerate(reversed(parts))), 2)
    try:
        return conv(v.rstrip("%"))
    except ValueError:
        return None


log = open(os.path.join(D, "log.txt"), errors="replace").read() if os.path.exists(os.path.join(D, "log.txt")) else ""
r["rc"] = int(open(os.path.join(D, "rc")).read().strip()) if os.path.exists(os.path.join(D, "rc")) else None
r["wall_s"] = tv(log, "Elapsed (wall clock) time (h:mm:ss or m:ss)")
r["user_s"] = tv(log, "User time (seconds)")
r["sys_s"] = tv(log, "System time (seconds)")
r["cpu_pct"] = tv(log, "Percent of CPU this job got", int)
r["peak_rss_kb"] = tv(log, "Maximum resident set size (kbytes)", int)
r["fs_inputs_MB"] = (lambda v: round(v * 512 / 1e6, 1) if v is not None else None)(tv(log, "File system inputs", int))

# R cells: one /usr/bin/time -v block per trait process in rlogs/<trait>.time
rt = sorted(glob.glob(os.path.join(D, "rlogs", "*.time")))
if rt:
    rss = [tv(open(f).read(), "Maximum resident set size (kbytes)", int) or 0 for f in rt]
    r["r_processes"] = len(rt)
    r["r_peak_rss_sum_kb"] = sum(rss)
    r["r_peak_rss_max_kb"] = max(rss)
    r["r_proc_wall_max_s"] = max(tv(open(f).read(), "Elapsed (wall clock) time (h:mm:ss or m:ss)") or 0 for f in rt)

# GPU memory samples
gm = os.path.join(D, "gpu_mem.txt")
if os.path.exists(gm):
    vals = [int(x) for x in re.findall(r"^\s*(\d+)\s*$", open(gm).read(), re.M)]
    base = None
    if os.path.exists(os.path.join(D, "gpu_mem_base.txt")):
        b = re.findall(r"\d+", open(os.path.join(D, "gpu_mem_base.txt")).read())
        base = int(b[0]) if b else None
    r["gpu_mem_samples"] = len(vals)
    r["gpu_mem_base_mib"] = base
    r["gpu_peak_mem_mib"] = (max(vals) - (base or 0)) if vals else None

# cache
cj = os.path.join(D, "cache.json")
if os.path.exists(cj):
    c = json.load(open(cj))
    r["cache_method"] = c["method"]
    r["cache_cold"] = c["cold"]
    r["cache_resident_after_MB"] = round(c["resident_after"] / 1e6, 2)
    r["cache_input_MB"] = round(c["bytes_total"] / 1e6, 1)
    r["fstypes"] = ",".join(c["fstypes"])

# step-2 log lines (C++ cells)
if r.get("path") in ("cpu", "gpu"):
    m = re.search(r"useGPU: refused, running on the CPU \((.*)\)", log)
    r["gpu_refused"] = m.group(1)[:200] if m else ""
    m = re.search(r"GPU coverage: (\d+) / (\d+) pairs", log)
    r["gpu_coverage"] = "%s/%s" % m.groups() if m else ""
    m = re.search(r"useGPU: (.*?), (\d+) MiB, (sm_\d+)", log)
    r["gpu_name"] = m.group(1).strip() if m else ""
    m = re.search(r"(\d+) MiB on the device", log)
    r["gpu_buffers_mib_log"] = int(m.group(1)) if m else None
    tested = [int(x) for x in re.findall(r"(?:\] |^)(\d+) markers were tested", log, re.M)]
    r["n_traits_logged"] = len(tested)
    r["n_markers_tested_max"] = max(tested) if tested else None
    m = re.search(r"\[TIMING\] 50_before_main_loop\s+\S+\s+total=([0-9.eE+-]+)s", log)
    r["startup_s"] = round(float(m.group(1)), 3) if m else None
    m = re.search(r"\[TIMING\] 60_after_main_loop\s+\S+\s+total=([0-9.eE+-]+)s", log)
    r["main_loop_s"] = round(float(m.group(1)) - (r["startup_s"] or 0), 3) if m else None
    r["n_firth_log"] = sum(int(x) for x in re.findall(r"Firth approx was applied to (\d+) markers", log))
    r["log_errors"] = len(re.findall(r"(?i)\berror\b|terminate called|Segmentation fault|std::bad_alloc", log))

# route dump
rf = sorted(glob.glob(os.path.join(D, "routes", "*.route")))
if rf:
    c = dict(n_pairs=0, n_need_spa=0, n_spa=0, n_lowmac=0, n_er=0, n_need_firth=0, n_firth=0, n_firth_conv=0, n_need_fast=0)
    from collections import Counter
    hist = Counter()
    for f in rf:
        hist.update(open(f, "rb").read()[0::9])     # route byte of every row
    for x, n in hist.items():
        if x == 0xFF:
            continue
        low, conv = (x >> 3) & 1, (x >> 4) & 1
        c["n_pairs"] += n
        c["n_need_spa"] += n * (x & 1)
        c["n_need_firth"] += n * ((x >> 1) & 1)
        c["n_need_fast"] += n * ((x >> 2) & 1)
        c["n_lowmac"] += n * low
        c["n_er"] += n * (low & conv)
        c["n_spa"] += n * (conv & (1 - low))
        c["n_firth"] += n * ((x >> 5) & 1)
        c["n_firth_conv"] += n * ((x >> 6) & 1)
    r.update(c)

# result files: md5 only
outs = sorted(glob.glob(os.path.join(D, "out", "*.txt")))
with open(os.path.join(D, "md5.txt"), "w") as f:
    for o in outs:
        f.write("%s  %s\n" % (md5(o), os.path.basename(o)))
r["n_result_files"] = len(outs)
r["md5_all"] = hashlib.md5(open(os.path.join(D, "md5.txt"), "rb").read()).hexdigest() if outs else ""
if outs and not rf:
    # no route dump (R cells, and the CPU single-trait path that P = 1 takes): count from the
    # result files -- rows and Is.SPA = true; ER and Firth counts are not available this way
    npairs = nspa = 0
    for o in outs:
        with open(o) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            k = hdr.index("Is.SPA") if "Is.SPA" in hdr else None
            for line in fh:
                npairs += 1
                if k is not None and line.split("\t")[k].lower() == "true":
                    nspa += 1
    r["n_pairs"], r["n_spa_is_spa_col"] = npairs, nspa
    r["counts_source"] = "result files"
elif rf:
    r["counts_source"] = "route dump"
json.dump(r, open(os.path.join(D, "cell.json"), "w"), indent=1, sort_keys=True)
print("%-44s rc=%s wall=%ss rss=%sMB gpu=%sMiB pairs=%s spa=%s er=%s firth=%s cache=%s" % (
    r["cell"], r["rc"], r["wall_s"], (r["peak_rss_kb"] or 0) // 1024, r.get("gpu_peak_mem_mib"),
    r.get("n_pairs"), r.get("n_spa"), r.get("n_er"), r.get("n_firth"), r.get("cache_method")))
