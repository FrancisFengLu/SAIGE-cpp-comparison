#!/usr/bin/env python3
"""stage_table.py RUNDIR [RUNDIR ...] [--labels a,b,..] [--csv OUT.csv]

One table, stage | pairs | wall s | thread-s, one column group per run, for step-2 runs of
SAIGE-cpp-comparison (multi-trait loop, CPU path or GPU path). RUNDIR holds log.txt
(stdout+stderr incl. [TIMING] marks and /usr/bin/time -v) and optionally routes/*.route
(SAIGE_STEP2_ROUTE_DUMP). A path to a log file is accepted too.

What each source gives:
  * [TIMING] marks (stderr, every build): startup and main-loop wall.
  * [phase] lines (PHASE_TIMING=1 builds only): thread-s ("cpu-s" in the log; it is the sum of
    per-thread omp_get_wtime intervals, NOT getrusage CPU time) and call counts per slot.
    CPU path: every slot. GPU path: only the slots inside the scalar getMarkerPval
    (score recompute / SPA / ER / Firth / 4 vec allocs) -- the CPU scalar calls the GPU
    path still makes; the read/GEMM/finalize slots stay 0 there.
  * [mt breakdown] (CPU path): output write wall.
  * [gpu breakdown] / [gpu device time] / [gpu pipeline] / [gpu busy] / [gpu overlap],
    "device SPA:", "device ER:", "device Firth:", "gpuSparse:", "gate:" (GPU path, every build).
  * per-trait "N markers were tested (...)" lines, "Firth approx was applied" lines.
  * routes/*.route: { u8 route, f64 gateP } per row; bit0 needSPA, bit1 needFirth, bit2 needFast,
    bit3 ER/low MAC, bit4 isSPAConverge, bit5 is_Firth, bit6 is_FirthConverge; 0xFF = no row.

wall s for CPU-path stages inside the block-parallel region is thread-s / nThreads (the
loop's own convention): it is the stage's share of the region wall, not a measured interval.
With gpuOverlap on, GPU-path stage walls are per-thread walls of concurrent workers and do
NOT add up to the main loop; use a gpuOverlap: false run for attribution.
"""
import argparse, glob, os, re, struct, sys

NUM = r"([-+0-9.eE]+|nan|inf)"


def fnum(s):
    try:
        return float(s)
    except Exception:
        return float("nan")


def read_log(path):
    txt = open(path, errors="replace").read()
    return txt


def parse(path):
    log = path if os.path.isfile(path) else os.path.join(path, "log.txt")
    rundir = os.path.dirname(os.path.abspath(log))
    t = read_log(log)
    r = {"log": log}
    # ---- [TIMING] marks ----
    marks = {m.group(1): fnum(m.group(2)) for m in re.finditer(r"\[TIMING\] (\S+)\s+\+\S+\s+total=" + NUM + "s", t)}
    r["marks"] = marks
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", t)
    if m:
        parts = [float(x) for x in m.group(1).split(":")]
        r["elapsed"] = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
    m = re.search(r"User time \(seconds\): " + NUM, t); r["user"] = fnum(m.group(1)) if m else None
    m = re.search(r"System time \(seconds\): " + NUM, t); r["sys"] = fnum(m.group(1)) if m else None
    m = re.search(r"Percent of CPU this job got: (\d+)%", t); r["pcpu"] = int(m.group(1)) if m else None
    m = re.search(r"File system inputs: (\d+)", t); r["fsin_MB"] = int(m.group(1)) * 512 / 1e6 if m else None
    m = re.search(r"nThreads:\s+(\d+)", t); r["nthreads"] = int(m.group(1)) if m else None
    # ---- [phase] (PHASE_TIMING builds) ----
    ph = {}
    for m in re.finditer(r"\[phase\]\s+(\S.*?\S)\s+cpu-s\s+" + NUM + r"\s+wall-s\s+" + NUM + r"\s+calls (\d+)", t):
        ph[m.group(1).strip()] = (fnum(m.group(2)), int(m.group(4)))
    m = re.search(r"\[phase\] threads (\d+)", t)
    r["pt_threads"] = int(m.group(1)) if m else None
    cnt = {m.group(1).strip(): int(m.group(2)) for m in re.finditer(r"\[phase\] count (.+?)\s+(\d+)\s*$", t, re.M)}
    r["phase"], r["pcount"] = ph, cnt
    # ---- CPU path ----
    m = re.search(r"\[mt breakdown\] output write " + NUM + " s", t)
    r["mt_write"] = fnum(m.group(1)) if m else None
    m = re.search(r"MT batch coverage: (\d+) / (\d+) pairs", t)
    r["mt_cov"] = (int(m.group(1)), int(m.group(2))) if m else None
    # ---- GPU path ----
    r["gpu"] = "[gpu breakdown]" in t
    m = re.search(r"\[gpu breakdown\] read\+QC\+stage " + NUM + r" s, device call " + NUM +
                  r" s, tail\+finalize " + NUM + r" s, device SPA \+ post " + NUM +
                  r" s, device Firth \+ post " + NUM + r"(?: s, device ER \+ post " + NUM + r")? s, output write " + NUM, t)
    if m:
        g = [fnum(x) if x is not None else 0.0 for x in m.groups()]
        r["gb"] = dict(zip(["read", "device", "tail", "spa", "firth", "er", "write"], g))
    m = re.search(r"\[gpu device time\] H2D " + NUM + r" s, decode " + NUM + r" s, GEMM " + NUM +
                  r" s, popcount " + NUM + r" s, D2H " + NUM, t)
    if m:
        r["gdev"] = dict(zip(["h2d", "decode", "gemm", "popc", "d2h"], [fnum(x) for x in m.groups()]))
    m = re.search(r"\[gpu pipeline\] main loop " + NUM + r" s; gpuPrefetch (\w+)(?:.*?reader wall " + NUM +
                  r" s on \d+ thread\(s\), consumer waited " + NUM + ")?", t)
    if m:
        r["gpipe"] = {"loop": fnum(m.group(1)), "prefetch": m.group(2),
                      "reader": fnum(m.group(3)) if m.group(3) else None,
                      "readwait": fnum(m.group(4)) if m.group(4) else None}
    m = re.search(r"\[gpu busy\] device calls (\d+), union " + NUM + r" s, sum " + NUM, t)
    if m:
        r["gbusy"] = fnum(m.group(2))
    r["overlap"] = "[gpu overlap]" in t
    m = re.search(r"\[gpu overlap\].*?waited " + NUM + r" s for the scan, post " + NUM +
                  r" s \(waited " + NUM + r" s for the SPA worker", t)
    if m:
        r["ov"] = {"wait_scan": fnum(m.group(1)), "post": fnum(m.group(2)), "wait_spa": fnum(m.group(3))}
    m = re.search(r"gate: needSPA (\d+), needFirth (\d+), needFast (\d+), ER \(MAC <= \d+\) (\d+) pairs", t)
    if m:
        r["ggate"] = dict(zip(["needSPA", "needFirth", "needFast", "lowMAC"], map(int, m.groups())))
    m = re.search(r"device SPA: (\d+) pairs solved in " + NUM + r" s of kernel time.*?; (\d+) kept the device result, (\d+) asked for Firth", t)
    if m:
        r["dspa"] = (int(m.group(1)), fnum(m.group(2)), int(m.group(3)), int(m.group(4)))
    m = re.search(r"device ER: (\d+) pairs in " + NUM + r" s of kernel time.*?; (\d+) finished on the device, (\d+) asked", t)
    if m:
        r["der"] = (int(m.group(1)), fnum(m.group(2)), int(m.group(3)), int(m.group(4)))
    m = re.search(r"device Firth: (\d+) pairs fitted in " + NUM + r" s of kernel time.*?; (\d+) of the (\d+) Firth fits", t)
    if m:
        r["dfirth"] = (int(m.group(1)), fnum(m.group(2)), int(m.group(3)), int(m.group(4)))
    m = re.search(r"gpuSparse: cross-term kernel " + NUM + r" s over \d+ slots; (\d+) trait\(s\) with the sparse first pass.*?; (\d+) fast-test recompute pairs finished on the device; (\d+) CPU scalar calls", t)
    if m:
        r["dsparse"] = (fnum(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))
    m = re.search(r"GPU coverage: (\d+) / (\d+) pairs", t)
    r["gpu_cov"] = (int(m.group(1)), int(m.group(2))) if m else None
    # ---- sparse-Sigma solves (s2_block_solve.cpp reportAll, every build, both paths) ----
    r["sparse_solve"] = {}
    for m in re.finditer(r"^\s+(PCG|block inverse):\s+(\d+) solves, " + NUM + " s total", t, re.M):
        r["sparse_solve"][m.group(1)] = (int(m.group(2)), fnum(m.group(3)))
    # ---- per-trait lines (both paths) ----
    nb = nf = nt = 0
    for m in re.finditer(r"\] (\d+) markers were tested \((.*?)\)\.", t):
        nt += int(m.group(1))
        s = m.group(2)
        a = re.search(r"(\d+) (?:batched|on the GPU)", s); nb += int(a.group(1)) if a else 0
        a = re.search(r"(\d+) via the scalar", s); nf += int(a.group(1)) if a else 0
    r["tested"], r["scalar"] = nt, nf
    r["firth_applied"] = sum(int(m.group(1)) for m in re.finditer(r"Firth approx was applied to (\d+) markers", t))
    # ---- route dump ----
    r["routes"] = None
    rf = sorted(glob.glob(os.path.join(rundir, "routes", "*.route")))
    if rf:
        c = dict(rows=0, needSPA=0, needFirth=0, needFast=0, lowMAC=0, ER_exec=0, SPA_conv=0, Firth=0, FirthConv=0, gated=0)
        for f in rf:
            b = open(f, "rb").read()
            for k in range(0, len(b), 9):
                x = b[k]
                if x == 0xFF:
                    continue
                c["rows"] += 1
                c["needSPA"] += x & 1; c["needFirth"] += (x >> 1) & 1; c["needFast"] += (x >> 2) & 1
                low = (x >> 3) & 1; conv = (x >> 4) & 1
                c["lowMAC"] += low
                c["ER_exec"] += low & conv
                c["SPA_conv"] += conv & (1 - low)
                c["Firth"] += (x >> 5) & 1; c["FirthConv"] += (x >> 6) & 1
                gp = struct.unpack_from("<d", b, k + 1)[0]
                c["gated"] += (gp == gp)
        r["routes"] = c
    return r


def rows_for(r):
    """List of (stage, pairs, wall_s, thread_s, note)."""
    out = []
    mk = r["marks"]
    nt = r.get("pt_threads") or r.get("nthreads") or 1
    ph, pc, rt = r["phase"], r["pcount"], r["routes"] or {}

    def P(name):
        v = ph.get(name)
        return v if v else (None, None)

    if "50_before_main_loop" in mk:
        out.append(("startup (main -> loop: yaml, null models, reader)", None, mk["50_before_main_loop"], None, ""))
    if r["gpu"]:
        gb = r.get("gb", {}); gd = r.get("gdev", {}); gp = r.get("gpipe", {})
        ov = r["overlap"]
        out.append(("read+QC+stage (host)", r["tested"] or None, gb.get("read"), None,
                    "prefetch reader wall, overlapped; consumer waited %.3f s" % gp["readwait"]
                    if gp.get("prefetch") == "on" and gp.get("readwait") is not None else ""))
        dsum = sum(gd.values()) if gd else None
        out.append(("score scan: device call (H2D+decode+GEMM+popc+D2H)", rt.get("rows") or r["tested"] or None,
                    gb.get("device"), None, "event sum %.3f s (streams overlap)" % dsum if dsum is not None else ""))
        out.append(("host tail + gate + finalize (incl. CPU scalar calls)", rt.get("gated") or None, gb.get("tail"), None, ""))
        if r.get("dsparse"):
            ds = r["dsparse"]
            out.append(("sparse cross-term kernel / device fast-test recompute", ds[2], ds[0], None,
                        "%d trait(s) sparse first pass; %d CPU calls given device variance" % (ds[1], ds[3])))
        if r.get("dspa"):
            d = r["dspa"]
            out.append(("SPA: device + post", d[0], gb.get("spa"), None,
                        "kernel %.3f s; kept %d, Firth->CPU %d" % (d[1], d[2], d[3])))
        if r.get("der"):
            d = r["der"]
            out.append(("ER: device + post", d[0], gb.get("er"), None, "kernel %.3f s; finished on device %d" % (d[1], d[2])))
        if r.get("dfirth"):
            d = r["dfirth"]
            out.append(("Firth: device + post", d[0], gb.get("firth"), None,
                        "kernel %.3f s; %d of %d fits on device" % (d[1], d[2], d[3])))
        out.append(("write", None, gb.get("write"), None, ""))
        if ph:
            for nm, lab in [("fb: score recompute", "CPU scalar: score recompute"), ("fb: SPA", "CPU scalar: SPA"),
                            ("fb: ER", "CPU scalar: ER"), ("fb: Firth fit", "CPU scalar: Firth")]:
                ts, n = P(nm)
                if n:
                    out.append(("  " + lab, n, ts / nt, ts, "inside the rows above (wall = thread-s/threads)"))
        loop = gp.get("loop")
        if ov and r.get("ov"):
            o = r["ov"]
            out.append(("  main thread waited for scan / SPA worker", None, o["wait_scan"] + o["wait_spa"], None,
                        "gpuOverlap on: rows above run concurrently"))
        if r.get("gbusy") is not None:
            out.append(("  GPU busy (union of device calls)", None, r["gbusy"], None, ""))
        if not ov and loop is not None and gb:
            resid = loop - sum(gb.get(k, 0) for k in ["device", "tail", "spa", "firth", "er", "write"]) \
                - (gp.get("readwait") or 0 if gp.get("prefetch") == "on" else gb.get("read", 0))
            out.append(("  residual (glue between phases)", None, resid, None,
                        "main loop - device - tail - SPA - Firth - ER - write - read(wait)"))
        out.append(("main loop", None, loop, None, "gpuOverlap " + ("on" if ov else "off")))
        if loop is not None and "50_before_main_loop" in mk and "60_after_main_loop" in mk:
            out.append(("device setup (CUDA ctx, resident data, self-checks)", None,
                        mk["60_after_main_loop"] - mk["50_before_main_loop"] - loop, None,
                        "[TIMING] 50->60 minus [gpu pipeline] main loop"))
    else:
        if ph:
            rd, _ = P("read+QC+impute"); ge, _ = P("batch GEMM"); fi, _ = P("finalize(total)")
            ga, _ = P("gate"); fa, nfa = P("fallback(total)")
            out.append(("read+QC+impute", r["tested"] or None, rd / nt, rd, ""))
            out.append(("score scan: batch GEMM", pc.get("pairs reaching gate"), ge / nt, ge, ""))
            tail = fi - fa
            out.append(("host tail: gate + finalize (stat, AF, slots)", pc.get("pairs reaching gate"), tail / nt, tail,
                        "gate %.2f thread-s, AF gather %.2f" % (ga, P("fin: AF gather")[0] or 0)))
            out.append(("scalar fallback (total)", nfa, fa / nt, fa, ""))
            for nm, lab in [("fb: score recompute", "score recompute"), ("fb: SPA", "SPA"),
                            ("fb: ER", "ER"), ("fb: Firth fit", "Firth")]:
                ts, n = P(nm)
                out.append(("  " + lab, n, ts / nt, ts, ""))
            rest = fa - sum((P(x)[0] or 0) for x in ["fb: score recompute", "fb: SPA", "fb: ER", "fb: Firth fit"])
            out.append(("  rest (index build, allocs, plumbing)", None, rest / nt, rest, ""))
            pw, _ = P("par-region wall(master)")
            busy = rd + ge + fi
            out.append(("  barrier idle", None, pw - busy / nt, pw * nt - busy, "par-region wall %.3f s" % pw))
        out.append(("write", None, r.get("mt_write"), None, ""))
        if "50_before_main_loop" in mk and "60_after_main_loop" in mk:
            out.append(("main loop", None, mk["60_after_main_loop"] - mk["50_before_main_loop"], None, ""))
    for kind, (n, ts) in r.get("sparse_solve", {}).items():
        out.append(("  sparse Sigma^-1 solves (%s)" % kind, n, ts / nt, ts,
                    "thread-s inside scalar score recompute (wall = thread-s/threads)"))
    if r.get("elapsed") is not None:
        out.append(("TOTAL elapsed", None, r["elapsed"], (r["user"] or 0) + (r["sys"] or 0),
                    "thread-s here = user+sys CPU s; %s%% CPU; fs inputs %.0f MB" % (r["pcpu"], r["fsin_MB"] or 0)))
    return out


def counts_for(r):
    """Pair counts: gate decisions and stage triggers, from log lines and route dump."""
    c = {}
    rt = r["routes"] or {}
    c["pairs tested (output rows)"] = r["tested"] or rt.get("rows")
    c["pairs reaching batch gate"] = rt.get("gated") if rt else r["pcount"].get("pairs reaching gate")
    gg = r.get("ggate") or {}
    pc = r["pcount"] if not r["gpu"] else {}   # PT gate counters stay 0 on the GPU path
    if rt:                                     # the route dump is the same on both paths: prefer it
        gg = {"needSPA": rt["needSPA"], "needFirth": rt["needFirth"], "needFast": rt["needFast"], "lowMAC": rt["lowMAC"]}
    c["gate needSPA"] = gg.get("needSPA", pc.get("gate needSPA", rt.get("needSPA")))
    c["gate needFirth"] = gg.get("needFirth", pc.get("gate needFirth", rt.get("needFirth")))
    c["gate needFast (sparse/fast-test recompute)"] = gg.get("needFast", pc.get("gate needFast", rt.get("needFast")))
    c["MAC <= ER cutoff"] = gg.get("lowMAC", pc.get("binary pairs MAC<=ERcut", rt.get("lowMAC")))
    c["routed to scalar path"] = r["scalar"]
    if rt:
        c["route: SPA ran+converged (non-ER)"] = rt["SPA_conv"]
        c["route: ER executed"] = rt["ER_exec"]
        c["route: Firth fitted / converged"] = "%d / %d" % (rt["Firth"], rt["FirthConv"])
    c["Firth applied (log)"] = r["firth_applied"]
    return c


# canonical stage names for the side-by-side view (CPU row label / GPU row label -> one key)
CANON = {
    "startup (main -> loop: yaml, null models, reader)": "startup",
    "read+QC+impute": "read+QC", "read+QC+stage (host)": "read+QC",
    "score scan: batch GEMM": "score scan", "score scan: device call (H2D+decode+GEMM+popc+D2H)": "score scan",
    "host tail: gate + finalize (stat, AF, slots)": "host tail",
    "host tail + gate + finalize (incl. CPU scalar calls)": "host tail",
    "sparse cross-term kernel / device fast-test recompute": "sparse / fast-test recompute (device)",
    "scalar fallback (total)": "CPU scalar fallback (total)",
    "  score recompute": "  scalar: score recompute", "  CPU scalar: score recompute": "  scalar: score recompute",
    "  SPA": "SPA (CPU scalar)", "  CPU scalar: SPA": "SPA (CPU scalar)",
    "SPA: device + post": "SPA (device + post)",
    "  ER": "ER (CPU scalar)", "  CPU scalar: ER": "ER (CPU scalar)", "ER: device + post": "ER (device + post)",
    "  Firth": "Firth (CPU scalar)", "  CPU scalar: Firth": "Firth (CPU scalar)",
    "Firth: device + post": "Firth (device + post)",
    "  rest (index build, allocs, plumbing)": "  scalar: rest", "  barrier idle": "  barrier idle (CPU)",
    "  main thread waited for scan / SPA worker": "  overlap: main waited",
    "  GPU busy (union of device calls)": "  GPU busy (union)",
    "  residual (glue between phases)": "  residual (GPU glue)",
    "  sparse Sigma^-1 solves (PCG)": "  sparse Sigma^-1 solves (PCG)",
    "  sparse Sigma^-1 solves (block inverse)": "  sparse Sigma^-1 solves (block inverse)",
    "device setup (CUDA ctx, resident data, self-checks)": "device setup",
    "write": "write", "main loop": "main loop", "TOTAL elapsed": "TOTAL elapsed",
}
ORDER = ["startup", "read+QC", "score scan", "host tail", "sparse / fast-test recompute (device)",
         "CPU scalar fallback (total)", "  scalar: score recompute", "  sparse Sigma^-1 solves (PCG)",
         "  sparse Sigma^-1 solves (block inverse)", "SPA (device + post)", "SPA (CPU scalar)",
         "ER (device + post)", "ER (CPU scalar)", "Firth (device + post)", "Firth (CPU scalar)",
         "  scalar: rest", "  barrier idle (CPU)", "write", "  overlap: main waited", "  GPU busy (union)",
         "  residual (GPU glue)", "main loop", "device setup", "TOTAL elapsed"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--labels", default="")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    def deflab(x):
        x = os.path.normpath(x[:-len("/log.txt")] if x.endswith("/log.txt") else x)
        b = os.path.basename(x)
        return (os.path.basename(os.path.dirname(x)) + "/" + b) if re.match(r"rep\d+$", b) else b
    labels = a.labels.split(",") if a.labels else [deflab(x) for x in a.runs]
    R = [parse(x) for x in a.runs]

    def f(v, w=9, d=3):
        if v is None:
            return " " * (w - 1) + "-"
        if isinstance(v, str):
            return v.rjust(w)
        if isinstance(v, int):
            return str(v).rjust(w)
        return ("%.*f" % (d, v)).rjust(w)

    csvrows = []
    for lab, r in zip(labels, R):
        print("=== %s  (%s path%s, %s threads)  %s" % (lab, "GPU" if r["gpu"] else "CPU",
              (", gpuOverlap " + ("on" if r["overlap"] else "off")) if r["gpu"] else "",
              r.get("pt_threads") or r.get("nthreads"), r["log"]))
        print("%-56s %9s %9s %10s  %s" % ("stage", "pairs", "wall s", "thread-s", "note"))
        for st, n, w, ts, note in rows_for(r):
            print("%-56s %s %s %s  %s" % (st, f(n), f(w), f(ts, 10), note))
            csvrows.append((lab, "time", st, n, w, ts, note))
        print("%-56s %9s" % ("pair counts", ""))
        for k, v in counts_for(r).items():
            print("  %-54s %s" % (k, f(v, 13)))
            csvrows.append((lab, "count", k, v, None, None, ""))
        print()
    if len(R) > 1:
        # compact side by side: wall s per stage key that is comparable across paths
        def key(r):
            d = {}
            for st, n, w, ts, note in rows_for(r):
                d[CANON.get(st, st)] = (n, w)
            return d
        print("=== side by side (wall s / pairs)")
        keys = []
        for r in R:
            for st, *_ in rows_for(r):
                st = CANON.get(st, st)
                if st not in keys:
                    keys.append(st)
        keys.sort(key=lambda k: ORDER.index(k) if k in ORDER else len(ORDER))
        print("%-56s" % "stage" + "".join(" %22s" % l[:22] for l in labels))
        K = [key(r) for r in R]
        for st in keys:
            cells = []
            for d in K:
                if st in d:
                    n, w = d[st]
                    cells.append(" %22s" % ((f(w, 9) + " / " + (str(n) if n is not None else "-")).strip()))
                else:
                    cells.append(" %22s" % "")
            print("%-56s" % st + "".join(cells))
    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["run", "kind", "stage", "pairs", "wall_s", "thread_s", "note"])
            w.writerows(csvrows)


if __name__ == "__main__":
    main()
