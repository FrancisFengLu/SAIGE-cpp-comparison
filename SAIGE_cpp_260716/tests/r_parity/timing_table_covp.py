#!/usr/bin/env python3
"""timing_table_covp.py OUT: the timing_covp.sh runs under OUT/<label>/p<p>/r<rep>/log.txt as a
Markdown table -- per label (base / lib64 / own64 / lib32) and p: SPA pairs, SPA kernel seconds and
us per pair, Firth pairs, Firth kernel seconds and us per pair (min over the repetitions, and the
spread), plus the wall time of the run."""
import sys, os, re, glob
OUT = sys.argv[1]
def parse(L):
    d = {}
    for l in open(L):
        m = re.search(r"device SPA: (\d+) pairs solved in ([\d.e+-]+) s of kernel time \(([\d.e+-]+) us/pair\)", l)
        if m: d["spaN"], d["spaS"], d["spaUs"] = int(m.group(1)), float(m.group(2)), float(m.group(3))
        m = re.search(r"device Firth: (\d+) pairs fitted in ([\d.e+-]+) s of kernel time \(([\d.e+-]+) us/pair\)", l)
        if m: d["fiN"], d["fiS"], d["fiUs"] = int(m.group(1)), float(m.group(2)), float(m.group(3))
        m = re.search(r"gpuSpa: device setup failed \((.*)\)", l)
        if m: d["spaOff"] = m.group(1)
        m = re.search(r"Total time: ([\d.]+)|total ([\d.]+) s|wall ([\d.]+)", l)
    return d
print("| run | p | SPA pairs | SPA kernel s (min / max) | SPA us per pair (min) | Firth pairs | Firth kernel s (min / max) | Firth us per pair (min) |")
print("|---|--:|--:|--:|--:|--:|--:|--:|")
for lab in ["base", "lib64", "own64", "lib32"]:
    for pd in sorted(glob.glob(f"{OUT}/{lab}/p*"), key=lambda x: int(x.rsplit("p", 1)[1])):
        p = pd.rsplit("p", 1)[1]
        reps = [parse(L) for L in sorted(glob.glob(f"{pd}/r*/log.txt"))]
        reps = [r for r in reps if r]
        if not reps: continue
        if "spaOff" in reps[0]:
            print(f"| {lab} | {p} | device SPA off: {reps[0]['spaOff'][:80]} | | | | | |"); continue
        sN = reps[0].get("spaN", 0); sS = [r["spaS"] for r in reps if "spaS" in r]; sU = [r["spaUs"] for r in reps if "spaUs" in r]
        fN = reps[0].get("fiN", 0); fS = [r["fiS"] for r in reps if "fiS" in r]; fU = [r["fiUs"] for r in reps if "fiUs" in r]
        f = lambda v: f"{min(v):.4f} / {max(v):.4f}" if v else "-"
        g = lambda v: f"{min(v):.1f}" if v else "-"
        print(f"| {lab} | {p} | {sN} | {f(sS)} | {g(sU)} | {fN} | {f(fS)} | {g(fU)} |")
