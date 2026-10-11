#!/usr/bin/env python3
"""sparse_table.py: every case under /opt/saige/logs/gpuprep/sparse, the GPU run against the CPU run of the
same binary (block inverse), worst over the traits; route bytes compared exactly; the device log lines."""
import os, sys, json, subprocess, struct, math, re
S = os.environ.get("SPARSEB", "/opt/saige/logs/gpuprep/sparse")
GD = os.environ.get("GDIR", "gpu")
CMP = os.environ.get("CMP_R", "/opt/saige/logs/rdefaults/cmp_r.py")
def cmp(a, b, label):
    if not (os.path.exists(a) and os.path.exists(b)): return None
    return json.loads(subprocess.check_output([sys.executable, CMP, a, b, "--label", label, "--json"]).decode())
def fmt(r):
    if r is None: return "| (missing) |||||||||"
    m = r["maxrel"]
    return (f"| {r['label']} | {r['rows']} | {m['p.value']:.1e} | {m['BETA']:.1e} | {m['SE']:.1e} | {m['Tstat']:.1e} | {m['var']:.1e} "
            f"| {r['beyondPrint']['p.value']}/{r['beyondPrint']['SE']}/{r['beyondPrint']['var']} | {r['cross5e8']}/{r['cross1e5']} "
            f"| {r['isSPAmismatch']} ({r['nSPA_A']}) | {r['firthMismatch']} ({r['firthA']}/{r['firthB']}) |")
hdr = ("| case | rows | p rel | BETA rel | SE rel | Tstat rel | var rel | beyond print p/SE/var | cross 5e-8/1e-5 | Is.SPA mism (SPA rows) | Firth mism (A/B rows) |\n"
       "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
def worst(rs):
    rs = [r for r in rs if r]
    if not rs: return None
    w = dict(rs[0]); w["rows"] = sum(r["rows"] for r in rs)
    for k in ["p.value", "BETA", "SE", "Tstat", "var"]:
        w["maxrel"][k] = max(r["maxrel"][k] for r in rs); w["beyondPrint"][k] = sum(r["beyondPrint"][k] for r in rs)
    for k in ["cross5e8", "cross1e5", "isSPAmismatch", "nSPA_A", "firthMismatch", "firthA", "firthB"]:
        w[k] = sum(r[k] for r in rs)
    return w
def tested(logpath):
    d = {}
    if not os.path.exists(logpath): return d
    for l in open(logpath):
        m = re.match(r"\[(\w+)\] (\d+) markers were tested", l.strip())
        if m: d[m.group(1)] = int(m.group(2))
    return d
cases = sorted(d for d in os.listdir(S) if os.path.isdir(os.path.join(S, d)) and os.path.exists(os.path.join(S, d, GD, "log.txt")))
print("## sparse GRM: GPU path vs CPU path (same binary, block inverse), worst over the traits\n")
print(hdr)
notes = []
for c in cases:
    G, C = f"{S}/{c}/{GD}", f"{S}/{c}/cpu"
    traits = sorted(f[:-4] for f in os.listdir(G) if f.endswith(".txt") and f != "log.txt")
    rs = [cmp(f"{G}/{t}.txt", f"{C}/{t}.txt", f"{c} {t}") for t in traits]
    w = worst(rs)
    if w: w["label"] = c
    print(fmt(w))
    tg, tc = tested(f"{G}/log.txt"), tested(f"{C}/log.txt")
    bad = [t for t in traits if tg.get(t) != tc.get(t)]
    nd = nb = 0
    for t in traits:
        a, b = f"{G}/routes/{t}.route", f"{C}/routes/{t}.route"
        if os.path.exists(a) and os.path.exists(b):
            da, db = open(a, "rb").read(), open(b, "rb").read()
            n = min(len(da), len(db)) // 9; nb += n; nd += abs(len(da) - len(db)) // 9
            for i in range(n):
                if da[9 * i] != db[9 * i]: nd += 1
    lines = [l.strip() for l in open(f"{G}/log.txt") if l.strip().startswith(("gpuSparse:", "gpuDeviceStats: ", "GPU coverage", "gate:", "useGPU: refused"))]
    notes.append((c, bad, nb, nd, [l[:220] for l in lines if "pairs took" in l or "the sparse statistic" in l or "cross-term" in l or "refused" in l or l.startswith("GPU coverage")]))
print("\n## markers tested per trait, route bytes, device lines\n")
for c, bad, nb, nd, lines in notes:
    print(f"- {c}: tested-count mismatches {len(bad)}{' ' + str(bad) if bad else ''}; route bytes compared {nb}, differing {nd}")
    for l in lines: print(f"  - {l}")
