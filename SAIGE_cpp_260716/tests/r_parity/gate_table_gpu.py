#!/usr/bin/env python3
"""gate_table_gpu.py OUT [datasets...]: the GPU runs under OUT/<ds>/<cfg>/gpu against R and against the
rdefaults single-trait path (worst over the dataset's traits), plus the per-trait marker counts."""
import sys, os, json, subprocess, re
OUT = sys.argv[1]
B = os.environ.get("RDEF", "/opt/saige/logs/rdefaults")
CMP = os.path.join(B, "cmp_r.py")
DS = sys.argv[2:] or ["audit", "bvs", "qt12"]
traits = {"audit": ["b1", "b2", "b3", "b4"], "bvs": ["b1", "b2", "b3", "b4"], "qt12": ["q1", "q2", "b1", "b2"]}
cfgs = {"audit": ["def", "defF", "adj", "adjF"], "bvs": ["def", "defF", "adj", "adjF"], "qt12": ["def", "adj", "defF", "adjF"]}
def cmp(a, b, label):
    if not (os.path.exists(a) and os.path.exists(b)): return None
    return json.loads(subprocess.check_output([sys.executable, CMP, a, b, "--label", label, "--json"]).decode())
def fmt(r):
    if r is None: return "| (missing) |||||||||"
    m = r["maxrel"]
    return (f"| {r['label']} | {r['rows']} | {m['p.value']:.1e} | {m['BETA']:.1e} | {m['SE']:.1e} | {m['Tstat']:.1e} | {m['var']:.1e} "
            f"| {r['beyondPrint']['p.value']}/{r['beyondPrint']['SE']}/{r['beyondPrint']['var']} | {r['cross5e8']}/{r['cross1e5']} "
            f"| {r['isSPAmismatch']} ({r['nSPA_A']}) | {r['firthMismatch']} ({r['firthA']}/{r['firthB']}) |")
hdr = ("| run | rows | p rel | BETA rel | SE rel | Tstat rel | var rel | beyond print p/SE/var | cross 5e-8/1e-5 | Is.SPA mism (SPA rows) | Firth mism (A/B rows) |\n"
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
def tested(logpath, trait=None):
    if not os.path.exists(logpath): return None
    for l in open(logpath):
        m = re.match(r"(\[(\w+)\] )?(\d+) markers were tested", l.strip())
        if m and (trait is None or m.group(2) == trait): return int(m.group(3))
    return None
print("## GPU path (new binary) vs R SAIGE 1.5.2 and vs the single-trait path\n")
print(hdr)
cnt = []
for ds in DS:
    for cfg in cfgs[ds]:
        rsR, rsS = [], []
        for t in traits[ds]:
            if cfg.endswith("F") and t.startswith("q"): continue
            rsR.append(cmp(f"{OUT}/{ds}/{cfg}/gpu/{t}.txt", f"{B}/R/{ds}/{cfg}/{t}.txt", f"{ds} {cfg} gpu {t}"))
            rsS.append(cmp(f"{OUT}/{ds}/{cfg}/gpu/{t}.txt", f"{B}/cpp/{ds}/{cfg}/single/{t}.txt", f"{ds} {cfg} gpu {t}"))
            ng = tested(f"{OUT}/{ds}/{cfg}/gpu/log.txt", t); ns = tested(f"{B}/cpp/{ds}/{cfg}/single/log_{t}.txt")
            cnt.append((ds, cfg, t, ng, ns))
        for rs, lab in [(rsR, "vs R"), (rsS, "vs single")]:
            w = worst(rs)
            if w: w["label"] = f"{ds} {cfg} gpu {lab}"
            print(fmt(w))
print("\n## markers tested per trait (GPU / single)\n")
bad = [c for c in cnt if c[3] != c[4]]
print(f"{len(cnt)} (dataset, cfg, trait) cells; differing counts: {len(bad)}" + (": " + str(bad) if bad else ""))
print("\n## device / host log lines\n")
for ds in DS:
    for cfg in cfgs[ds]:
        L = f"{OUT}/{ds}/{cfg}/gpu/log.txt"
        if not os.path.exists(L): continue
        for l in open(L):
            s = l.strip()
            if s.startswith("gpuDeviceStats: ") and "pairs took" in s or s.startswith("gpuPrep:") or s.startswith("GPU coverage") or s.startswith("gpuSparse: the sparse statistic") or s.startswith("gpuSparse: cross-term"):
                print(f"- {ds} {cfg}: {s[:230]}")
