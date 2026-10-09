#!/usr/bin/env python3
"""gate_table.py [datasets...]: every C++ path x configuration x trait against R, and the multi-trait /
GPU paths against the single-trait path. Prints Markdown tables (worst over the traits of a dataset)."""
import sys, os, json, subprocess, glob
B = os.environ.get("RDEF", "/opt/saige/logs/rdefaults")
CMP = os.path.join(B, "cmp_r.py")
DS = sys.argv[1:] or ["audit", "bvs", "qt12"]
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
        w["maxrel"][k] = max(r["maxrel"][k] for r in rs)
        w["beyondPrint"][k] = sum(r["beyondPrint"][k] for r in rs)
    for k in ["cross5e8", "cross1e5", "isSPAmismatch", "nSPA_A", "firthMismatch", "firthA", "firthB"]:
        w[k] = sum(r[k] for r in rs)
    return w
detail = []
print("## C++ vs R SAIGE 1.5.2 (worst over the dataset's traits; rel = max relative difference)\n")
print(hdr)
for ds in DS:
    for cfg in cfgs[ds]:
        for path in ["single", "multi", "gpu"]:
            rs = []
            for t in traits[ds]:
                if cfg.endswith("F") and t.startswith("q"): continue
                r = cmp(f"{B}/cpp/{ds}/{cfg}/{path}/{t}.txt", f"{B}/R/{ds}/{cfg}/{t}.txt", f"{ds} {cfg} {path} {t}")
                rs.append(r); detail.append(r)
            w = worst(rs)
            if w: w["label"] = f"{ds} {cfg} {path} vs R"
            print(fmt(w))
print("\n## multi-trait and GPU paths vs the single-trait path (same binary)\n")
print(hdr)
for ds in DS:
    for cfg in cfgs[ds]:
        for path in ["multi", "gpu"]:
            rs = []
            for t in traits[ds]:
                if cfg.endswith("F") and t.startswith("q"): continue
                rs.append(cmp(f"{B}/cpp/{ds}/{cfg}/{path}/{t}.txt", f"{B}/cpp/{ds}/{cfg}/single/{t}.txt", f"{ds} {cfg} {path} {t}"))
            w = worst(rs)
            if w: w["label"] = f"{ds} {cfg} {path} vs single"
            print(fmt(w))
if "--detail" in os.environ.get("GATE_OPTS", ""):
    print("\n## per trait\n"); print(hdr)
    for r in detail: print(fmt(r))
