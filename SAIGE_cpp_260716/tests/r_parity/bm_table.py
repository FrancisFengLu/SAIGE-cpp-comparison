#!/usr/bin/env python3
"""bm_table.py: the bm gate table (bingpu_test bm, 8 binary traits with 4 missing-phenotype patterns).
Per configuration (worst over the 8 traits): the new binary's GPU multi-trait run vs its CPU single-trait
runs, and vs the base (origin/rdefaults) binary's GPU multi-trait run (host tail); plus the route-byte
comparison of the two GPU runs and the device-stats log lines."""
import os, sys, json, subprocess, glob
B = "/opt/saige/logs/ownstats/bm"
CMP = "/opt/saige/logs/ownstats/rparity/cmp_r.py"
traits = ["bm%d" % k for k in range(1, 9)]
cfgs = ["def", "defF", "adj", "adjF"]
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
print("## bm (N = 50,000, 26,060 markers, 8 binary traits, 4 sample lists): new GPU multi-trait vs CPU single-trait and vs the base GPU run (host tail)\n")
print(hdr)
for cfg in cfgs:
    for lab, other in [("vs single", lambda t: f"{B}/{cfg}/single/{t}.txt"), ("vs base gpu", lambda t: f"{B}/{cfg}/gpu_base/{t}.txt")]:
        rs = [cmp(f"{B}/{cfg}/gpu_new/{t}.txt", other(t), f"bm {cfg} gpu {t}") for t in traits]
        w = worst(rs)
        if w: w["label"] = f"bm {cfg} gpu_new {lab}"
        print(fmt(w))
print("\n## route bytes (SAIGE_STEP2_ROUTE_DUMP) of the two GPU runs, and the device-stats lines\n")
for cfg in cfgs:
    # each record is { u8 route, f64 gateP }: the route byte is compared exactly, the
    # gate p (the batch kernel's pre-SPA p: device erfc vs host Boost) by relative difference
    import struct, math
    nd = 0; nb = 0; gmax = 0.0; gdiff = 0; gn = 0
    for t in traits:
        a = f"{B}/{cfg}/gpu_new/routes/{t}.route"; b = f"{B}/{cfg}/gpu_base/routes/{t}.route"
        if os.path.exists(a) and os.path.exists(b):
            da = open(a, "rb").read(); db = open(b, "rb").read()
            n = min(len(da), len(db)) // 9
            nb += n; nd += abs(len(da) - len(db)) // 9
            for i in range(n):
                ra, rb = da[9*i], db[9*i]
                if ra != rb: nd += 1
                pa = struct.unpack("<d", da[9*i+1:9*i+9])[0]; pb = struct.unpack("<d", db[9*i+1:9*i+9])[0]
                if math.isnan(pa) and math.isnan(pb): continue
                gn += 1
                if pa != pb:
                    gdiff += 1
                    r = abs(pa - pb) / max(abs(pa), abs(pb))
                    if r > gmax: gmax = r
    lines = {}
    for lab in ["gpu_new", "gpu_base"]:
        L = open(f"{B}/{cfg}/{lab}/log.txt").read().splitlines()
        lines[lab] = [l.strip() for l in L if l.strip().startswith("gpuDeviceStats: ") and ("pairs took" in l or "off (" in l or "own sample list" in l)]
    print(f"- {cfg}: route bytes compared {nb}, differing {nd}; gate p compared {gn}, differing {gdiff}, max rel {gmax:.1e}")
    for lab in ["gpu_new", "gpu_base"]:
        for l in lines[lab]: print(f"  - {lab}: {l[:200]}")
