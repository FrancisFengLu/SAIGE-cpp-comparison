#!/usr/bin/env python3
"""Compare two step-2 single-variant output sets pair for pair (binary traits
included) -- the exactness gate of /opt/saige/data/bingpu_test (README.md §7).

A and B are run directories: <dir>/out/ holds <trait>.txt (text) and/or
<trait>.txt.sgs (+ <first>.txt.markers.sgs; outputFormat: sgs, fp64 columns),
<dir>/routes/<trait>.route is the optional route dump
(SAIGE_STEP2_ROUTE_DUMP: {u8 route, f64 gateP} per (marker, trait) in input
marker order, 0xFF = no output row). A bare directory of *.txt files also works.
When both text and sgs exist for a side, sgs is used for the numbers (the
fp64 BETA / SE / Tstat / var / AF / counts) and text only for what sgs lacks.

Rows are aligned by (trait, MarkerID); rows present on one side only are
reported, never silently dropped. Per field: max |d|, max relative |d|, pairs
differing. p.value and p.value.NA: max |d(-log10 p)| (log-domain strings like
"1.0E-350" are handled), pairs crossing 5e-8 and 1e-5 in each direction.
Routing: per route bit (needSPA, needFirth, needFast, ER, SPAconv, Firth,
FirthConv) the count on each side and the pair-for-pair mismatches; without
route files, Is.SPA and "SPA changed p" (p.value != p.value.NA) from the output.
--coverage prints one side's route counts only (the self-gate of a reference).

usage: s2ref_cmp.py A B [--tol-logp 1e-8] [--tol-rel 0] [--traits t1,t2]
       s2ref_cmp.py --coverage A
Exit 1 when exact columns differ, routing differs, rows are missing on one
side, or a p-value tolerance is exceeded.
"""
import argparse, glob, math, os, struct, sys
import numpy as np

RT = np.dtype([("r", "u1"), ("p", "<f8")])
BITS = [("needSPA", 1), ("needFirth", 2), ("needFast", 4), ("ER", 8),
        ("SPAconv", 16), ("Firth", 32), ("FirthConv", 64)]
EXACT = ["CHR", "POS", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate",
         "imputationInfo", "N", "N_case", "N_ctrl", "N_case_hom", "N_case_het",
         "N_ctrl_hom", "N_ctrl_het", "Is.SPA"]
NUM = ["BETA", "SE", "Tstat", "var", "AF_case", "AF_ctrl", "p.value", "p.value.NA",
       "BETA_c", "SE_c", "Tstat_c", "var_c", "p.value_c", "p.value.NA_c"]
PCOLS = ["p.value", "p.value.NA", "p.value_c", "p.value.NA_c"]

# ----------------------------------------------------------------- sgs reader
NAME = {1: "BETA", 2: "SE", 3: "Tstat", 4: "var", 5: "p.value", 6: "p.value.NA", 7: "Is.SPA",
        8: "BETA_c", 9: "SE_c", 10: "Tstat_c", 11: "var_c", 12: "p.value_c", 13: "p.value.NA_c",
        14: "AF_case", 15: "AF_ctrl", 16: "N_case", 17: "N_ctrl",
        18: "N_case_hom", 19: "N_case_het", 20: "N_ctrl_hom", 21: "N_ctrl_het", 22: "N"}
F64 = {1, 2, 3, 4, 8, 9, 10, 11, 14, 15, 18, 19, 20, 21}
PV = {5, 6, 12, 13}
U32 = {16, 17, 22}
U8 = {7}


class R:
    def __init__(s, b): s.b = b; s.i = 0
    def take(s, n): q = s.b[s.i:s.i + n]; s.i += n; return q
    def u8(s): return s.take(1)[0]
    def u32(s): return struct.unpack("<I", s.take(4))[0]
    def st(s): n = s.u32(); return s.take(n).decode()
    def sstr(s):
        n = s.u8()
        if n == 255: n = s.u32()
        return s.take(n).decode()
    def peek_u32(s): return struct.unpack("<I", s.b[s.i:s.i + 4])[0]


def rd_f64(r, n):
    e = r.u8()
    if e == 0: return np.frombuffer(r.take(8 * n), "<f8").copy()
    if e == 1: return np.full(n, struct.unpack("<d", r.take(8))[0])
    if e == 3: return np.frombuffer(r.take(4 * n), "<f4").astype("f8")
    if e == 4: return np.full(n, struct.unpack("<f", r.take(4))[0], dtype="f8")
    raise RuntimeError("f64 enc %d" % e)


def rd_pod(r, n, w, dt):
    e = r.u8()
    if e == 0: return np.frombuffer(r.take(w * n), dt).copy()
    if e == 1: return np.full(n, np.frombuffer(r.take(w), dt)[0])
    raise RuntimeError("pod enc %d" % e)


def rd_str(r, n):
    e = r.u8()
    if e == 1: return [r.sstr()] * n
    if e == 0: return [r.sstr() for _ in range(n)]
    raise RuntimeError("str enc %d" % e)


def rd_pval(r, n):
    """returns the p-value strings as the text writer printed them"""
    e = r.u8()
    if e == 2: d = np.frombuffer(r.take(8 * n), "<f8").copy()
    elif e == 5: d = np.frombuffer(r.take(4 * n), "<f4").astype("f8")
    else: raise RuntimeError("pval enc %d" % e)
    out = ["%.6E" % v for v in d]
    for _ in range(r.u32()):
        i = r.u32(); out[i] = r.sstr()
    return out


def read_sgs_markers(path):
    r = R(open(path, "rb").read())
    assert r.take(8) == b"SAIGESGM", path
    r.u32(); flags = r.u32()
    info = "imputationInfo" if flags & 1 else "MissingRate"
    cols = {k: [] for k in ("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", info)}
    while r.peek_u32() == 0x314B4C42:
        r.u32(); n = r.u32()
        for k in ("CHR", "POS", "MarkerID", "Allele1", "Allele2"): cols[k] += rd_str(r, n)
        for k in ("AC_Allele2", "AF_Allele2", info): cols[k].append(rd_f64(r, n))
    assert r.peek_u32() == 0x21444E45, "no END! in " + path
    for k in ("AC_Allele2", "AF_Allele2", info): cols[k] = np.concatenate(cols[k]) if cols[k] else np.zeros(0)
    return cols, info


def read_sgs_trait(path, markers, info):
    """-> dict column -> np.array over the PRESENT rows (as the text file has them)"""
    r = R(open(path, "rb").read())
    assert r.take(8) == b"SAIGESGT", path
    r.u32(); r.u32(); r.st(); name = r.st(); r.st(); r.u8(); r.u8(); r.st(); r.st()
    ncol = r.u32(); codes = [r.u8() for _ in range(ncol)]
    out = {NAME[c]: [] for c in codes}
    for k in ("AC_Allele2", "AF_Allele2", info): out[k] = []
    present_all = []
    off = 0
    while r.peek_u32() == 0x314B4C42:
        r.u32(); n = r.u32(); bf = r.u32()
        present = np.frombuffer(r.take(n), "u1").astype(bool) if bf & 1 else np.ones(n, bool)
        ac = rd_f64(r, n) if bf & 2 else markers["AC_Allele2"][off:off + n]
        af = rd_f64(r, n) if bf & 4 else markers["AF_Allele2"][off:off + n]
        mi = rd_f64(r, n) if bf & 8 else markers[info][off:off + n]
        out["AC_Allele2"].append(ac[present]); out["AF_Allele2"].append(af[present]); out[info].append(mi[present])
        for c in codes:
            if c in F64: v = rd_f64(r, n)
            elif c in PV: v = np.array(rd_pval(r, n), dtype=object)
            elif c in U32: v = rd_pod(r, n, 4, "<u4")
            elif c in U8: v = rd_pod(r, n, 1, "u1")
            else: raise RuntimeError("col %d" % c)
            out[NAME[c]].append(v[present])
        present_all.append(present); off += n
    assert r.peek_u32() == 0x21444E45, "no END! in " + path
    present_all = np.concatenate(present_all) if present_all else np.zeros(0, bool)
    res = {k: np.concatenate(v) for k, v in out.items()}
    for k in ("CHR", "POS", "MarkerID", "Allele1", "Allele2"):
        res[k] = np.array(markers[k], dtype=object)[present_all]
    if "Is.SPA" in res: res["Is.SPA"] = np.where(res["Is.SPA"] != 0, "true", "false").astype(object)
    return name, res


# ----------------------------------------------------------------- text reader
def read_text(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f]
    cols = {}
    for k, h in enumerate(head):
        v = [row[k] for row in rows]
        if h in ("CHR", "POS", "MarkerID", "Allele1", "Allele2", "Is.SPA") or h in PCOLS:
            cols[h] = np.array(v, dtype=object)
        else:
            cols[h] = np.array([float(x) if x not in ("NA", "") else np.nan for x in v])
    return cols


# ----------------------------------------------------------------- run loader
def load_run(d, traits=None):
    out = d if os.path.isdir(os.path.join(d, "out")) is False else os.path.join(d, "out")
    sgs_traits = sorted(p for p in glob.glob(os.path.join(out, "*.txt.sgs")) if not p.endswith(".markers.sgs"))
    txt = sorted(glob.glob(os.path.join(out, "*.txt")))
    runs = {}
    if sgs_traits:
        mk = glob.glob(os.path.join(out, "*.markers.sgs"))
        assert len(mk) == 1, "expected one *.markers.sgs in %s" % out
        markers, info = read_sgs_markers(mk[0])
        for p in sgs_traits:
            t = os.path.basename(p)[:-8]
            name, cols = read_sgs_trait(p, markers, info)
            runs[t] = cols
    for p in txt:
        t = os.path.basename(p)[:-4]
        if t not in runs: runs[t] = read_text(p)
    if traits: runs = {t: v for t, v in runs.items() if t in traits}
    routes = {}
    rd = os.path.join(d, "routes")
    if os.path.isdir(rd):
        for p in glob.glob(os.path.join(rd, "*.route")):
            routes[os.path.basename(p)[:-6]] = np.fromfile(p, dtype=RT)
    src = "sgs" if sgs_traits else "text"
    return runs, routes, src


# ----------------------------------------------------------------- p helpers
def p_to_neglog10(s):
    """'1.234560E-05' -> 4.9...; log-domain underflow forms 'x.yE-NNN' parse the same way;
    'NA' -> nan. Values are clamped so a printed 0 maps to +inf."""
    if s is None or s == "NA" or s == "": return float("nan")
    try:
        v = float(s)
    except ValueError:
        # '%.1fE%d' underflow form may carry exponents a double cannot (e.g. 1.0E-350)
        m, e = s.upper().split("E")
        return -(math.log10(float(m)) + int(e))
    if v > 0: return -math.log10(v)
    if v == 0:
        m, e = s.upper().split("E") if "E" in s.upper() else (s, "0")
        return float("inf")
    return float("nan")


def p_float(s):
    try: return float(s)
    except (ValueError, TypeError):
        try:
            m, e = s.upper().split("E"); return float(m) * 10.0 ** int(e)
        except Exception: return float("nan")


def route_counts(rt):
    """counts per bit over rows that have an output (route != 0xFF)"""
    has = rt["r"] != 0xFF
    r = rt["r"][has]
    c = {name: int((r & bit != 0).sum()) for name, bit in BITS}
    c["rows"] = int(has.sum()); c["dropped"] = int((~has).sum())
    c["SPAconv&Firth"] = int(((r & 16 != 0) & (r & 32 != 0)).sum())
    c["ER&SPAconv"] = int(((r & 8 != 0) & (r & 16 != 0)).sum())
    return c


def fmt_counts(c):
    return ("rows %d (dropped %d)  needSPA %d  needFirth %d  needFast %d  ER(MAC<=cut) %d  "
            "SPAconv %d  Firth %d  FirthConv %d  [ER&SPAconv %d]" %
            (c["rows"], c["dropped"], c["needSPA"], c["needFirth"], c["needFast"], c["ER"],
             c["SPAconv"], c["Firth"], c["FirthConv"], c["ER&SPAconv"]))


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("A"); ap.add_argument("B", nargs="?")
    ap.add_argument("--coverage", action="store_true", help="print route coverage of A only")
    ap.add_argument("--tol-logp", type=float, default=1e-8, help="max allowed |d(-log10 p)| (default 1e-8)")
    ap.add_argument("--tol-rel", type=float, default=0.0, help="max allowed relative |d| on BETA/SE/Tstat/var (default 0 = exact)")
    ap.add_argument("--traits", help="comma list of trait names to restrict to")
    ap.add_argument("--max-print", type=int, default=5)
    a = ap.parse_args()
    traits = set(a.traits.split(",")) if a.traits else None

    A, RA, srcA = load_run(a.A, traits)
    if a.coverage or not a.B:
        tot = None
        for t in sorted(A):
            if t in RA:
                c = route_counts(RA[t]); print("%-10s %s" % (t, fmt_counts(c)))
                tot = c if tot is None else {k: tot[k] + c[k] for k in c}
            else:
                sp = A[t].get("Is.SPA"); pv, pn = A[t].get("p.value"), A[t].get("p.value.NA")
                print("%-10s rows %d  Is.SPA true %d  SPA changed p %d  (no route file)" %
                      (t, len(A[t]["MarkerID"]), int((sp == "true").sum()) if sp is not None else -1,
                       int((pv != pn).sum()) if pv is not None else -1))
        if tot: print("%-10s %s" % ("TOTAL", fmt_counts(tot)))
        sys.exit(0)

    B, RB, srcB = load_run(a.B, traits)
    print("A: %s (%s, %d traits)   B: %s (%s, %d traits)" % (a.A, srcA, len(A), a.B, srcB, len(B)))
    if srcA != srcB:
        # One side is text, whose float columns are printed with "%.6g" (6
        # significant digits); the sgs side holds the fp64 values. Compare at
        # the text's precision by passing the sgs side through the same format.
        q6 = np.vectorize(lambda v: float("%.6g" % v) if not np.isnan(v) else v, otypes=[float])
        for runs in ((A,) if srcA == "sgs" else (B,)):
            for t in runs:
                for h, v in runs[t].items():
                    if isinstance(v, np.ndarray) and v.dtype.kind == "f": runs[t][h] = q6(v)
        print("mixed sources: fp64 columns of the sgs side rounded to the text's %.6g before comparing")
    bad = 0
    only = sorted(set(A) ^ set(B))
    if only: print("traits on one side only: %s" % only); bad += 1
    agg = {}; dlog = {}; cross = {}
    exact_mis = {}; nrows = 0; miss_rows = 0
    rbits_a = {n: 0 for n, _ in BITS}; rbits_b = dict(rbits_a); rbits_mis = dict(rbits_a); rroute_pairs = 0; rbyte_mis = 0
    isspa_mis = spachg_mis = 0; isspa_a = isspa_b = 0; chg_a = chg_b = 0
    gate = {"n": 0, "rel": 0.0, "dlog": 0.0, "cross": {5e-8: [0, 0], 1e-5: [0, 0]}}
    printed = 0
    for t in sorted(set(A) & set(B)):
        ca, cb = A[t], B[t]
        ida, idb = ca["MarkerID"], cb["MarkerID"]
        pa = {m: i for i, m in enumerate(ida)}; pb = {m: i for i, m in enumerate(idb)}
        common = [m for m in ida if m in pb]
        oa = [m for m in ida if m not in pb]; ob = [m for m in idb if m not in pa]
        if oa or ob:
            miss_rows += len(oa) + len(ob); bad += 1
            print("  %s: rows only in A: %d (e.g. %s), only in B: %d (e.g. %s)" % (t, len(oa), oa[:3], len(ob), ob[:3]))
        ia = np.array([pa[m] for m in common], dtype=int); ib = np.array([pb[m] for m in common], dtype=int)
        nrows += len(common)
        if list(ida) != list(idb) and not (oa or ob):
            print("  %s: same row set, different order" % t)
        # exact columns
        for h in EXACT:
            if h not in ca or h not in cb: continue
            x, y = ca[h][ia], cb[h][ib]
            if x.dtype == object or y.dtype == object:
                ne = np.array([str(u) != str(v) for u, v in zip(x, y)])
            else:
                ne = ~((x == y) | (np.isnan(x) & np.isnan(y)))
            n = int(ne.sum())
            if n:
                exact_mis[h] = exact_mis.get(h, 0) + n
                if printed < a.max_print:
                    k = int(np.nonzero(ne)[0][0]); printed += 1
                    print("  EXACT %s %s %s: A %s  B %s" % (t, common[k], h, x[k], y[k]))
        # numeric columns
        for h in NUM:
            if h not in ca or h not in cb: continue
            if h in PCOLS:
                sa, sb = ca[h][ia], cb[h][ib]
                la = np.array([p_to_neglog10(s) for s in sa]); lb = np.array([p_to_neglog10(s) for s in sb])
                va = np.array([p_float(s) for s in sa]); vb = np.array([p_float(s) for s in sb])
                both_nan = np.isnan(la) & np.isnan(lb)
                one_nan = np.isnan(la) != np.isnan(lb)
                with np.errstate(invalid="ignore"):
                    d = np.where(both_nan, 0.0, np.abs(la - lb))
                    d = np.where(np.isinf(la) & np.isinf(lb), 0.0, d)
                d[one_nan] = np.inf
                m = dlog.setdefault(h, [0.0, 0])
                m[0] = max(m[0], float(np.nanmax(d)) if d.size else 0.0); m[1] += int((d > 0).sum())
                cr = cross.setdefault(h, {5e-8: [0, 0], 1e-5: [0, 0]})
                for thr in cr:
                    cr[thr][0] += int(((va < thr) & ~(vb < thr)).sum())   # significant in A only
                    cr[thr][1] += int((~(va < thr) & (vb < thr)).sum())   # significant in B only
                if (d > 0).any() and printed < a.max_print:
                    k = int(np.argmax(d)); printed += 1
                    print("  %s %s %s: A %s  B %s  |d log10| %.3e" % (h, t, common[k], sa[k], sb[k], d[k]))
                x, y = va, vb
            else:
                x, y = ca[h][ia].astype(float), cb[h][ib].astype(float)
            both = np.isnan(x) & np.isnan(y)
            d = np.where(both, 0.0, np.abs(x - y))
            den = np.maximum(np.abs(x), np.abs(y))
            with np.errstate(invalid="ignore", divide="ignore"):
                rel = np.where(den > 0, d / np.where(den > 0, den, 1.0), 0.0)
            m = agg.setdefault(h, [0.0, 0.0, 0])
            m[0] = max(m[0], float(np.nanmax(d)) if d.size else 0.0)
            m[1] = max(m[1], float(np.nanmax(rel)) if rel.size else 0.0)
            m[2] += int((d > 0).sum())
        # routing from the outputs
        if "Is.SPA" in ca and "Is.SPA" in cb:
            sa, sb = ca["Is.SPA"][ia], cb["Is.SPA"][ib]
            isspa_a += int((sa == "true").sum()); isspa_b += int((sb == "true").sum())
            isspa_mis += int((sa != sb).sum())
        if "p.value.NA" in ca and "p.value.NA" in cb:
            c1 = ca["p.value"][ia] != ca["p.value.NA"][ia]; c2 = cb["p.value"][ib] != cb["p.value.NA"][ib]
            chg_a += int(c1.sum()); chg_b += int(c2.sum()); spachg_mis += int((c1 != c2).sum())
        # routing from the route dumps (input order, all markers)
        if t in RA and t in RB:
            ra, rb = RA[t], RB[t]
            if len(ra) != len(rb):
                print("  %s: route dump length %d vs %d" % (t, len(ra), len(rb))); bad += 1
            else:
                rroute_pairs += len(ra)
                xa, xb = ra["r"], rb["r"]
                ne = xa != xb; rbyte_mis += int(ne.sum())
                if ne.any() and printed < a.max_print:
                    k = int(np.nonzero(ne)[0][0]); printed += 1
                    print("  ROUTE %s input row %d: A %#04x B %#04x" % (t, k, int(xa[k]), int(xb[k])))
                va_ = xa != 0xFF; vb_ = xb != 0xFF
                for n, bit in BITS:
                    ba = (xa & bit != 0) & va_; bb = (xb & bit != 0) & vb_
                    rbits_a[n] += int(ba.sum()); rbits_b[n] += int(bb.sum()); rbits_mis[n] += int((ba != bb).sum())
                ga, gb = ra["p"], rb["p"]
                m = ~np.isnan(ga) & ~np.isnan(gb)
                if (np.isnan(ga) != np.isnan(gb)).any():
                    print("  %s: gateP presence differs on %d pairs" % (t, int((np.isnan(ga) != np.isnan(gb)).sum()))); bad += 1
                if m.any():
                    x, y = ga[m], gb[m]
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rel = np.abs(x - y) / np.maximum(np.maximum(np.abs(x), np.abs(y)), 1e-300)
                        lx = np.where(x > 0, -np.log10(np.where(x > 0, x, 1)), np.abs(x) / math.log(10))
                        ly = np.where(y > 0, -np.log10(np.where(y > 0, y, 1)), np.abs(y) / math.log(10))
                    gate["n"] += int(m.sum()); gate["rel"] = max(gate["rel"], float(rel.max()))
                    gate["dlog"] = max(gate["dlog"], float(np.abs(lx - ly).max()))
                    for thr in gate["cross"]:
                        gate["cross"][thr][0] += int(((x < thr) & ~(y < thr)).sum())
                        gate["cross"][thr][1] += int((~(x < thr) & (y < thr)).sum())

    print("rows compared: %d (traits %d); rows on one side only: %d" % (nrows, len(set(A) & set(B)), miss_rows))
    print("exact columns: " + ("all identical" if not exact_mis else ", ".join("%s %d differ" % kv for kv in exact_mis.items())))
    for h in NUM:
        if h not in agg: continue
        m = agg[h]
        line = "%-12s max|d| %.3e  max rel %.3e  pairs differing %d" % (h, m[0], m[1], m[2])
        if h in dlog:
            cr = cross[h]
            line += "  max|d(-log10 p)| %.3e  cross 5e-8 A-only/B-only %d/%d  cross 1e-5 %d/%d" % (
                dlog[h][0], cr[5e-8][0], cr[5e-8][1], cr[1e-5][0], cr[1e-5][1])
        print(line)
    print("routing (outputs): Is.SPA true A %d B %d, mismatch %d; SPA changed p A %d B %d, mismatch %d" %
          (isspa_a, isspa_b, isspa_mis, chg_a, chg_b, spachg_mis))
    if rroute_pairs:
        print("routing (route dumps, %d pairs incl. dropped markers): route bytes differ on %d pairs" % (rroute_pairs, rbyte_mis))
        for n, _ in BITS:
            print("   %-10s A %7d  B %7d  mismatch %d" % (n, rbits_a[n], rbits_b[n], rbits_mis[n]))
        if gate["n"]:
            print("   batch-kernel pre-SPA p (fp64, %d pairs): max rel %.3e  max|d(-log10 p)| %.3e  cross 5e-8 %d/%d  cross 1e-5 %d/%d" %
                  (gate["n"], gate["rel"], gate["dlog"], gate["cross"][5e-8][0], gate["cross"][5e-8][1],
                   gate["cross"][1e-5][0], gate["cross"][1e-5][1]))
    else:
        print("routing (route dumps): not available on both sides")
    if exact_mis: bad += 1
    if isspa_mis or spachg_mis or rbyte_mis: bad += 1
    for h in dlog:
        if dlog[h][0] > a.tol_logp: bad += 1
    if gate["n"] and gate["dlog"] > a.tol_logp: bad += 1
    for h in ("BETA", "SE", "Tstat", "var"):
        if h in agg and agg[h][1] > a.tol_rel: bad += 1
    print("RESULT: " + ("PASS" if bad == 0 else "FAIL"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
