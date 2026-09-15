#!/usr/bin/env python3
"""Tolerance comparator for saige-null (step 1) outputs.

Compares ONE trait's outputs between a reference run and a test run, item by
item, against the thresholds in thresholds.yaml, and exits 0 (pass) / 1 (fail)
/ 2 (usage or internal error). Typical use: reference = the trait run alone
(P=1, "solo" layout), test = the same trait inside a multi-phenotype run
("multi" layout).

  compare.py --ref SOLO_DIR --test MULTI_DIR --trait y3 [--json out.json]
  compare.py --ref SOLO_A --test SOLO_B            # solo vs solo (run-to-run noise)
  compare.py orphans --dir MULTI_DIR --traits y1,y2 # files no trait claims

Run-directory layouts (out_prefix = DIR/m, out_prefix_vr = DIR/mvr, which is
what run_gate.sh writes):

  solo   DIR/m/...  DIR/m.grm_diag.txt  DIR/m_cov_*.csv  DIR/mvr.*
  multi  DIR/m/T/...  DIR/m/T.grm_diag.txt  DIR/m/T_cov_*.csv  DIR/mvr_T.*

The log (DIR.log by default) is required: converged / iterations / the LOCO
line are only printed there. A missing log is a failure unless --no-log.

Every file of the trait is claimed by a rule; a file with no numeric rule is
compared byte-for-byte. Missing or extra files fail. Parse errors fail.
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_THRESHOLDS = os.path.join(HERE, "thresholds.yaml")

# ----------------------------------------------------------------------------
# thresholds
# ----------------------------------------------------------------------------


def load_thresholds(path):
    import yaml
    with open(path) as f:
        t = yaml.safe_load(f)
    if not isinstance(t, dict) or "items" not in t:
        raise ValueError(f"{path}: expected a mapping with an 'items' key")
    return t


class Rules:
    def __init__(self, t):
        self.t = t
        self.defaults = t.get("defaults", {}) or {}

    def get(self, key):
        items = self.t["items"]
        if key not in items:
            raise KeyError(f"thresholds: no rule for item '{key}'")
        r = items[key]
        if isinstance(r, str):
            r = {"metric": r}
        out = dict(self.defaults)
        out.update(r)
        return out


# ----------------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------------

_ARMA_TYPES = {
    "FN004": np.float32, "FN008": np.float64,
    "IS001": np.int8, "IU001": np.uint8, "IS002": np.int16, "IU002": np.uint16,
    "IS004": np.int32, "IU004": np.uint32, "IS008": np.int64, "IU008": np.uint64,
}


def read_arma(path):
    """Armadillo arma_binary Mat/Col/Row: 'ARMA_MAT_BIN_<T>\\n<rows> <cols>\\n<raw column-major>'."""
    with open(path, "rb") as f:
        magic = f.readline().decode("ascii", "replace").strip()
        dims = f.readline().decode("ascii", "replace").split()
        blob = f.read()
    m = re.fullmatch(r"ARMA_MAT_BIN_([A-Z]{2}\d{3})", magic)
    if not m or m.group(1) not in _ARMA_TYPES:
        raise ValueError(f"unsupported arma header '{magic}'")
    if len(dims) != 2:
        raise ValueError(f"bad arma dims line {dims}")
    r, c = int(dims[0]), int(dims[1])
    dt = np.dtype(_ARMA_TYPES[m.group(1)])
    if len(blob) != r * c * dt.itemsize:
        raise ValueError(f"arma payload {len(blob)} bytes, expected {r}x{c}x{dt.itemsize}")
    a = np.frombuffer(blob, dtype=dt).reshape((c, r)).T  # column-major -> (rows, cols)
    return m.group(1), a


_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$|^[+-]?(inf|nan|Inf|NaN|INF|NAN)$")


def parse_num(s):
    s = s.strip()
    if s in ("NA", "null", ""):
        return float("nan")
    return float(s)


def read_table(path, header=True, sep=None):
    """-> (header list or None, list of rows (list of str))."""
    with open(path) as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f]
    lines = [l for l in lines if l != ""]
    if sep is None:
        sep = "," if (lines and "," in lines[0] and "\t" not in lines[0]) else None
    rows = [l.split(sep) if sep else l.split() for l in lines]
    if header:
        if not rows:
            return [], []
        return rows[0], rows[1:]
    return None, rows


def parse_logs(path):
    """All '== SAIGE Null Fit Completed ==' blocks -> list of dicts."""
    blocks = []
    cur = None
    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("== SAIGE Null Fit Completed =="):
                cur = {}
                blocks.append(cur)
                continue
            if cur is None:
                continue
            if line.startswith("Phenotype: "):
                cur["phenotype"] = line[len("Phenotype: "):].strip()
            elif line.startswith("Converged: "):
                cur["converged"] = line[len("Converged: "):].strip()
            elif line.startswith("Iterations: "):
                cur["iterations"] = line[len("Iterations: "):].strip()
            elif line.startswith("LOCO: "):
                cur["loco_line"] = " ".join(line.split())
                cur = None
    return blocks


# ----------------------------------------------------------------------------
# file discovery
# ----------------------------------------------------------------------------


def walk(root):
    out = []
    for d, _, fs in os.walk(root):
        for f in fs:
            out.append(os.path.relpath(os.path.join(d, f), root))
    return sorted(out)


def detect_layout(d, trait):
    if os.path.isfile(os.path.join(d, "m", "nullmodel.json")):
        return "solo"
    if trait and os.path.isfile(os.path.join(d, "m", trait, "nullmodel.json")):
        return "multi"
    if trait and os.path.isdir(os.path.join(d, "m", trait)):
        return "multi"
    if os.path.isdir(os.path.join(d, "m")) and not trait:
        return "solo"
    return None


def collect(d, layout, trait):
    """-> {logical_name: relative path} for one trait, and a list of info notes."""
    files = {}
    if not os.path.isdir(d):
        return files
    for rel in walk(d):
        logical = None
        if layout == "solo":
            if rel.startswith("m/"):
                logical = "model/" + rel[2:]
            elif rel == "m.grm_diag.txt":
                logical = "grm_diag.txt"
            elif rel.startswith("m_cov_"):
                logical = "cov_" + rel[len("m_cov_"):]
            elif rel.startswith("mvr."):
                logical = "vr." + rel[len("mvr."):]
            else:
                logical = "unclaimed/" + rel
        else:
            t = trait
            if rel.startswith(f"m/{t}/"):
                logical = "model/" + rel[len(f"m/{t}/"):]
            elif rel == f"m/{t}.grm_diag.txt":
                logical = "grm_diag.txt"
            elif rel.startswith(f"m/{t}_cov_") and "/" not in rel[len(f"m/{t}_cov_"):]:
                logical = "cov_" + rel[len(f"m/{t}_cov_"):]
            elif rel.startswith(f"mvr_{t}.") and "/" not in rel:
                logical = "vr." + rel[len(f"mvr_{t}."):]
            else:
                continue  # another trait's file (see `orphans`)
        m = re.fullmatch(r"vr\.(\d+)markers\.SAIGE\.results\.txt", logical)
        if m:
            files["vr.markers.SAIGE.results.txt"] = rel
            files["__vr_n_markers__"] = int(m.group(1))
            continue
        files[logical] = rel
    return files


def claimed_by(rel, t):
    return (rel.startswith(f"m/{t}/") or rel == f"m/{t}.grm_diag.txt"
            or (rel.startswith(f"m/{t}_cov_") and "/" not in rel[len(f"m/{t}_cov_"):])
            or (rel.startswith(f"mvr_{t}.") and "/" not in rel))


# ----------------------------------------------------------------------------
# comparison core
# ----------------------------------------------------------------------------


class Report:
    def __init__(self):
        self.items = []

    def add(self, file, item, rule_key, rule, passed, dev=None, n=None, note=None, **extra):
        tol = rule.get("tol") if rule else None
        ratio = None
        if dev is not None and tol not in (None, 0):
            ratio = (dev / tol) if math.isfinite(dev) else float("inf")
        elif dev is not None and tol == 0:
            ratio = 0.0 if dev == 0 else float("inf")
        if rule and rule.get("metric") == "exact":
            ratio = 0.0 if passed else float("inf")
        rec = dict(file=file, item=item, rule=rule_key, metric=(rule or {}).get("metric"),
                   tol=tol, dev=dev, ratio=ratio, n=n, passed=bool(passed))
        if note:
            rec["note"] = note
        rec.update(extra)
        self.items.append(rec)
        return rec

    def fail(self, file, item, note):
        return self.add(file, item, "structure", {"metric": "exact"}, False, note=note)


def _fmt(x):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isnan(x):
            return "nan"
        if math.isinf(x):
            return "inf"
        return f"{x:.3g}"
    return str(x)


def num_compare(rep, file, item, rule_key, rule, a, b, labels=None, floor_vec=None):
    """a, b: arrays of equal shape (a = ref, b = test). 2-D arrays get a per-column floor.
    floor_vec (optional, same size as a.ravel()) replaces the RMS-based floor."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    col_rms = None
    if a.ndim == 2 and a.shape == b.shape and a.size:
        with np.errstate(invalid="ignore"):
            col_rms = np.sqrt(np.nanmean(np.where(np.isfinite(a), a, np.nan) ** 2, axis=0))
        col_rms = np.nan_to_num(np.broadcast_to(col_rms, a.shape)).ravel()
    a = a.ravel()
    b = b.ravel()
    metric = rule["metric"]
    if a.shape != b.shape:
        return rep.add(file, item, rule_key, rule, False, dev=float("inf"),
                       note=f"shape ref {a.shape} test {b.shape}")
    n = a.size
    if n == 0:
        return rep.add(file, item, rule_key, rule, True, dev=0.0, n=0)
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if np.any(nan_a != nan_b):
        i = int(np.nonzero(nan_a != nan_b)[0][0])
        return rep.add(file, item, rule_key, rule, False, dev=float("inf"), n=n,
                       note="NaN pattern differs", worst_index=_lab(labels, i),
                       ref_value=_v(a[i]), test_value=_v(b[i]))
    ok = ~nan_a
    inf_a, inf_b = np.isinf(a) & ok, np.isinf(b) & ok
    if np.any((inf_a | inf_b) & (a != b)):
        i = int(np.nonzero((inf_a | inf_b) & (a != b))[0][0])
        return rep.add(file, item, rule_key, rule, False, dev=float("inf"), n=n,
                       note="Inf differs", worst_index=_lab(labels, i),
                       ref_value=_v(a[i]), test_value=_v(b[i]))
    fin = ok & ~inf_a & ~inf_b
    diff = np.zeros(n)
    diff[fin] = np.abs(a[fin] - b[fin])
    if metric == "exact":
        bad = np.nonzero(fin & (a != b))[0]
        dev = float(diff.max()) if n else 0.0
        if bad.size:
            i = int(bad[np.argmax(diff[bad])])
            return rep.add(file, item, rule_key, rule, False, dev=dev, n=n,
                           note=f"{bad.size} of {n} values differ", worst_index=_lab(labels, i),
                           ref_value=_v(a[i]), test_value=_v(b[i]))
        return rep.add(file, item, rule_key, rule, True, dev=0.0, n=n)
    if metric == "rel":
        if floor_vec is not None:
            floor = np.maximum(np.asarray(floor_vec, dtype=np.float64).ravel(), float(rule.get("abs_floor", 0.0)))
        elif col_rms is not None:
            floor = np.maximum(float(rule.get("abs_floor", 0.0)), float(rule.get("floor_frac", 0.0)) * col_rms)
        else:
            rms = float(np.sqrt(np.mean(a[fin] ** 2))) if np.any(fin) else 0.0
            floor = max(float(rule.get("abs_floor", 0.0)), float(rule.get("floor_frac", 0.0)) * rms)
        slack = np.zeros(n)
        if rule.get("print_sig"):
            # values read back from text written with `print_sig` significant digits (C++ stream
            # precision; trailing zeros are dropped, so the string itself understates it): two
            # correctly rounded prints of true values x, y differ by at most |x - y| + 1 unit in the
            # last place, so that unit is subtracted before the relative test
            mag = np.maximum(np.abs(a), np.abs(b))
            with np.errstate(divide="ignore"):
                e = np.floor(np.log10(np.where(mag > 0, mag, 1.0)))
            slack = np.where(mag > 0, 10.0 ** (e - (int(rule["print_sig"]) - 1)), 0.0)
        scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), floor)
        scale[scale == 0] = 1.0
        d = np.zeros(n)
        d[fin] = np.maximum(diff[fin] - slack[fin], 0.0) / scale[fin]
        i = int(np.argmax(d))
        dev = float(d[i])
        tol = float(rule["tol"])
        nbad = int(np.sum(d > tol))
        return rep.add(file, item, rule_key, rule, dev <= tol, dev=dev, n=n,
                       note=(f"{nbad} of {n} values over tol" if nbad else None),
                       worst_index=_lab(labels, i), ref_value=_v(a[i]), test_value=_v(b[i]),
                       floor=float(np.broadcast_to(floor, (n,))[i]), max_abs_diff=float(diff.max()))
    if metric == "log10p":
        pa = np.clip(np.where(fin, a, 1.0), 1e-300, None)
        pb = np.clip(np.where(fin, b, 1.0), 1e-300, None)
        d = np.abs(np.log10(pa) - np.log10(pb))
        i = int(np.argmax(d))
        dev = float(d[i])
        tol = float(rule["tol"])
        return rep.add(file, item, rule_key, rule, dev <= tol, dev=dev, n=n,
                       worst_index=_lab(labels, i), ref_value=_v(a[i]), test_value=_v(b[i]))
    raise ValueError(f"unknown metric {metric}")


def _lab(labels, i):
    return labels[i] if labels is not None else i


def _v(x):
    x = float(x)
    return x if math.isfinite(x) else str(x)


def exact_compare(rep, file, item, a, b, note_fn=None):
    ok = a == b
    note = None
    if not ok:
        note = note_fn() if note_fn else f"ref {str(a)[:120]} | test {str(b)[:120]}"
    return rep.add(file, item, "exact", {"metric": "exact"}, ok, note=note)


# ---- per-file comparators --------------------------------------------------

ARMA_ITEM = {
    "mu": "arma.mu", "res": "arma.res", "V": "arma.V", "S_a": "arma.S_a",
    "XV": "arma.XV", "XVX": "arma.XVX", "XVX_inv": "arma.XVX_inv", "XXVX_inv": "arma.XXVX_inv",
    "XVX_inv_XV": "arma.XVX_inv_XV", "X": "arma.X", "y": "arma.y", "offset": "arma.offset",
}


def score_floor(rule, ref_root, fa, logical_dir):
    """S_a = colSums(X * res) sits near 0 at convergence, so a floor from S_a itself is noise.
    Scale it by the score's absolute mass instead: floor_j = sa_floor_frac * sum_i |X_ij res_i|
    (X from model/X.arma, res from the same directory as S_a)."""
    frac = rule.get("sa_floor_frac")
    xl, rl = "model/X.arma", logical_dir + "res.arma"
    if frac is None or xl not in fa or rl not in fa:
        return None
    _, X = read_arma(os.path.join(ref_root, fa[xl]))
    _, r = read_arma(os.path.join(ref_root, fa[rl]))
    if X.shape[0] != r.shape[0]:
        return None
    return float(frac) * np.abs(X * r.reshape(-1, 1)).sum(axis=0)


def cmp_arma(rep, rules, logical, pa, pb, ref_root=None, fa=None):
    base = os.path.basename(logical)[:-len(".arma")]
    ta, a = read_arma(pa)
    tb, b = read_arma(pb)
    if ta != tb:
        return rep.fail(logical, "type", f"arma type ref {ta} test {tb}")
    if a.shape != b.shape:
        return rep.fail(logical, "shape", f"shape ref {a.shape} test {b.shape}")
    if base in ARMA_ITEM:
        key = ARMA_ITEM[base]
        rule = rules.get(key)
        fv = None
        if base == "S_a" and fa is not None:
            fv = score_floor(rule, ref_root, fa, logical[:-len("S_a.arma")])
            if fv is not None and fv.size != a.size:
                fv = None
        return num_compare(rep, logical, base, key, rule, a, b, floor_vec=fv)
    return rep.fail(logical, "rule", f"no rule for arma file '{base}' (add one to thresholds.yaml)")


def cmp_sparse_grm(rep, rules, dir_logical, pa_loc, pa_val, pb_loc, pb_val):
    """sparseGRM_locationMat (2 x nnz) + sparseGRM_valueVec (nnz): compared as a set of triplets."""
    key = "arma.sparseGRM"
    rule = rules.get(key)
    _, la = read_arma(pa_loc)
    _, va = read_arma(pa_val)
    _, lb = read_arma(pb_loc)
    _, vb = read_arma(pb_val)
    file = dir_logical + "sparseGRM_{locationMat,valueVec}.arma"
    if la.shape[0] != 2 or lb.shape[0] != 2 or la.shape[1] != va.size or lb.shape[1] != vb.size:
        return rep.fail(file, "shape", f"loc ref {la.shape} val {va.shape}; loc test {lb.shape} val {vb.shape}")
    if va.size != vb.size:
        return rep.fail(file, "nnz", f"nnz ref {va.size} test {vb.size}")
    same_order = bool(np.array_equal(la, lb))
    oa = np.lexsort((va.ravel(), la[0], la[1]))
    ob = np.lexsort((vb.ravel(), lb[0], lb[1]))
    ra, rb = la[:, oa], lb[:, ob]
    if not np.array_equal(ra, rb):
        k = int(np.nonzero(np.any(ra != rb, axis=0))[0][0])
        return rep.fail(file, "positions",
                        f"triplet {k} after sorting: ref (r{ra[0, k]},c{ra[1, k]}) test (r{rb[0, k]},c{rb[1, k]})")
    rep.add(file, "positions", "exact", {"metric": "exact"}, True, n=int(va.size),
            note=("entry order identical" if same_order else "entry order differs (compared as a set)"))
    return num_compare(rep, file, "values", key, rule, va.ravel()[oa], vb.ravel()[ob])


NULLMODEL_EXACT = ["trait", "traitType", "n", "p", "loco", "lowmem_loco", "loco_chroms", "SPA_Cutoff",
                   "impute_method", "flagSparseGRM", "isFastTest", "isnoadjCov", "isCondition",
                   "is_Firth_beta", "pCutoffforFirth"]


def cmp_nullmodel(rep, rules, logical, pa, pb):
    A = json.load(open(pa))
    B = json.load(open(pb))
    ka, kb = set(A), set(B)
    if ka != kb:
        rep.fail(logical, "keys", f"ref-only {sorted(ka - kb)} test-only {sorted(kb - ka)}")
    known = set(NULLMODEL_EXACT) | {"theta", "alpha", "sampleIDs"}
    for k in NULLMODEL_EXACT:
        if k in A or k in B:
            exact_compare(rep, logical, k, A.get(k), B.get(k))
    for k in sorted((ka & kb) - known):
        exact_compare(rep, logical, k, A.get(k), B.get(k))
    ia, ib = A.get("sampleIDs"), B.get("sampleIDs")

    def ids_note():
        if not isinstance(ia, list) or not isinstance(ib, list):
            return "sampleIDs missing"
        if len(ia) != len(ib):
            return f"length ref {len(ia)} test {len(ib)}"
        j = next(i for i in range(len(ia)) if ia[i] != ib[i])
        return f"first difference at {j}: ref {ia[j]} test {ib[j]}"
    exact_compare(rep, logical, "sampleIDs", ia, ib, ids_note)
    for k, key in (("theta", "nullmodel.theta"), ("alpha", "nullmodel.alpha")):
        va, vb = A.get(k), B.get(k)
        if not isinstance(va, list) or not isinstance(vb, list) or len(va) != len(vb):
            rep.fail(logical, k, f"ref {va} test {vb}")
            continue
        num_compare(rep, logical, k, key, rules.get(key), va, vb,
                    labels=[f"{k}[{i}]" for i in range(len(va))])


def cmp_obj_nok(rep, rules, logical, pa, pb, ref_root=None, fa=None):
    A = json.load(open(pa))
    B = json.load(open(pb))

    def packs(J):
        out = {"baseline": J.get("baseline")}
        for e in J.get("loco", []) or []:
            out[f"chr{e.get('chrom')}"] = e.get("pack")
        return out
    pa_, pb_ = packs(A), packs(B)
    exact_compare(rep, logical, "packs", sorted(pa_), sorted(pb_))
    for name in sorted(set(pa_) & set(pb_)):
        x, y = pa_[name] or {}, pb_[name] or {}
        exact_compare(rep, logical, f"{name}.n,p", (x.get("n"), x.get("p")), (y.get("n"), y.get("p")))
        for fld in ("V", "S_a", "XVX", "XVX_inv"):
            u, v = x.get(fld), y.get(fld)
            if not isinstance(u, list) or not isinstance(v, list) or len(u) != len(v):
                rep.fail(logical, f"{name}.{fld}", "missing or length differs")
                continue
            key = f"arma.{fld}"
            rule = rules.get(key)
            fv = None
            if fld == "S_a" and fa is not None:
                d = "model/" if name == "baseline" else f"model/{name}/"
                fv = score_floor(rule, ref_root, fa, d)
                # obj_noK packs are built on the design the GLMM saw; with covariate_offset that is
                # the intercept-only X, i.e. the first column of the restored X.arma
                fv = fv[:len(u)] if (fv is not None and len(u) <= fv.size) else None
            num_compare(rep, logical, f"{name}.{fld}", key, rule, u, v, floor_vec=fv)


def cmp_grm_diag(rep, rules, logical, pa, pb):
    _, ra = read_table(pa, header=False)
    _, rb = read_table(pb, header=False)
    ta = [r[0] for r in ra]
    tb = [r[0] for r in rb]
    if len(ta) != len(tb):
        return rep.fail(logical, "length", f"ref {len(ta)} test {len(tb)} lines")
    return num_compare(rep, logical, "diag", "grm_diag", rules.get("grm_diag"),
                       [parse_num(x) for x in ta], [parse_num(x) for x in tb])


def cmp_cov_csv(rep, rules, logical, pa, pb):
    kind = logical[len("cov_"):-len(".csv")]
    key = f"cov.{kind}" if f"cov.{kind}" in rules.t["items"] else "cov.other"
    ha, ra = read_table(pa, sep=",")
    hb, rb = read_table(pb, sep=",")
    exact_compare(rep, logical, "header", ha, hb)
    if ha != hb:
        return
    if len(ra) != len(rb):
        return rep.fail(logical, "rows", f"ref {len(ra)} test {len(rb)}")
    for j, col in enumerate(ha):
        ca = [r[j] if j < len(r) else "" for r in ra]
        cb = [r[j] if j < len(r) else "" for r in rb]
        if all(_NUM_RE.match(x.strip()) or x.strip() == "NA" for x in ca + cb):
            num_compare(rep, logical, col, key, rules.get(key),
                        [parse_num(x) for x in ca], [parse_num(x) for x in cb])
        else:
            exact_compare(rep, logical, col, ca, cb,
                          lambda: _first_str_diff(ca, cb))


def _first_str_diff(ca, cb):
    j = next((i for i in range(min(len(ca), len(cb))) if ca[i] != cb[i]), None)
    return f"first difference at row {j}: ref {ca[j]} test {cb[j]}" if j is not None else "lengths differ"


def cmp_varratio(rep, rules, logical, pa, pb):
    _, ra = read_table(pa, header=False)
    _, rb = read_table(pb, header=False)
    sa = [(r[1], r[2]) for r in ra]
    sb = [(r[1], r[2]) for r in rb]
    exact_compare(rep, logical, "rows(type,bin)", sa, sb)
    if sa != sb:
        return
    ta, tb = [r[0] for r in ra], [r[0] for r in rb]
    num_compare(rep, logical, "ratio", "vr.ratio", rules.get("vr.ratio"),
                [parse_num(x) for x in ta], [parse_num(x) for x in tb],
                labels=[f"{t}/bin{k}" for t, k in sa])


PVAL_RE = re.compile(r"^(p\.value|p\.value\.NA|pval|p_value|P|p)$", re.I)


def cmp_vr_markers(rep, rules, logical, pa, pb, n_a, n_b):
    exact_compare(rep, logical, "n_markers_tested(file name)", n_a, n_b)
    ha, ra = read_table(pa)
    hb, rb = read_table(pb)
    exact_compare(rep, logical, "header", ha, hb)
    if ha != hb:
        return
    if len(ra) != len(rb):
        return rep.fail(logical, "rows", f"ref {len(ra)} test {len(rb)}")
    for j, col in enumerate(ha):
        ca = [r[j] for r in ra]
        cb = [r[j] for r in rb]
        if PVAL_RE.match(col):
            key = "vr_markers.pvalue"
        else:
            key = f"vr_markers.{col}"
            if key not in rules.t["items"]:
                key = "vr_markers.other"
        rule = rules.get(key)
        if all(_NUM_RE.match(x.strip()) or x.strip() == "NA" for x in ca + cb):
            num_compare(rep, logical, col, key, rule, [parse_num(x) for x in ca],
                        [parse_num(x) for x in cb],
                        labels=[f"row{i}" for i in range(len(ca))])
        else:
            exact_compare(rep, logical, col, ca, cb, lambda: _first_str_diff(ca, cb))


def cmp_bytes(rep, logical, pa, pb):
    same = open(pa, "rb").read() == open(pb, "rb").read()
    rep.add(logical, "bytes", "exact", {"metric": "exact"}, same,
            note=None if same else "no numeric rule for this file; bytes differ")


# ----------------------------------------------------------------------------


def compare_trait(ref, test, trait=None, ref_layout="auto", test_layout="auto",
                  ref_log=None, test_log=None, use_log=True, thresholds=DEFAULT_THRESHOLDS):
    rules = Rules(load_thresholds(thresholds))
    rep = Report()
    info = {}
    la = detect_layout(ref, trait) if ref_layout == "auto" else ref_layout
    lb = detect_layout(test, trait) if test_layout == "auto" else test_layout
    if la is None:
        la = "multi" if trait else "solo"
        rep.fail("(ref)", "layout", f"no outputs found under {ref} for trait {trait}")
    if lb is None:
        lb = "multi" if trait else "solo"
        rep.fail("(test)", "layout", f"no outputs found under {test} for trait {trait}")
    if "multi" in (la, lb) and not trait:
        raise SystemExit("compare.py: --trait is required when a side uses the multi layout")
    info.update(ref_layout=la, test_layout=lb)
    fa = collect(ref, la, trait)
    fb = collect(test, lb, trait)
    n_a = fa.pop("__vr_n_markers__", None)
    n_b = fb.pop("__vr_n_markers__", None)
    only_a = sorted(set(fa) - set(fb))
    only_b = sorted(set(fb) - set(fa))
    rep.add("(file set)", "files", "exact", {"metric": "exact"}, not only_a and not only_b,
            n=len(set(fa) | set(fb)),
            note=(f"ref-only {only_a} test-only {only_b}" if (only_a or only_b) else None))
    common = sorted(set(fa) & set(fb))
    if not common:
        rep.fail("(file set)", "files", "no common files to compare")

    # sparse GRM pairs are compared together
    handled = set()
    for lg in common:
        if lg.endswith("sparseGRM_locationMat.arma"):
            d = lg[:-len("sparseGRM_locationMat.arma")]
            lv = d + "sparseGRM_valueVec.arma"
            if lv in fa and lv in fb:
                try:
                    cmp_sparse_grm(rep, rules, d, os.path.join(ref, fa[lg]), os.path.join(ref, fa[lv]),
                                   os.path.join(test, fb[lg]), os.path.join(test, fb[lv]))
                except Exception as e:  # noqa: BLE001
                    rep.fail(lg, "parse", f"{type(e).__name__}: {e}")
                handled |= {lg, lv}

    for lg in common:
        if lg in handled:
            continue
        pa, pb = os.path.join(ref, fa[lg]), os.path.join(test, fb[lg])
        try:
            if lg.startswith("model/") and lg.endswith("nullmodel.json"):
                cmp_nullmodel(rep, rules, lg, pa, pb)
            elif lg.startswith("model/") and lg.endswith("obj_noK.json"):
                cmp_obj_nok(rep, rules, lg, pa, pb, ref, fa)
            elif lg.startswith("model/") and lg.endswith(".arma"):
                cmp_arma(rep, rules, lg, pa, pb, ref, fa)
            elif lg == "grm_diag.txt":
                cmp_grm_diag(rep, rules, lg, pa, pb)
            elif lg.startswith("cov_") and lg.endswith(".csv"):
                cmp_cov_csv(rep, rules, lg, pa, pb)
            elif lg == "vr.varianceRatio.txt":
                cmp_varratio(rep, rules, lg, pa, pb)
            elif lg == "vr.markers.SAIGE.results.txt":
                cmp_vr_markers(rep, rules, lg, pa, pb, n_a, n_b)
            else:
                cmp_bytes(rep, lg, pa, pb)
        except Exception as e:  # noqa: BLE001
            rep.fail(lg, "parse", f"{type(e).__name__}: {e}")

    # LOCO: chr directories must match (already covered by the file set, but say it)
    ca = sorted({lg.split("/")[1] for lg in fa if re.match(r"model/chr[^/]+/", lg)})
    cb = sorted({lg.split("/")[1] for lg in fb if re.match(r"model/chr[^/]+/", lg)})
    info["loco_chr_dirs"] = ca
    if ca or cb:
        exact_compare(rep, "model/chr*", "chr directories", ca, cb)

    if use_log:
        ref_log = ref_log or (ref.rstrip("/") + ".log")
        test_log = test_log or (test.rstrip("/") + ".log")
        blk = {}
        for side, path, layout in (("ref", ref_log, la), ("test", test_log, lb)):
            if not os.path.isfile(path):
                rep.fail("(log)", side, f"log not found: {path} (pass --no-log to skip)")
                continue
            bl = parse_logs(path)
            if layout == "solo":
                if len(bl) != 1:
                    rep.fail("(log)", side, f"solo log has {len(bl)} completed-fit blocks, expected 1")
                    continue
                b = bl[0]
                if trait and b.get("phenotype") not in (None, trait):
                    rep.fail("(log)", side, f"solo log is for phenotype {b.get('phenotype')}")
                    continue
            else:
                mb = [b for b in bl if b.get("phenotype") == trait]
                if len(mb) != 1:
                    rep.fail("(log)", side, f"{len(mb)} completed-fit blocks for phenotype {trait}")
                    continue
                b = mb[0]
            blk[side] = b
        if "ref" in blk and "test" in blk:
            for fld in ("converged", "iterations", "loco_line"):
                va, vb = blk["ref"].get(fld), blk["test"].get(fld)
                key = f"log.{fld}"
                rule = rules.get(key)
                ok = va is not None and va == vb
                rep.add("(log)", fld, key, rule, ok,
                        note=None if ok else f"ref {va!r} test {vb!r}", ref_value=va, test_value=vb)
            info["iterations"] = blk["ref"].get("iterations")
            info["converged"] = blk["ref"].get("converged")

    failed = [it for it in rep.items if not it["passed"]]
    graded = [it for it in rep.items if it.get("ratio") is not None]
    worst = max(graded, key=lambda it: it["ratio"]) if graded else None
    num_items = [it for it in rep.items if it.get("metric") in ("rel", "log10p")]
    worst_num = max(num_items, key=lambda it: it["ratio"]) if num_items else None
    return dict(trait=trait, ref=ref, test=test, thresholds=os.path.abspath(thresholds),
                passed=not failed, n_items=len(rep.items), n_failed=len(failed),
                worst=worst, worst_numeric=worst_num, info=info, items=rep.items)


def describe(it):
    s = f"{it['file']} :: {it['item']}"
    if it.get("metric") in ("rel", "log10p"):
        s += f"  {it['metric']} dev={_fmt(it['dev'])} tol={_fmt(it['tol'])}"
        if it.get("worst_index") is not None:
            s += f" at {it['worst_index']} (ref {_fmt(it.get('ref_value'))}, test {_fmt(it.get('test_value'))})"
    elif not it["passed"]:
        s += "  exact: differs"
    if it.get("note"):
        s += f"  [{it['note']}]"
    return s


def summary_text(res, verbose=False):
    lines = []
    head = "PASS" if res["passed"] else "FAIL"
    t = res["trait"] or "(solo)"
    w = res["worst_numeric"]
    ws = f"worst numeric: {describe(w)}" if w else "worst numeric: -"
    lines.append(f"[compare] {head} trait={t}: {res['n_items']} items, {res['n_failed']} failed; {ws}")
    for it in res["items"]:
        if not it["passed"]:
            lines.append("  FAIL " + describe(it))
    if verbose:
        for it in res["items"]:
            if it["passed"]:
                lines.append("  ok   " + describe(it))
    return "\n".join(lines)


def orphans(d, traits):
    out = []
    for rel in walk(d):
        if not any(claimed_by(rel, t) for t in traits):
            out.append(rel)
    return out


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "orphans":
        ap = argparse.ArgumentParser(prog="compare.py orphans")
        ap.add_argument("--dir", required=True)
        ap.add_argument("--traits", required=True)
        a = ap.parse_args(argv[1:])
        o = orphans(a.dir, [t for t in a.traits.split(",") if t])
        print(f"[orphans] {'PASS' if not o else 'FAIL'} {len(o)} file(s) not claimed by any trait"
              + (f": {o[:20]}" if o else ""))
        return 0 if not o else 1
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference run directory (usually the solo run)")
    ap.add_argument("--test", required=True, help="test run directory (usually the multi-phenotype run)")
    ap.add_argument("--trait", help="phenotype column; required for the multi layout")
    ap.add_argument("--ref-layout", default="auto", choices=["auto", "solo", "multi"])
    ap.add_argument("--test-layout", default="auto", choices=["auto", "solo", "multi"])
    ap.add_argument("--ref-log")
    ap.add_argument("--test-log")
    ap.add_argument("--no-log", action="store_true", help="do not require/compare logs")
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS)
    ap.add_argument("--json", help="write the full report here")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args(argv)
    try:
        res = compare_trait(a.ref, a.test, a.trait, a.ref_layout, a.test_layout,
                            a.ref_log, a.test_log, not a.no_log, a.thresholds)
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001
        print(f"[compare] ERROR {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    if a.json:
        os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1, default=str)
    if not a.quiet:
        print(summary_text(res, a.verbose))
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
