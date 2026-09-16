#!/usr/bin/env python3
"""Generate the phenotype files + manifests the multi-trait tolerance gate runs on.

  make_cases.py [--workdir /opt/saige/logs/mt_gate] [--scale mid|small|all] [--cases a,b] [--force]

Writes WORKDIR/data/<scale>/<case>.pheno.txt and <case>.manifest.json. The run
configs and gate-time checks live in cases/<case>.yaml; this script only makes
the data those configs point at, verifies that the data really contain the
scenario the case is for (QC / fill-value boundary markers, subset sizes, ...),
and records what the binary must print if it takes that scenario (predicted
GRM marker counts, first-5-sample genotype counts). A case whose data do not
meet its minimum coverage makes this script exit non-zero.

Everything is seeded (SEEDS below); rerunning reproduces the files byte for
byte. External inputs (genotypes, mid.mp32/block16/indep16 phenotypes, the
sparse GRMs) are recorded with md5 in each manifest.
"""
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import genoqc as gq  # noqa: E402

GENERATOR_VERSION = 2
DATA = "/opt/saige/data"
SPDATA = "/opt/saige/logs/missing_mt/step1/data"   # sparse GRM + sparse-signal phenotypes
SEEDS = dict(small_base=51001, block_small=51002, indep_small=51003, tiny=51004, qc=51005,
             fill=51006, weak=51007, quant=51008, loco=51009, sparse=51010, vrcate=51011)

QC_TRAIT_SEED_BUMP = {("small", "qA3"): 1, ("mid", "qA3"): 1}   # (scale, trait) -> int, see case_qc

SCALES = {
    "mid": dict(plink=f"{DATA}/mid", loco_bim_src=f"{DATA}/mid.bim", n_chr=5,
                x_src=f"{DATA}/mid.mp32.pheno.txt", base=f"{DATA}/mid.mp32.pheno.txt",
                block=f"{DATA}/mid.block16.pheno.txt", indep=f"{DATA}/mid.indep16.pheno.txt",
                score=f"{DATA}/mid.gs.sscore",
                sgrm=f"{SPDATA}/mid.sgrm2", sp_pheno=f"{SPDATA}/mid.sp.indep.pheno.txt"),
    "small": dict(plink=f"{DATA}/small", loco_plink=f"{DATA}/smallloco", n_chr=5,
                  x_src=f"{DATA}/small.pheno.txt", base=None, block=None, indep=None,
                  score=f"{DATA}/small.gs.sscore",
                  sgrm=f"{SPDATA}/small.sgrm2", sp_pheno=f"{SPDATA}/small.sp.block.pheno.txt"),
}
ALL_CASES = ["nomiss", "block", "indep", "tiny", "qc", "fill", "weak", "quant", "loco",
             "sparse", "vrcate"]


def md5(path, big_ok=False):
    st = os.stat(path)
    if st.st_size > 64 << 20 and not big_ok:
        return f"size={st.st_size},mtime={int(st.st_mtime)}"
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    return hdr, {h: [r[i] for r in rows] for i, h in enumerate(hdr)}


def set_sha1(idx):
    return hashlib.sha1(np.sort(np.asarray(idx, np.int64)).tobytes()).hexdigest()[:16]


class Ctx:
    """Per-scale shared state: genotypes, full-sample QC, IIDs, covariates, base traits."""

    def __init__(self, scale, workdir):
        self.scale = scale
        self.cfg = SCALES[scale]
        self.out = os.path.join(workdir, "data", scale)
        os.makedirs(self.out, exist_ok=True)
        self.workdir = workdir
        self._G = None
        self._full = None
        self.iids = [l.split()[1] for l in open(self.cfg["plink"] + ".fam")]
        self.N = len(self.iids)
        hdr, cols = read_tsv(self.cfg["x_src"])
        assert cols["IID"] == self.iids, "covariate file not in FAM order"
        self.x1, self.x2 = cols["x1"], cols["x2"]
        self.inputs = {self.cfg["plink"] + ext: md5(self.cfg["plink"] + ext) for ext in (".bed", ".bim", ".fam")}
        self.inputs[self.cfg["x_src"]] = md5(self.cfg["x_src"])

    @property
    def G(self):
        if self._G is None:
            print(f"[{self.scale}] decoding genotypes ...", flush=True)
            self._G = gq.decode(self.cfg["plink"])
            self.M = self._G.shape[0]
            vr = os.path.join(self.workdir, "data", "vrdraw")
            if not os.path.exists(vr):
                cxx = os.environ.get("CXX", "g++")
                subprocess.run([cxx, "-O2", "-o", vr, os.path.join(HERE, "vrdraw.cpp")], check=True)
            self.drawn = gq.drawn_mask(self.M, vr)
        return self._G

    @property
    def full(self):
        if self._full is None:
            self._full = gq.qc(self.G, np.arange(self.N), self.drawn)
        return self._full

    def qc(self, idx):
        return gq.qc(self.G, idx, self.drawn)

    def first5(self, idx, st):
        """Genotype counts (0/1/2, missing filled) over GRM markers for the set's first 5 samples,
        as printed by output_grm_diagonal ('Sample k: 0=a, 1=b, 2=c')."""
        idx = np.sort(np.asarray(idx))[:5]
        J = np.nonzero(st["passQC"])[0]
        out = []
        for i in idx:
            g = self.G[J, i].astype(np.int64)
            m = g == 3
            g[m] = st["fill"][J][m]
            out.append([int((g == 0).sum()), int((g == 1).sum()), int((g == 2).sum())])
        return out

    def base_traits(self):
        """Binary traits y1..y32 on everyone: mid.mp32 for mid, simulated for small."""
        if self.cfg["base"]:
            hdr, cols = read_tsv(self.cfg["base"])
            assert cols["IID"] == self.iids
            self.inputs[self.cfg["base"]] = md5(self.cfg["base"])
            return {k: v for k, v in cols.items() if k.startswith("y")}
        rng = np.random.default_rng(SEEDS["small_base"])
        g = self.score()
        x1 = np.array([float(v) for v in self.x1])
        x2 = np.array([float(v) for v in self.x2])
        out = {}
        for k in range(1, 33):
            w = 0.2 + 0.02 * k
            gk = math.sqrt(w) * g + math.sqrt(1 - w) * rng.standard_normal(self.N)
            eta = -0.3 + 0.9 * gk + 0.3 * x1 - 0.2 * x2
            y = (rng.random(self.N) < 1 / (1 + np.exp(-eta))).astype(int)
            out[f"y{k}"] = [str(v) for v in y]
        return out

    def score(self):
        sc = {}
        for l in open(self.cfg["score"]):
            if l.startswith("#"):
                continue
            a = l.split()
            sc[a[0]] = float(a[1])
        g = np.array([sc[i] for i in self.iids])
        self.inputs[self.cfg["score"]] = md5(self.cfg["score"])
        return (g - g.mean()) / g.std()


def strong_binary(ctx, rng, eff=2.5):
    """Binary trait with a large genetic effect (liability h2 ~ 0.65 from the plink2 genetic
    score). Small subsets of mp32-like traits fit tau1 = 0, and with tau1 = 0 the GRM never
    enters the fit (VR = 1, mu = GLM): a GRM bug on such a trait is invisible except in
    grm_diag. These traits keep tau1 > 0 down to n ~ 600 on mid (checked: nondegenerate)."""
    g = ctx.score()
    x1 = np.array([float(v) for v in ctx.x1])
    x2 = np.array([float(v) for v in ctx.x2])
    eta = -0.3 + eff * g + 0.3 * x1 - 0.2 * x2
    return [str(v) for v in (rng.random(ctx.N) < 1 / (1 + np.exp(-eta))).astype(int)]


def mask(vals, keep_idx):
    keep = np.zeros(len(vals), bool)
    keep[np.asarray(keep_idx, np.int64)] = True
    return [v if keep[i] else "NA" for i, v in enumerate(vals)]


def nonmissing(vals):
    return np.array([i for i, v in enumerate(vals) if v != "NA"], np.int64)


def write_case(ctx, case, traits, extra):
    """traits: list of (name, values[str], source note). Writes pheno + manifest."""
    path = os.path.join(ctx.out, f"{case}.pheno.txt")
    names = [t[0] for t in traits]
    with open(path, "w") as f:
        f.write("\t".join(["IID"] + names + ["x1", "x2"]) + "\n")
        for i, iid in enumerate(ctx.iids):
            f.write("\t".join([iid] + [t[1][i] for t in traits] + [ctx.x1[i], ctx.x2[i]]) + "\n")
    tr = {}
    sets = {}
    for name, vals, src in traits:
        idx = nonmissing(vals)
        h = set_sha1(idx)
        sets.setdefault(h, []).append(name)
        tr[name] = dict(n=int(len(idx)), set=h, source=src)
    man = dict(case=case, scale=ctx.scale, generator_version=GENERATOR_VERSION,
               seed=SEEDS.get(case), pheno=path, pheno_md5=md5(path), traits=tr,
               distinct_sample_sets=len(sets), sets={h: v for h, v in sets.items()},
               inputs=dict(sorted(ctx.inputs.items())))
    man.update(extra)
    with open(os.path.join(ctx.out, f"{case}.manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    return man


def predictions(ctx, traits, which=None):
    """grm marker count + first-5 counts per trait (on the dense-GRM genotype file)."""
    pred = {"grm_markers": {}, "first5_counts": {}}
    cache = {}
    for name, vals, _ in traits:
        if which is not None and name not in which:
            continue
        idx = nonmissing(vals)
        h = set_sha1(idx)
        if h not in cache:
            st = ctx.qc(idx)
            cache[h] = (int(st["passQC"].sum()), ctx.first5(idx, st))
        pred["grm_markers"][name], pred["first5_counts"][name] = cache[h]
    return pred


def require(cond, msg, problems):
    print(("  ok   " if cond else "  MISS ") + msg, flush=True)
    if not cond:
        problems.append(msg)


# ----------------------------------------------------------------------------
# cases
# ----------------------------------------------------------------------------


def case_nomiss(ctx, problems):
    base = ctx.base_traits()
    traits = [(f"b{k}", base[f"y{k}"], f"base:y{k}") for k in range(1, 5)]
    require(all(len(nonmissing(v)) == ctx.N for _, v, _ in traits), "all 4 traits on every sample", problems)
    return write_case(ctx, "nomiss", traits, dict(predicted=predictions(ctx, traits)))


def case_block(ctx, problems):
    if ctx.cfg["block"]:
        hdr, cols = read_tsv(ctx.cfg["block"])
        assert cols["IID"] == ctx.iids
        ctx.inputs[ctx.cfg["block"]] = md5(ctx.cfg["block"])
        traits = [(c, cols[c], f"mid.block16:{c}") for c in ("y1", "y2", "y3", "y9", "y10", "y11")]
    else:
        base = ctx.base_traits()
        rng = np.random.default_rng(SEEDS["block_small"])
        panel = rng.choice(ctx.N, int(0.6 * ctx.N), replace=False)
        traits = [(f"y{k}", base[f"y{k}"], f"base:y{k}") for k in (1, 2, 3)]
        traits += [(f"y{k}", mask(base[f"y{k}"], panel), f"base:y{k} on 60% panel") for k in (9, 10, 11)]
    ns = sorted({len(nonmissing(v)) for _, v, _ in traits})
    require(len(ns) == 2 and ns[1] == ctx.N, f"two sample sets, n={ns}", problems)
    return write_case(ctx, "block", traits, dict(predicted=predictions(ctx, traits)))


def case_indep(ctx, problems):
    if ctx.cfg["indep"]:
        hdr, cols = read_tsv(ctx.cfg["indep"])
        assert cols["IID"] == ctx.iids
        ctx.inputs[ctx.cfg["indep"]] = md5(ctx.cfg["indep"])
        traits = [(c, cols[c], f"mid.indep16:{c}") for c in ("y1", "y2", "y3", "y4")]
    else:
        base = ctx.base_traits()
        rng = np.random.default_rng(SEEDS["indep_small"])
        traits = []
        for k in (1, 2, 3, 4):
            drop = rng.choice(ctx.N, int(0.05 * ctx.N), replace=False)
            keep = np.setdiff1d(np.arange(ctx.N), drop)
            traits.append((f"y{k}", mask(base[f"y{k}"], keep), f"base:y{k}, 5% dropped independently"))
    hs = {set_sha1(nonmissing(v)) for _, v, _ in traits}
    require(len(hs) == 4, f"4 distinct sample sets (got {len(hs)})", problems)
    return write_case(ctx, "indep", traits, dict(predicted=predictions(ctx, traits)))


def case_tiny(ctx, problems):
    base = ctx.base_traits()
    rng = np.random.default_rng(SEEDS["tiny"])
    s5 = rng.choice(ctx.N, int(0.05 * ctx.N), replace=False)
    s10 = rng.choice(ctx.N, int(0.10 * ctx.N), replace=False)
    traits = [("tF", base["y5"], "base:y5 all samples"),
              ("t5a", mask(strong_binary(ctx, rng), s5), "strong binary on a random 5%"),
              ("t5b", mask(strong_binary(ctx, rng), s5), "strong binary on the same 5%"),
              ("t10", mask(strong_binary(ctx, rng), s10), "strong binary on an independent random 10%")]
    ns = [len(nonmissing(v)) for _, v, _ in traits]
    require(ns == [ctx.N, len(s5), len(s5), len(s10)], f"n = {ns}", problems)
    pred = predictions(ctx, traits)
    fullc = pred["grm_markers"]["tF"]
    d = {k: pred["grm_markers"][k] - fullc for k in ("t5a", "t10")}
    require(all(v != 0 for v in d.values()), f"subset GRM marker counts differ from full: {d}", problems)
    return write_case(ctx, "tiny", traits, dict(predicted=pred))


def case_qc(ctx, problems):
    """Five traits on four subsets whose GRM marker lists differ from the one of the UNION of all
    the case's samples, in both directions and for both reasons (MAF, missing rate)."""
    G, full, N = ctx.G, ctx.full, ctx.N
    rng = np.random.default_rng(SEEDS["qc"])
    nA = int(0.02 * N)
    miss_cnt = full["miss"]
    cand = np.nonzero(full["passQC"])[0]
    order = cand[np.argsort(-miss_cnt[cand], kind="stable")]
    q4, q2a, q2b = int(order[0]), int(order[1]), int(order[2])
    missers = lambda j: np.nonzero(G[j] == 3)[0]
    B = missers(q4)                                      # (iv): every sample missing at q4
    pool = np.setdiff1d(np.arange(N), B)                 # A-sets never contain a q4 misser
    # (ii) A2: 17% of it missing at q2a, 17% at q2b
    k2 = int(math.ceil(0.17 * nA))
    ma = rng.permutation(np.intersect1d(missers(q2a), pool))[:k2]
    mb = rng.permutation(np.setdiff1d(np.intersect1d(missers(q2b), pool), ma))[:k2]
    rest = rng.permutation(np.setdiff1d(pool, np.concatenate([ma, mb])))
    A2 = np.concatenate([ma, mb, rest[:nA - len(ma) - len(mb)]])
    pool1 = np.setdiff1d(pool, A2)
    # (i) A1: non-carriers of Z1 (full MAF 2%-5%, GRM markers in full)
    z1c = np.nonzero(full["passQC"] & (full["maf"] >= 0.02) & (full["maf"] < 0.05))[0]
    Z1 = np.sort(rng.choice(z1c, 12, replace=False))
    mac_z1 = np.zeros(N, np.int64)
    for j in Z1:
        g = G[j].astype(np.int64)
        minor = np.where(g == 3, 0, g) if full["alt"][j] < 0.5 else np.where(g == 3, 0, 2 - g)
        mac_z1 += minor
    c1 = pool1[mac_z1[pool1] == 0]
    A1 = rng.choice(c1, nA, replace=False)
    pool3 = np.setdiff1d(pool1, A1)
    # (iii) A3: enriched for carriers of Z3 (full MAF 0.2%-0.8%, NOT GRM markers in full)
    z3c = np.nonzero(~full["passQC"] & ~full["passVR"] & (full["maf"] >= 0.002) & (full["maf"] < 0.008)
                     & (full["mr"] <= 0.15))[0]
    Z3 = np.sort(rng.choice(z3c, 10, replace=False))
    target = int(math.ceil(0.014 * 2 * nA))              # minor alleles wanted per Z3 marker
    chosen = []
    for j in Z3:
        g = G[j].astype(np.int64)
        minor = np.where(g == 3, 0, g) if full["alt"][j] < 0.5 else np.where(g == 3, 0, 2 - g)
        carriers = rng.permutation(pool3[minor[pool3] > 0])
        have = int(minor[np.array(chosen, np.int64)].sum()) if chosen else 0
        for c in carriers:
            if have >= target:
                break
            if c not in chosen:
                chosen.append(int(c)); have += int(minor[c])
    chosen = np.array(chosen, np.int64)
    rest3 = rng.permutation(np.setdiff1d(pool3, chosen))
    A3 = np.concatenate([chosen, rest3[:max(0, nA - len(chosen))]])

    def trait(k, name, idx, note):
        # one RNG stream per trait: at n = 140-1000 whether AI-REML ends at tau1 = 0 is mostly
        # sampling luck, so a trait that lands on 0 gets its own stream bumped (QC_TRAIT_SEED_BUMP)
        # without changing the others; the gate's nondegenerate check verifies the result
        bump = QC_TRAIT_SEED_BUMP.get((ctx.scale, name), 0)
        r = np.random.default_rng([SEEDS["qc"], k, bump])
        return (name, mask(strong_binary(ctx, r), idx), note + (f" (trait seed bump {bump})" if bump else ""))
    traits = [trait(1, "qA1", A1, "strong binary on A1 (Z1 non-carriers)"),
              trait(2, "qA1b", A1, "strong binary on A1"),
              trait(3, "qA2", A2, "strong binary on A2 (missers of q2a/q2b)"),
              trait(4, "qA3", A3, "strong binary on A3 (Z3 carriers)"),
              trait(5, "qB", B, "strong binary on B (every q4 misser)")]
    U = np.unique(np.concatenate([A1, A2, A3, B]))
    sU = ctx.qc(U)
    st = {"A1": ctx.qc(A1), "A2": ctx.qc(A2), "A3": ctx.qc(A3), "B": ctx.qc(B)}
    cls = {}
    for k, s in st.items():
        pu, px = sU["passQC"], s["passQC"]
        cls[k] = dict(
            n=int(s["n"]), grm_markers=int(px.sum()),
            i_passU_failX_maf=int((pu & ~px & (s["maf"] < 0.01)).sum()),
            ii_passU_failX_miss=int((pu & ~px & (s["mr"] > 0.15)).sum()),
            iii_failU_passX_maf=int((~pu & px & (sU["maf"] < 0.01)).sum()),
            iv_failU_passX_miss=int((~pu & px & (sU["mr"] > 0.15)).sum()),
            v_vrU_grmX=int((sU["passVR"] & px).sum()),
            passU_failX_total=int((pu & ~px).sum()), failU_passX_total=int((~pu & px).sum()))
    des = dict(
        q4=q4, q2=[q2a, q2b], Z1=Z1.tolist(), Z3=Z3.tolist(), union_n=int(len(U)),
        union_grm_markers=int(sU["passQC"].sum()),
        designed_i=int((sU["passQC"][Z1] & (st["A1"]["maf"][Z1] < 0.01)).sum()),
        designed_ii=int((sU["passQC"][[q2a, q2b]] & (st["A2"]["mr"][[q2a, q2b]] > 0.15)).sum()),
        designed_iii=int((~sU["passQC"][Z3] & (sU["maf"][Z3] < 0.01) & st["A3"]["passQC"][Z3]).sum()),
        designed_iv=int(sum(bool(~sU["passQC"][q4] and sU["mr"][q4] > 0.15 and st[k]["passQC"][q4])
                            for k in ("A1", "A2", "A3"))),
        q4_missrate_union=float(sU["mr"][q4]), q4_missrate_B=float(st["B"]["mr"][q4]),
        classes=cls)
    print(json.dumps({k: v for k, v in des.items() if k not in ("Z1", "Z3")}, indent=1))
    require(des["designed_i"] >= 5, f"(i) pass in union, MAF<0.01 in A1: {des['designed_i']} of 12 designed", problems)
    require(des["designed_ii"] >= 1, f"(ii) pass in union, missing>0.15 in A2: {des['designed_ii']} of 2", problems)
    require(des["designed_iii"] >= 3, f"(iii) MAF<0.01 in union, pass in A3: {des['designed_iii']} of 10", problems)
    require(des["designed_iv"] >= 1, f"(iv) missing>0.15 in union, pass in A1/A2/A3: {des['designed_iv']} of 3", problems)
    return write_case(ctx, "qc", traits, dict(design=des, predicted=predictions(ctx, traits)))


def case_fill(ctx, problems):
    """A 10% subset chosen so that many markers WITH missing genotypes in it get a different
    fill value round(2 f_pre) than on the union (= everyone, because tF is on all samples)."""
    G, full, N = ctx.G, ctx.full, ctx.N
    rng = np.random.default_rng(SEEDS["fill"])
    nS = int(0.10 * N)
    a = 2.0 * full["altpre"].astype(np.float64)
    w = 0.06
    near = full["passQC"] & ((np.abs(a - 0.5) < w) | (np.abs(a - 1.5) < w))
    J = np.nonzero(near)[0]
    if len(J) > 800:
        J = np.sort(rng.choice(J, 800, replace=False))
    b = np.where(np.abs(a[J] - 0.5) < w, 0.5, 1.5)
    dirn = np.where(a[J] >= b, -1.0, 1.0)                # push 2f across the boundary
    f = full["altpre"][J].astype(np.float64)
    sd = np.sqrt(np.maximum(2 * f * (1 - f), 1e-6))
    score = np.zeros(N)
    nmissJ = np.zeros(N)
    for k, j in enumerate(J):
        g = G[j].astype(np.float64)
        m = g == 3
        score += np.where(m, 0.0, dirn[k] * (g - a[j]) / sd[k])
        nmissJ += m
    score += 6.0 * nmissJ                                 # also prefer samples missing at J
    anchors = np.argsort(-nmissJ[:200], kind="stable")[:5]   # first 5 samples: most missing at J
    lo = int(anchors.max())
    cand = np.arange(lo + 1, N)
    top = cand[np.argsort(-(score[cand] + rng.gumbel(0, 1, len(cand))), kind="stable")[:nS - 5]]
    S = np.sort(np.concatenate([anchors, top]))
    sS = ctx.qc(S)
    flip = (sS["fill"] != full["fill"]) & sS["passQC"] & (sS["miss"] > 0)
    first5_S = ctx.first5(S, sS)
    # the same 5 samples with the UNION's fill value substituted (what a bug would print)
    stU = dict(sS)
    stU["fill"] = full["fill"]
    first5_Ufill = ctx.first5(S, stU)
    des = dict(nS=int(nS), candidate_markers=int(len(J)),
               flip_markers=int(flip.sum()),
               flip_markers_grm_in_both=int((flip & full["passQC"]).sum()),
               flip_cells=int(sS["miss"][flip].sum()),
               flip_u1_s0=int((flip & (full["fill"] == 1) & (sS["fill"] == 0)).sum()),
               flip_u0_s1=int((flip & (full["fill"] == 0) & (sS["fill"] == 1)).sum()),
               flip_u2_s1=int((flip & (full["fill"] == 2) & (sS["fill"] == 1)).sum()),
               flip_u1_s2=int((flip & (full["fill"] == 1) & (sS["fill"] == 2)).sum()),
               max_missrate_on_flipped=float(sS["mr"][flip].max()) if flip.any() else 0.0,
               first5_S=first5_S, first5_with_union_fill=first5_Ufill,
               first5_sensitive=int(sum(x != y for x, y in zip(first5_S, first5_Ufill))),
               grm_markers_S=int(sS["passQC"].sum()), grm_markers_union=int(full["passQC"].sum()))
    print(json.dumps(des, indent=1))
    need = 50 if ctx.scale == "mid" else 20
    require(des["flip_markers"] >= need, f"markers with missing genotypes whose fill differs subset vs union: "
            f"{des['flip_markers']} (>= {need}), {des['flip_cells']} genotype cells", problems)
    require(des["first5_sensitive"] >= 1, f"first-5-sample counts would change under the union fill for "
            f"{des['first5_sensitive']} of 5 samples", problems)
    base = ctx.base_traits()
    traits = [("fF", base["y16"], "base:y16 all samples"),
              ("fS", mask(strong_binary(ctx, rng), S), "strong binary on the fill-boundary 10%"),
              ("fS2", mask(strong_binary(ctx, rng), S), "strong binary on the same 10%")]
    return write_case(ctx, "fill", traits, dict(design=des, predicted=predictions(ctx, traits)))


def case_weak(ctx, problems):
    rng = np.random.default_rng(SEEDS["weak"])
    x1 = np.array([float(v) for v in ctx.x1])
    x2 = np.array([float(v) for v in ctx.x2])
    panel = rng.choice(ctx.N, int(0.6 * ctx.N), replace=False)
    traits = []
    for name in ("wF1", "wF2", "wS"):
        eta = -0.2 + 0.3 * x1 - 0.2 * x2                   # no genetic component at all
        y = (rng.random(ctx.N) < 1 / (1 + np.exp(-eta))).astype(int)
        v = [str(t) for t in y]
        if name == "wS":
            v = mask(v, panel)
        traits.append((name, v, "binary, no genetic effect" + (" on 60% panel" if name == "wS" else "")))
    return write_case(ctx, "weak", traits, dict(predicted=predictions(ctx, traits)))


def case_quant(ctx, problems):
    rng = np.random.default_rng(SEEDS["quant"])
    g = ctx.score()
    x1 = np.array([float(v) for v in ctx.x1])
    x2 = np.array([float(v) for v in ctx.x2])
    panel = rng.choice(ctx.N, int(0.6 * ctx.N), replace=False)
    traits = []
    for name in ("qtF1", "qtF2", "qtS"):
        e = rng.standard_normal(ctx.N)
        q = math.sqrt(0.6) * g + 0.3 * x1 + 0.2 * x2 + math.sqrt(0.4) * e
        v = [f"{t:.5f}" for t in q]
        if name == "qtS":
            v = mask(v, panel)
        traits.append((name, v, "quantitative, h2~0.6 from gs score" + (" on 60% panel" if name == "qtS" else "")))
    return write_case(ctx, "quant", traits, dict(predicted=predictions(ctx, traits)))


def loco_plink(ctx):
    if "loco_plink" in ctx.cfg:
        return ctx.cfg["loco_plink"]
    dst = os.path.join(ctx.out, "midloco")
    for ext in (".bed", ".fam"):
        if not os.path.lexists(dst + ext):
            os.symlink(ctx.cfg["plink"] + ext, dst + ext)
    n_chr = ctx.cfg["n_chr"]
    lines = open(ctx.cfg["loco_bim_src"]).read().splitlines()
    M = len(lines)
    per = M // n_chr
    with open(dst + ".bim", "w") as f:
        for i, l in enumerate(lines):
            a = l.split("\t")
            a[0] = str(min(i // per, n_chr - 1) + 1)
            f.write("\t".join(a) + "\n")
    return dst


def case_loco(ctx, problems):
    base = ctx.base_traits()
    rng = np.random.default_rng(SEEDS["loco"])
    panel = rng.choice(ctx.N, int(0.6 * ctx.N), replace=False)
    lp = loco_plink(ctx)
    chroms = sorted({l.split()[0] for l in open(lp + ".bim")}, key=int)
    require(len(chroms) >= 2, f"LOCO genotype file has {len(chroms)} chromosomes: {chroms}", problems)
    traits = [("lF1", base["y19"], "base:y19"), ("lF2", base["y20"], "base:y20"),
              ("lS", mask(base["y21"], panel), "base:y21 on 60% panel")]
    return write_case(ctx, "loco", traits, dict(loco_plink=lp, loco_chroms=[int(c) for c in chroms],
                                                predicted=predictions(ctx, traits)))


def case_sparse(ctx, problems):
    sp = ctx.cfg["sp_pheno"]
    hdr, cols = read_tsv(sp)
    assert cols["IID"] == ctx.iids
    for p in (sp, ctx.cfg["sgrm"] + ".mtx", ctx.cfg["sgrm"] + ".ids"):
        ctx.inputs[p] = md5(p)
    rng = np.random.default_rng(SEEDS["sparse"])
    panel = rng.choice(ctx.N, int(0.6 * ctx.N), replace=False)
    if ctx.scale == "mid":        # sb1..sb4, each 5% missing independently
        both = np.intersect1d(nonmissing(cols["sb1"]), nonmissing(cols["sb2"]))
        traits = [("sA", mask(cols["sb1"], both), "mid.sp.indep:sb1 on nonmiss(sb1)&nonmiss(sb2)"),
                  ("sB", mask(cols["sb2"], both), "mid.sp.indep:sb2 on the same set"),
                  ("sC", cols["sb3"], "mid.sp.indep:sb3 as is"),
                  ("sD", mask(cols["sb4"], np.intersect1d(panel, nonmissing(cols["sb4"]))),
                   "mid.sp.indep:sb4 on 60% panel")]
    else:                         # small.sp.block: sb1..sb4 full, sb5..sb8 on a 60% block
        drop = rng.choice(ctx.N, int(0.05 * ctx.N), replace=False)
        traits = [("sA", cols["sb1"], "small.sp.block:sb1"), ("sB", cols["sb2"], "small.sp.block:sb2"),
                  ("sC", cols["sb5"], "small.sp.block:sb5 (block)"),
                  ("sD", mask(cols["sb3"], np.setdiff1d(np.arange(ctx.N), drop)), "small.sp.block:sb3, 5% dropped")]
    return write_case(ctx, "sparse", traits, dict(sparse_grm=ctx.cfg["sgrm"] + ".mtx",
                                                  sparse_grm_ids=ctx.cfg["sgrm"] + ".ids",
                                                  predicted=predictions(ctx, traits)))


def case_vrcate(ctx, problems):
    """Categorical VR with two MAC bins that both get markers. The loader only keeps VR-eligible
    markers with MAC >= 20, so R's default bins (10.5, 20.5] / >20.5 leave bin 1 fed by MAC == 20
    markers only (empty or a handful). The case YAML sets the bin edge per scale instead
    (mid 1000.5, small 200.5 = 0.01 * 2N): a 10% subset has MAF <= 0.1 in bin 1, a 20% subset
    MAF <= 0.05, the full sample MAF <= 0.01."""
    base = ctx.base_traits()
    for p in (ctx.cfg["sgrm"] + ".mtx", ctx.cfg["sgrm"] + ".ids"):
        ctx.inputs[p] = md5(p)
    rng = np.random.default_rng(SEEDS["vrcate"])
    s1 = rng.choice(ctx.N, int(0.10 * ctx.N), replace=False)
    s2 = rng.choice(ctx.N, int(0.20 * ctx.N), replace=False)
    edge = 0.01 * 2 * ctx.N
    traits = [("cF", base["y22"], "base:y22"),
              ("cS1", mask(strong_binary(ctx, rng), s1), "strong binary on 10%"),
              ("cS1b", mask(strong_binary(ctx, rng), s1), "strong binary on the same 10%"),
              ("cS2", mask(strong_binary(ctx, rng), s2), "strong binary on 20%")]
    info = {}
    for name, idx in (("cF", np.arange(ctx.N)), ("cS1", s1), ("cS2", s2)):
        st = ctx.qc(idx)
        pool = st["passVR"]
        info[name] = dict(vr_pool=int(pool.sum()), bin1=int((pool & (st["mac"] <= edge)).sum()),
                          bin2=int((pool & (st["mac"] > edge)).sum()), mac_edge=edge)
    print(json.dumps(info))
    for name, v in info.items():
        require(v["bin1"] >= 5 and v["bin2"] >= 30,
                f"{name}: VR-eligible markers bin1 (MAC<={edge:g}) {v['bin1']}, bin2 {v['bin2']}", problems)
    return write_case(ctx, "vrcate", traits, dict(sparse_grm=ctx.cfg["sgrm"] + ".mtx",
                                                  sparse_grm_ids=ctx.cfg["sgrm"] + ".ids", vr_pool=info,
                                                  predicted=predictions(ctx, traits)))


CASE_FN = {c: globals()[f"case_{c}"] for c in ALL_CASES}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--scale", default="all", choices=["mid", "small", "all"])
    ap.add_argument("--cases", default=",".join(ALL_CASES))
    ap.add_argument("--force", action="store_true", help="regenerate even if the manifest exists")
    a = ap.parse_args()
    scales = ["small", "mid"] if a.scale == "all" else [a.scale]
    cases = [c for c in a.cases.split(",") if c]
    bad = [c for c in cases if c not in CASE_FN]
    if bad:
        raise SystemExit(f"unknown case(s): {bad}; known: {ALL_CASES}")
    failures = {}
    for sc in scales:
        ctx = Ctx(sc, a.workdir)
        for c in cases:
            man = os.path.join(ctx.out, f"{c}.manifest.json")
            if os.path.exists(man) and not a.force:
                m = json.load(open(man))
                if m.get("generator_version") == GENERATOR_VERSION and m.get("coverage_ok"):
                    print(f"[{sc}/{c}] up to date")
                    continue
            print(f"[{sc}/{c}] generating", flush=True)
            problems = []
            m = CASE_FN[c](ctx, problems)
            m["coverage_ok"] = not problems
            m["coverage_problems"] = problems
            with open(man, "w") as f:
                json.dump(m, f, indent=1)
            if problems:
                failures[f"{sc}/{c}"] = problems
    if failures:
        print("COVERAGE NOT MET:", json.dumps(failures, indent=1))
        return 1
    print("all requested cases generated with coverage met")
    return 0


if __name__ == "__main__":
    sys.exit(main())
