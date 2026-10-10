#!/usr/bin/env python3
# fam32: N=4000 in 1000 families of 4 (two parents, two children drawn by Mendelian
# transmission), M=5000 (4500 common with ALT AF uniform(0.02,0.98) -> about half flipped,
# 500 rare with MAC 1..8 among the founders), 1% missing genotypes, 12 covariates
# (8 continuous, 4 binary), 16 binary traits (prevalence 2-50%) and 16 quantitative
# traits, EVERY trait with its own missing-phenotype pattern (5-25% NA), a shared
# polygenic component and a family random effect (so a sparse GRM has 4 x 4 blocks).
import numpy as np, sys
out = sys.argv[1] if len(sys.argv) > 1 else "/opt/saige/logs/gpuprep/data/fam32"
rng = np.random.default_rng(20261010)
F = 1000; N = 4 * F; Mc = 4500; Mr = 500; M = Mc + Mr
af = np.concatenate([rng.uniform(0.02, 0.98, Mc), np.zeros(Mr)])
# haplotypes of the 2000 founders: H[founder, hap, marker]
H = np.zeros((2 * F, 2, M), dtype=np.int8)
for j in range(Mc): H[:, :, j] = rng.random((2 * F, 2)) < af[j]
for j in range(Mc, M):
    mac = rng.integers(1, 9)
    idx = rng.choice(2 * F * 2, mac, replace=False)
    H[idx // 2, idx % 2, j] = 1
G = np.zeros((N, M), dtype=np.int8)
for f in range(F):
    pa, ma = H[2 * f], H[2 * f + 1]
    G[4 * f] = pa.sum(0); G[4 * f + 1] = ma.sum(0)
    for c in range(2):
        tp = rng.integers(0, 2, M); tm = rng.integers(0, 2, M)
        G[4 * f + 2 + c] = pa[tp, np.arange(M)] + ma[tm, np.arange(M)]
perm = rng.permutation(M); G = G[:, perm]; israre = (perm >= Mc)
miss = rng.random((N, M)) < 0.01
X = np.column_stack([rng.normal(size=(N, 8)), (rng.random((N, 4)) < [0.3, 0.4, 0.5, 0.6]).astype(float)])
bet = np.array([0.3, -0.2, 0.15, 0.1, -0.1, 0.05, 0.2, -0.15, 0.4, -0.3, 0.2, 0.1])
Gs = (G - G.mean(0)) / np.maximum(G.std(0), 1e-9)
common = np.where(~israre)[0]; rare = np.where(israre)[0]
poly = Gs[:, common[:800]] @ rng.normal(0, np.sqrt(0.25 / 800), 800)
famid = np.repeat(np.arange(F), 4)
causal = list(common[[5, 150, 1200, 2600, 3900]]) + list(rare[[2, 40]])
beta = np.array([0.30, -0.25, 0.20, 0.35, 0.18, 2.0, 1.6])
sig = (G[:, causal] - G[:, causal].mean(0)) @ beta
def binary(prev, l):
    lo, hi = -15, 15
    for _ in range(80):
        mid = (lo + hi) / 2; p = 1 / (1 + np.exp(-(mid + l)))
        if p.mean() < prev: lo = mid
        else: hi = mid
    return rng.binomial(1, 1 / (1 + np.exp(-(mid + l))))
prev = [0.02, 0.05, 0.05, 0.1, 0.1, 0.15, 0.2, 0.2, 0.25, 0.3, 0.3, 0.35, 0.4, 0.45, 0.5, 0.5]
Y = {}
for k in range(16):
    u = rng.normal(0, 0.6, F)[famid]
    lin = 0.8 * poly + X @ bet + sig * (0.6 + 0.05 * k) + u
    Y[f"b{k + 1}"] = binary(prev[k], lin)
for k in range(16):
    u = rng.normal(0, 0.6, F)[famid]
    lin = poly + X @ bet + sig * (0.7 + 0.04 * k) + u
    e = rng.normal(0, 1.0, N) if k % 4 != 3 else (rng.chisquare(3, N) - 3) / np.sqrt(6)
    Y[f"q{k + 1}"] = 0.9 * lin + e
# every trait its own missing pattern
MISS = {}
for i, name in enumerate(list(Y)):
    frac = 0.05 + 0.20 * (i / 31)
    MISS[name] = rng.random(N) < frac
Gm = G.copy(); Gm[miss] = -1
lut = np.array([0b11, 0b10, 0b00, 0b01])  # dosage 0,1,2,missing -> PLINK code
nb = (N + 3) // 4
with open(f"{out}/g.bed", "wb") as f:
    f.write(bytes([0x6c, 0x1b, 0x01]))
    for j in range(M):
        g = Gm[:, j].astype(int); g[g < 0] = 3
        c = lut[g]; pad = np.zeros(nb * 4, dtype=np.uint8); pad[:N] = c
        pad = pad.reshape(nb, 4); byte = pad[:, 0] | (pad[:, 1] << 2) | (pad[:, 2] << 4) | (pad[:, 3] << 6)
        f.write(byte.astype(np.uint8).tobytes())
with open(f"{out}/g.bim", "w") as f:
    for j in range(M):
        f.write(f"{1 if j < M // 2 else 2}\tsnp{j}\t0\t{1000 + j * 100}\tA\tG\n")
with open(f"{out}/g.fam", "w") as f:
    for i in range(N): f.write(f"s{i}\ts{i}\t0\t0\t1\t-9\n")
names = list(Y)
with open(f"{out}/pheno.txt", "w") as f:
    f.write("IID\t" + "\t".join(names) + "\t" + "\t".join(f"x{k}" for k in range(1, 13)) + "\n")
    for i in range(N):
        vals = []
        for nm in names:
            if MISS[nm][i]: vals.append("NA")
            elif nm.startswith("b"): vals.append(str(int(Y[nm][i])))
            else: vals.append(f"{Y[nm][i]:.6f}")
        f.write(f"s{i}\t" + "\t".join(vals) + "\t"
                + "\t".join(f"{X[i, k]:.6f}" if k < 8 else f"{int(X[i, k])}" for k in range(12)) + "\n")
mac = np.minimum(G.sum(0), 2 * N - G.sum(0))
print("N", N, "M", M, "ALT AF>0.5:", (G.mean(0) / 2 > 0.5).sum(), "missing cells:", miss.sum(),
      "MAC<=4:", (mac <= 4).sum())
print("prevalence:", {nm: round(float(np.nanmean(np.where(MISS[nm], np.nan, Y[nm]))), 3) for nm in names if nm.startswith("b")})
print("NA per trait:", {nm: int(MISS[nm].sum()) for nm in names})
print("distinct missing patterns:", len({MISS[nm].tobytes() for nm in names}))
