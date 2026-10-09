#!/usr/bin/env python3
# qt12: N=4000, M=5000 (4500 common with ALT AF uniform(0.02,0.98) -> about half flipped,
# 500 rare MAC 1..8), 1% missing genotypes, 12 covariates (8 continuous, 4 binary),
# 2 quantitative + 2 binary traits, two of them with missing phenotypes, planted signals.
import numpy as np, sys
out = sys.argv[1] if len(sys.argv) > 1 else "/opt/saige/logs/rdefaults/data/qt12"
rng = np.random.default_rng(20261009)
N = 4000; Mc = 4500; Mr = 500; M = Mc + Mr
af = np.concatenate([rng.uniform(0.02, 0.98, Mc), np.zeros(Mr)])
G = np.zeros((N, M), dtype=np.int8)
for j in range(Mc): G[:, j] = rng.binomial(2, af[j], N)
for j in range(Mc, M):
    mac = rng.integers(1, 9); idx = rng.choice(N, mac, replace=False); G[idx, j] = 1
perm = rng.permutation(M); G = G[:, perm]; israre = (perm >= Mc)
miss = rng.random((N, M)) < 0.01
X = np.column_stack([rng.normal(size=(N, 8)), (rng.random((N, 4)) < [0.3, 0.4, 0.5, 0.6]).astype(float)])
bet = np.array([0.3, -0.2, 0.15, 0.1, -0.1, 0.05, 0.2, -0.15, 0.4, -0.3, 0.2, 0.1])
Gs = (G - G.mean(0)) / np.maximum(G.std(0), 1e-9)
common = np.where(~israre)[0]; rare = np.where(israre)[0]
poly = Gs[:, common[:800]] @ rng.normal(0, np.sqrt(0.25 / 800), 800)
causal = list(common[[5, 150, 1200, 2600, 3900]]) + list(rare[[2, 40]])
beta = np.array([0.30, -0.25, 0.20, 0.35, 0.18, 2.0, 1.6])
lin = poly + X @ bet + (G[:, causal] - G[:, causal].mean(0)) @ beta
q1 = lin + rng.normal(0, 1.0, N)
q2 = 0.7 * lin + rng.normal(0, 1.3, N)
def binary(prev, l):
    lo, hi = -15, 15
    for _ in range(80):
        mid = (lo + hi) / 2; p = 1 / (1 + np.exp(-(mid + l)))
        if p.mean() < prev: lo = mid
        else: hi = mid
    return rng.binomial(1, 1 / (1 + np.exp(-(mid + l))))
b1 = binary(0.10, lin); b2 = binary(0.35, lin)
# missing phenotypes: q2 15%, b1 10%
mq2 = rng.random(N) < 0.15; mb1 = rng.random(N) < 0.10
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
with open(f"{out}/pheno.txt", "w") as f:
    f.write("IID\tq1\tq2\tb1\tb2\t" + "\t".join(f"x{k}" for k in range(1, 13)) + "\n")
    for i in range(N):
        f.write(f"s{i}\t{q1[i]:.6f}\t{'NA' if mq2[i] else f'{q2[i]:.6f}'}\t{'NA' if mb1[i] else b1[i]}\t{b2[i]}\t"
                + "\t".join(f"{X[i, k]:.6f}" if k < 8 else f"{int(X[i, k])}" for k in range(12)) + "\n")
mac = np.minimum(G.sum(0), 2 * N - G.sum(0))
print("N", N, "M", M, "ALT AF>0.5:", (G.mean(0) / 2 > 0.5).sum(), "missing cells:", miss.sum(),
      "MAC<=4:", (mac <= 4).sum(), "prev b1/b2:", b1.mean(), b2.mean(), "NA q2/b1:", mq2.sum(), mb1.sum())
print("causal:", ["snp%d" % c for c in causal])
