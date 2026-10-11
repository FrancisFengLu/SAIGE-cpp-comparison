#!/usr/bin/env python3
"""gen_covp.py OUTDIR [N] [M] [seed]

covp: the covariate-count gate's data (tests/r_parity/r_runs_covp.sh). N unrelated samples,
M markers (90% common with ALT AF uniform(0.02, 0.98), so about half are flipped; 10% rare with
MAC 1..8), 1% missing genotypes, 40 covariates -- x1 sex (Bernoulli 1/2), x2 age-like, x3..x40
PC-like, all continuous ones N(0, 1), effects shrinking with the index -- and 4 binary traits
(prevalence 5 / 10 / 20 / 40%): b1 and b2 complete (one sample list), b3 with 10% and b4 with
20% missing phenotypes (their own lists). The null models are fitted with the first p - 1
covariates, p = 4, 9, 13, 24, 40 (R step 1, r_runs_covp.sh). Default N = M = 5000; the timing
set covt is the same generator at N = M = 20,000."""
import numpy as np, sys, os
out = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
M = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 20261010
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
Mc = int(M * 0.9); Mr = M - Mc
af = np.concatenate([rng.uniform(0.02, 0.98, Mc), np.zeros(Mr)])
G = np.zeros((N, M), dtype=np.int8)
for j in range(Mc): G[:, j] = rng.binomial(2, af[j], N)
for j in range(Mc, M):
    mac = rng.integers(1, 9); idx = rng.choice(N, mac, replace=False); G[idx, j] = 1
perm = rng.permutation(M); G = G[:, perm]; israre = (perm >= Mc)
miss = rng.random((N, M)) < 0.01
P = 40
X = np.column_stack([(rng.random(N) < 0.5).astype(float), rng.normal(size=(N, P - 1))])
bet = np.concatenate([[0.4, 0.5], 0.3 * np.array([(-1) ** k / np.sqrt(k + 1) for k in range(P - 2)])])
common = np.where(~israre)[0]; rare = np.where(israre)[0]
Gp = G[:, common[:800]].astype(float)
Gs = (Gp - Gp.mean(0)) / np.maximum(Gp.std(0), 1e-9)
poly = Gs @ rng.normal(0, np.sqrt(0.25 / 800), 800)
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
prev = [0.05, 0.10, 0.20, 0.40]
Y = {}
for k in range(4):
    lin = poly + X @ bet + sig * (0.7 + 0.1 * k) + rng.normal(0, 0.3, N)
    Y[f"b{k + 1}"] = binary(prev[k], lin)
MISS = {"b1": np.zeros(N, bool), "b2": np.zeros(N, bool), "b3": rng.random(N) < 0.10, "b4": rng.random(N) < 0.20}
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
    f.write("IID\t" + "\t".join(names) + "\t" + "\t".join(f"x{k}" for k in range(1, P + 1)) + "\n")
    for i in range(N):
        vals = ["NA" if MISS[nm][i] else str(int(Y[nm][i])) for nm in names]
        f.write(f"s{i}\t" + "\t".join(vals) + "\t" + f"{int(X[i, 0])}\t"
                + "\t".join(f"{X[i, k]:.6f}" for k in range(1, P)) + "\n")
mac = np.minimum(G.sum(0), 2 * N - G.sum(0))
print("N", N, "M", M, "ALT AF>0.5:", (G.mean(0) / 2 > 0.5).sum(), "missing cells:", miss.sum(),
      "MAC<=4:", (mac <= 4).sum())
print("prevalence:", {nm: round(float(np.nanmean(np.where(MISS[nm], np.nan, Y[nm]))), 3) for nm in names})
print("NA per trait:", {nm: int(MISS[nm].sum()) for nm in names})
print("causal:", ["snp%d" % c for c in causal])
