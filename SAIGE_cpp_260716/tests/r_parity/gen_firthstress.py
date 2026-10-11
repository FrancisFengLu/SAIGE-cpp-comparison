#!/usr/bin/env python3
"""gen_firthstress.py OUTDIR [N] [seed]

firthstress: the Firth-status gate's data (tests/r_parity/r_runs_firth.sh, cpp_runs_firth.sh,
firth_table.py). N unrelated samples (default 20,000); 2,000 common markers (ALT AF uniform(0.05,
0.95), the GRM) and 5,000 rare ones with MAC uniform in 1..30; 0.5% missing genotypes; 3
covariates (x1 sex, x2 age-like, x3 PC-like); 8 binary traits with few cases (about 100-300 each,
prevalence 0.5-1.5%), b7 with 5% and b8 with 10% missing phenotypes (their own sample lists).
40% of the rare markers are "signal" markers: each picks one trait and draws its carriers from
that trait's cases with probability rho ~ uniform(0.2, 1), so with 1 to 30 carriers and 100-300
cases many Firth fits (p <= 0.05) see quasi-separation (all carriers cases) and run to maxit.
"""
import numpy as np, sys, os
out = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261011
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
Mc, Mr = 2000, 5000
M = Mc + Mr
P = 3
X = np.column_stack([(rng.random(N) < 0.5).astype(float), rng.normal(size=(N, P - 1))])
bet = np.array([0.3, 0.4, 0.25])
# common markers and the polygenic component
G = np.zeros((N, M), dtype=np.int8)
af = rng.uniform(0.05, 0.95, Mc)
for j in range(Mc): G[:, j] = rng.binomial(2, af[j], N)
Gp = G[:, :500].astype(float)
Gs = (Gp - Gp.mean(0)) / np.maximum(Gp.std(0), 1e-9)
poly = Gs @ rng.normal(0, np.sqrt(0.3 / 500), 500)
# traits: shared polygenic + covariates + own noise, intercept bisected to the target prevalence
target = [150, 150, 120, 200, 300, 150, 100, 250]
names = [f"b{k + 1}" for k in range(8)]
Y = {}
for k, nm in enumerate(names):
    lin = poly + X @ bet + rng.normal(0, 0.5, N)
    prev = target[k] / N
    lo, hi = -20, 5
    for _ in range(100):
        mid = (lo + hi) / 2; p = 1 / (1 + np.exp(-(mid + lin)))
        if p.mean() < prev: lo = mid
        else: hi = mid
    Y[nm] = rng.binomial(1, 1 / (1 + np.exp(-(mid + lin))))
MISS = {nm: np.zeros(N, bool) for nm in names}
MISS["b7"] = rng.random(N) < 0.05
MISS["b8"] = rng.random(N) < 0.10
# rare markers: MAC 1..30; 40% signal (carriers enriched in one trait's cases)
cases = {nm: np.where(Y[nm] == 1)[0] for nm in names}
nsig = 0; sigrec = []
for j in range(Mc, M):
    mac = int(rng.integers(1, 31))
    if rng.random() < 0.4:
        nm = names[int(rng.integers(0, 8))]
        rho = rng.uniform(0.2, 1.0)
        pool = cases[nm]
        chosen = set()
        while len(chosen) < mac:
            if rng.random() < rho and len(chosen) < len(pool):
                chosen.add(int(rng.choice(pool)))
            else:
                chosen.add(int(rng.integers(0, N)))
        idx = np.array(sorted(chosen))
        nsig += 1; sigrec.append((j, nm, mac, round(float(rho), 2)))
    else:
        idx = rng.choice(N, mac, replace=False)
    G[idx, j] = 1
perm = rng.permutation(M); G = G[:, perm]
miss = rng.random((N, M)) < 0.005
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
        f.write(f"1\tsnp{perm[j]}\t0\t{1000 + j * 100}\tA\tG\n")
with open(f"{out}/g.fam", "w") as f:
    for i in range(N): f.write(f"s{i}\ts{i}\t0\t0\t1\t-9\n")
with open(f"{out}/pheno.txt", "w") as f:
    f.write("IID\t" + "\t".join(names) + "\t" + "\t".join(f"x{k}" for k in range(1, P + 1)) + "\n")
    for i in range(N):
        vals = ["NA" if MISS[nm][i] else str(int(Y[nm][i])) for nm in names]
        f.write(f"s{i}\t" + "\t".join(vals) + "\t" + f"{int(X[i, 0])}\t"
                + "\t".join(f"{X[i, k]:.6f}" for k in range(1, P)) + "\n")
with open(f"{out}/signal.txt", "w") as f:
    f.write("marker\ttrait\tmac\trho\n")
    for j, nm, mac, rho in sigrec: f.write(f"snp{j}\t{nm}\t{mac}\t{rho}\n")
mac = np.minimum(G.sum(0), 2 * N - G.sum(0))
print("N", N, "M", M, "common", Mc, "rare", Mr, "signal", nsig, "missing cells:", miss.sum(),
      "MAC<=4:", int((mac <= 4).sum()), "MAC<=30:", int((mac <= 30).sum()))
print("cases:", {nm: int(np.sum(Y[nm][~MISS[nm]])) for nm in names})
print("NA per trait:", {nm: int(MISS[nm].sum()) for nm in names})
