#!/usr/bin/env python3
"""Synthetic edge-case markers on the real sample set (N=50000): ultra-low MAC
(1..6), carriers all cases / all controls / mixed (complete and quasi-complete
separation), a monomorphic marker (singular Fisher information), markers with
missing calls (mean imputation, MAC-gated zeroing), one common marker and one
flipped (AF>0.5) marker. Writes a .bed/.bim/.fam and a pair directory in the
same layout as make_pairs.py (no production output -> beta_out = nan).
Traits: c01_1 (1% cases) and y1 (52% cases). Usage: gen_edge.py OUTDIR"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from arma_io import load_arma
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
N = 50000; rng = np.random.default_rng(11)
traits = [('c01_1', '/opt/saige/logs/spagaps/step1/out/m/c01_1'), ('y1', '/opt/saige/logs/binsplit/models/spa/m/y1')]
Y = {t: load_arma(d + '/y.arma') for t, d in traits}
cases = {t: np.nonzero(Y[t] == 1)[0] for t, _ in traits}; ctrls = {t: np.nonzero(Y[t] == 0)[0] for t, _ in traits}
cols = []; names = []
def add(name, G): cols.append(G.astype(np.int8)); names.append(name)
for t, _ in traits:
    for mac in [1, 2, 3, 4, 5, 6]:
        for kind in ['allcase', 'allctrl', 'mixed']:
            G = np.zeros(N, np.int8)
            if kind == 'allcase': idx = rng.choice(cases[t], mac, replace=False)
            elif kind == 'allctrl': idx = rng.choice(ctrls[t], mac, replace=False)
            else: idx = np.concatenate([rng.choice(cases[t], (mac + 1) // 2, replace=False), rng.choice(ctrls[t], mac // 2, replace=False)])
            G[idx] = 1; add(f'{t}_mac{mac}_{kind}', G)
    G = np.zeros(N, np.int8); idx = rng.choice(cases[t], 3, replace=False); G[idx] = 2; add(f'{t}_hom2x3_allcase', G)
add('mono', np.zeros(N, np.int8))
for mac in [3, 8]:
    G = np.zeros(N, np.int8); G[rng.choice(N, mac, replace=False)] = 1; G[rng.choice(N, 500, replace=False)] = -1; add(f'miss500_mac{mac}', G)
G = rng.binomial(2, 0.3, N).astype(np.int8); add('common_af0.3', G)
G = rng.binomial(2, 0.8, N).astype(np.int8); G[rng.choice(N, 300, replace=False)] = -1; add('flip_af0.8_miss300', G)
G = rng.binomial(2, 0.001, N).astype(np.int8); add('rare_af0.001', G)
M = len(cols); B = (N + 3) // 4
code = {0: 3, 1: 2, 2: 0, -1: 1}
bed = bytearray([0x6c, 0x1b, 0x01])
for G in cols:
    c = np.array([code[int(v)] for v in G], np.uint8); w = np.zeros(B, np.uint8)
    for j in range(4): w[: (N - j + 3) // 4] |= (c[j::4] << (2 * j)).astype(np.uint8)
    bed += bytes(w)
open(out + '/edge.bed', 'wb').write(bed)
with open(out + '/edge.bim', 'w') as f:
    for i, n in enumerate(names): f.write(f'1\t{n}\t0\t{i+1}\tA\tB\n')
with open(out + '/edge.fam', 'w') as f:
    for i in range(N): f.write(f'per{i} per{i} 0 0 0 -9\n')
# pair dir: every marker x every trait
pairs = [(ti, m) for ti in range(len(traits)) for m in range(M)]
import subprocess
with open(out + '/pairs.tsv', 'w') as f:
    f.write('pair\ttrait\tmarker\tp\tbeta_out\tse_out\ttstat_out\tisSPA\taf\tflip\n')
    for k, (ti, m) in enumerate(pairs): f.write(f'{k}\t{ti}\t{m}\tnan\tnan\tnan\tnan\tNA\tnan\t0\n')
with open(out + '/traits.tsv', 'w') as f:
    for ti, (t, d) in enumerate(traits): f.write(f'{ti}\t{t}\t{d}\n')
# LUT via the same rule as make_pairs.py (reuse its code path by importing is awkward; replicate)
lut = np.zeros((M, 4)); idx = np.arange(N); bedm = np.frombuffer(bytes(bed), np.uint8)
for m in range(M):
    col = bedm[3 + m * B: 3 + (m + 1) * B]; cd = (col[idx >> 2] >> (2 * (idx & 3))) & 3; cnt = np.bincount(cd, minlength=4)
    nmiss = cnt[1]; nonmiss = N - nmiss; altc = cnt[2] + 2 * cnt[0]; af_pre = altc / nonmiss / 2.0 if nonmiss else 0.0
    flip = af_pre > 0.5; af_flip = 1 - af_pre if flip else af_pre; imputeG = 2 * af_flip
    maf_pre = min(af_pre, 1 - af_pre); mac_clean = maf_pre * N * (1 - nmiss / N) * 2 + imputeG * nmiss
    doclean = mac_clean <= 10.0
    d0 = np.array([2.0, imputeG, 1.0, 0.0])
    if flip: d0[[0, 2, 3]] = 2.0 - d0[[0, 2, 3]]
    if doclean: d0[np.abs(d0) < 0.2] = 0.0
    lut[m] = d0
np.save(out + '/lut.npy', lut); np.save(out + '/markers.npy', np.arange(M, dtype=np.int32))
print(f'{M} markers x {len(traits)} traits = {len(pairs)} pairs ->', out); print('\n'.join(f'{i} {n}' for i, n in enumerate(names)))
