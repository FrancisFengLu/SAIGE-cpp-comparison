#!/usr/bin/env python3
"""Build the (marker, trait) pair list that the production step-2 run sent to Firth,
plus the per-marker code->dosage table (flip / mean imputation / MAC-gated clean),
replicating imputeGenoAndFlip (UTIL.cpp) for hard-call PLINK, alt-first.

Usage: make_pairs.py --bed X.bed --bim X.bim --N N --out DIR --cut 0.05 \
          trait=model_dir=output.txt [trait=model_dir=output.txt ...]
Writes DIR/pairs.tsv (pair, trait_idx, marker_idx, p, beta_out, se_out, flip),
       DIR/traits.tsv (trait_idx, name, model_dir), DIR/lut.npy (M x 4 float64 for ALL markers),
       DIR/markers.npy (int32 marker idx of all markers used).
"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--bed', required=True); ap.add_argument('--bim', required=True)
ap.add_argument('--N', type=int, required=True); ap.add_argument('--out', required=True)
ap.add_argument('--cut', type=float, default=0.05)
ap.add_argument('--routes', default=None, help='route-dump dir: keep rows whose route byte has bit 5 (Firth was run) instead of p <= cut')
ap.add_argument('--zerod_cutoff', type=float, default=0.2); ap.add_argument('--zerod_mac', type=float, default=10.0)
ap.add_argument('specs', nargs='+')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
N = a.N; B = (N + 3) // 4
ids = [l.split()[1] for l in open(a.bim)]
M = len(ids); pos = {k: i for i, k in enumerate(ids)}
bed = np.memmap(a.bed, dtype=np.uint8, mode='r')
assert bed.shape[0] == 3 + M * B, 'bed size'
idx = np.arange(N)
rows = []; traits = []
for ti, spec in enumerate(a.specs):
    name, mdir, outf = spec.split('=')
    traits.append((ti, name, mdir))
    route = None
    if a.routes:
        rr = np.fromfile(os.path.join(a.routes, name + '.route'), dtype=np.dtype([('r', 'u1'), ('g', '<f8')]))
        assert len(rr) == M, (len(rr), M); route = rr['r']
    with open(outf) as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            t = l.rstrip('\n').split('\t')
            p = float(t[c['p.value']])
            keep = (p <= a.cut) if route is None else bool(route[pos[t[c['MarkerID']]]] & 0x20)
            if keep:
                rows.append((ti, pos[t[c['MarkerID']]], p, float(t[c['BETA']]), float(t[c['SE']]), float(t[c['Tstat']]), t[c['Is.SPA']], float(t[c['AF_Allele2']])))
rows.sort()
mk = np.array(sorted(set(r[1] for r in rows)), dtype=np.int32)
lut = np.zeros((M, 4)); flipv = np.zeros(M, dtype=np.int8)
for m in mk:
    col = np.asarray(bed[3 + m * B: 3 + (m + 1) * B])
    code = (col[idx >> 2] >> (2 * (idx & 3))) & 3
    cnt = np.bincount(code, minlength=4)
    nmiss = cnt[1]; nonmiss = N - nmiss
    altc = cnt[2] + 2 * cnt[0]                      # alt-first: A1 (bim col 5) counted, code 00 = hom A1
    af_pre = altc / nonmiss / 2.0 if nonmiss > 0 else 0.0
    flip = af_pre > 0.5
    af_flip = 1 - af_pre if flip else af_pre
    imputeG = 2.0 * af_flip                          # impute_method mean
    maf_pre = min(af_pre, 1 - af_pre)
    mac_clean = maf_pre * N * (1 - nmiss / N) * 2 + imputeG * nmiss
    doclean = a.zerod_cutoff > 0 and mac_clean <= a.zerod_mac
    d0 = np.array([2.0, imputeG, 1.0, 0.0])          # code 0,1,2,3 -> dosage (pre-flip; code 1 missing)
    if flip: d0[[0, 2, 3]] = 2.0 - d0[[0, 2, 3]]
    if doclean: d0[np.abs(d0) < a.zerod_cutoff] = 0.0
    lut[m] = d0; flipv[m] = flip
with open(os.path.join(a.out, 'pairs.tsv'), 'w') as f:
    f.write('pair\ttrait\tmarker\tp\tbeta_out\tse_out\ttstat_out\tisSPA\taf\tflip\n')
    for k, r in enumerate(rows):
        f.write(f'{k}\t{r[0]}\t{r[1]}\t{r[2]:.6E}\t{r[3]!r}\t{r[4]!r}\t{r[5]!r}\t{r[6]}\t{r[7]!r}\t{int(flipv[r[1]])}\n')
with open(os.path.join(a.out, 'traits.tsv'), 'w') as f:
    for t in traits: f.write('%d\t%s\t%s\n' % t)
np.save(os.path.join(a.out, 'lut.npy'), lut); np.save(os.path.join(a.out, 'markers.npy'), mk)
print(f'{len(rows)} pairs, {len(mk)} distinct markers, {len(traits)} traits, flipped markers: {int(flipv[mk].sum())}, N={N}, M={M}')
