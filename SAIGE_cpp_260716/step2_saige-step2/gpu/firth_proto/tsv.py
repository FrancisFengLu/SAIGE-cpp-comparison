"""Tiny TSV -> dict-of-numpy-columns reader (no pandas on this box)."""
import numpy as np
def read_tsv(path):
    with open(path) as f:
        h = f.readline().rstrip('\n').split('\t'); rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    out = {}
    for j, k in enumerate(h):
        col = [r[j] for r in rows]
        try: out[k] = np.array([float(v[11:-1] if v.startswith('np.float64(') else v) for v in col])
        except ValueError: out[k] = np.array(col)
    return out
def join(a, b, keys=('pair',)):
    ka = {tuple(int(a[k][i]) for k in keys): i for i in range(len(a[keys[0]]))}
    ia, ib = [], []
    for i in range(len(b[keys[0]])):
        t = tuple(int(b[k][i]) for k in keys)
        if t in ka: ia.append(ka[t]); ib.append(i)
    ia, ib = np.array(ia), np.array(ib)
    return {k: v[ia] for k, v in a.items()} | {k: v[ib] for k, v in b.items() if k not in a}
