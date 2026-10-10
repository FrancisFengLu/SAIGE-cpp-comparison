#!/usr/bin/env python3
"""subset_model.py <src model dir> <dst model dir> <seed> <drop fraction>

A null-model directory with a random subset of its samples DROPPED, for timing runs in
which every trait must have its own sample list (missing phenotypes): the per-sample
arrays (mu / res / y / V / offset, X / XXVX_inv rows, XV / XVX_inv_XV columns) lose
those rows and nullmodel.json's n / sampleIDs follow. The p x p constants are copied
as they are, so the model is no longer a consistent fit -- the numbers it produces are
finite but meaningless -- which is fine for measuring the step-2 pipeline (every
code path is the same) and not for anything else."""
import sys, os, json, struct, numpy as np, shutil
src, dst, seed, frac = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
os.makedirs(dst, exist_ok=True)
def rd(p):
    with open(p, "rb") as f:
        h = f.readline().decode().strip(); r, c = map(int, f.readline().decode().split())
        assert h == "ARMA_MAT_BIN_FN008", (p, h)
        a = np.frombuffer(f.read(), dtype="<f8", count=r * c).reshape((c, r)).T   # column-major
    return np.ascontiguousarray(a)
def wr(p, a):
    r, c = a.shape
    with open(p, "wb") as f:
        f.write(f"ARMA_MAT_BIN_FN008\n{r} {c}\n".encode())
        f.write(np.asfortranarray(a.astype("<f8")).tobytes(order="F"))
nm = json.load(open(os.path.join(src, "nullmodel.json")))
n = int(nm["n"]); ids = nm["sampleIDs"]; assert len(ids) == n
rng = np.random.default_rng(seed)
keep = np.sort(rng.choice(n, int(round(n * (1.0 - frac))), replace=False))
for f in os.listdir(src):
    p = os.path.join(src, f)
    if not f.endswith(".arma"):
        if f != "nullmodel.json": shutil.copy(p, os.path.join(dst, f))
        continue
    a = rd(p)
    if a.shape[0] == n: a = a[keep, :]
    elif a.shape[1] == n: a = a[:, keep]
    wr(os.path.join(dst, f), a)
nm["n"] = int(len(keep)); nm["sampleIDs"] = [ids[i] for i in keep]
json.dump(nm, open(os.path.join(dst, "nullmodel.json"), "w"))
print(f"{dst}: n {n} -> {len(keep)}")
