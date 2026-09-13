#!/usr/bin/env python3
"""Simulate a rare-variant PLINK1 set on the SAME 50k samples as /opt/saige/data/mid.

Why: `mid` is plink2 --dummy with a ~uniform MAF spectrum, so at maxMAF<=0.01 a
400-marker gene keeps 1-9 markers and NOTHING has MAC<=10.  The ultra-rare
collapse branch and the binary ER branch are therefore unreachable on `mid`.
This writes a WES-like spectrum so both branches are entered.

Output: <prefix>.bed/.bim/.fam  and a group file with 3 annotations.
Variant IDs are spelled "chr:pos:A1:A2" so that they match on BOTH sides:
R merges the group file's SNP against markerInfo$ID (bim column 2) and ID2,
the C++ side matches against chr:pos:ref:alt in either allele order.
"""
import numpy as np
import sys

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/opt/saige/data/rare"
FAM_SRC = "/opt/saige/data/mid.fam"
rng = np.random.default_rng(20260913)

fam = [l.rstrip("\n") for l in open(FAM_SRC) if l.strip()]
N = len(fam)

# (count, MAC low, MAC high) -- MAC of the counted allele (bim A1 == Allele2)
SPEC = [(600, 1, 4),      # ER path (binary, MAC<=4) + ultra-rare
        (900, 5, 10),     # ultra-rare collapse (MAC<=10)
        (600, 11, 20),    # just above the collapse cut
        (600, 21, 100),   # maxMAF 0.001 band (MAC<=100)
        (300, 101, 1000)] # maxMAF 0.01 band
M = sum(s[0] for s in SPEC)

macs = np.concatenate([rng.integers(lo, hi + 1, n) for n, lo, hi in SPEC])
rng.shuffle(macs)

# PLINK1 bed codes: 00 hom-A1, 01 missing, 10 het, 11 hom-A2.
# dosage of A1 (= SAIGE Allele2) is 2 / NA / 1 / 0.
HOM_A2, HET, MISS, HOM_A1 = 0b11, 0b10, 0b01, 0b00
nbytes = (N + 3) // 4
bed = bytearray([0x6c, 0x1b, 0x01])
bim = []
for j in range(M):
    mac = int(macs[j])
    codes = np.full(N, HOM_A2, dtype=np.uint8)
    nhom = mac // 40                       # a few homozygous carriers
    nhet = mac - 2 * nhom
    idx = rng.choice(N, size=nhet + nhom, replace=False)
    codes[idx[:nhet]] = HET
    codes[idx[nhet:]] = HOM_A1
    if j % 2 == 0:                          # half the markers carry missingness
        nmiss = rng.binomial(N, 0.005)
        midx = rng.choice(N, size=nmiss, replace=False)
        codes[midx] = MISS
    packed = np.zeros(nbytes, dtype=np.uint8)
    c = codes.astype(np.uint8)
    for k in range(4):
        sl = c[k::4]
        packed[:len(sl)] |= sl << (2 * k)
    bed += packed.tobytes()
    pos = j + 1
    bim.append("1\t1:%d:A:B\t0\t%d\tA\tB" % (pos, pos))

open(PREFIX + ".bed", "wb").write(bytes(bed))
open(PREFIX + ".bim", "w").write("\n".join(bim) + "\n")
open(PREFIX + ".fam", "w").write("\n".join(fam) + "\n")

# group file: 30 genes x 100 markers, annotations cycling
ANN = ["lof", "missense", "synonymous"]
lines = []
for g in range(M // 100):
    ids = ["1:%d:A:B" % (g * 100 + i + 1) for i in range(100)]
    ann = [ANN[i % 3] for i in range(100)]
    lines.append("RGENE%03d\tvar\t%s" % (g, "\t".join(ids)))
    lines.append("RGENE%03d\tanno\t%s" % (g, "\t".join(ann)))
open(PREFIX + ".group.txt", "w").write("\n".join(lines) + "\n")

print("N=%d M=%d  MAC<=4: %d  MAC<=10: %d  MAC<=20: %d  MAC<=100: %d"
      % (N, M, (macs <= 4).sum(), (macs <= 10).sum(), (macs <= 20).sum(),
         (macs <= 100).sum()))
print("wrote", PREFIX + ".{bed,bim,fam,group.txt}")
