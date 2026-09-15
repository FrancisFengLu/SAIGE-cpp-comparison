"""Genotype decode + saige-null's per-marker QC, in numpy.

Mirrors tools/bed_reader/marker_decoder.cpp (decode_marker) as called by
init_global_geno with setminMAC_VarianceRatio(20, -1, true):

  g            = number of A1 alleles (BED 00 -> 2, 10 -> 1, 11 -> 0, 01 -> missing)
  altFreq_pre  = sum(g over non-missing) / (2 * n_nonmiss)            (float32)
  fill         = std::round(2.0f * altFreq_pre)                       (0, 1 or 2)
  alleleCount += fill * numMissing
  altFreq      = alleleCount / (2 * n)                                (float32)
  maf          = min(altFreq, 1 - altFreq);  mac = min(AC, 2n - AC)
  passQC       = maf >= min_maf && missRate <= max_miss
  passVR       = VR enabled && drawn(marker) && mac >= 20; passVR => !passQC

All sample sets are 0-based FAM row indices. The count predicted here is
checked at gate time against the "<k> markers with MAF >= ..." log line of the
solo run, so a divergence between this model and the binary is a case failure,
not a silent mis-design.
"""
import numpy as np

F32 = np.float32


def bed_dims(prefix):
    n = sum(1 for _ in open(prefix + ".fam"))
    m = sum(1 for _ in open(prefix + ".bim"))
    return n, m


def decode(prefix):
    """int8 array (M, N): 0/1/2 = A1 count, 3 = missing."""
    N, M = bed_dims(prefix)
    nb = (N + 3) // 4
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)
    assert raw[0] == 0x6C and raw[1] == 0x1B and raw[2] == 0x01, "not a SNP-major BED"
    raw = raw[3:].reshape(M, nb)
    lut = np.zeros((256, 4), np.int8)
    code = {0: 2, 1: 3, 2: 1, 3: 0}
    for b in range(256):
        for j in range(4):
            lut[b, j] = code[(b >> (2 * j)) & 3]
    G = np.empty((M, N), np.int8)
    step = 2000
    for m0 in range(0, M, step):
        G[m0:m0 + step] = lut[raw[m0:m0 + step]].reshape(-1, nb * 4)[:, :N]
    return G


def counts(G, idx, chunk=2000):
    """Per-marker (allele count over non-missing, missing count) on samples idx."""
    idx = np.sort(np.asarray(idx, dtype=np.int64))
    M = G.shape[0]
    ac = np.zeros(M, np.int64)
    miss = np.zeros(M, np.int64)
    for m0 in range(0, M, chunk):
        sub = G[m0:m0 + chunk][:, idx]
        mm = sub == 3
        miss[m0:m0 + chunk] = mm.sum(1)
        ac[m0:m0 + chunk] = np.where(mm, 0, sub).sum(1, dtype=np.int64)
    return ac, miss, len(idx)


def stats(ac, miss, n, drawn=None, min_maf=0.01, max_miss=0.15, vr_min_mac=20):
    nn = n - miss
    with np.errstate(divide="ignore", invalid="ignore"):
        altpre = np.where(nn > 0, ac.astype(F32) / (nn * 2).astype(F32), F32(0)).astype(F32)
    # std::round on the float32 product, half away from zero (values are >= 0)
    fill = np.floor((F32(2.0) * altpre).astype(np.float64) + 0.5).astype(np.int64)
    ac2 = ac + fill * miss
    alt = (ac2.astype(F32) / F32(n * 2)).astype(F32)
    maf = np.minimum(alt, F32(1) - alt)
    mac = np.minimum(ac2, 2 * n - ac2)
    mr = (miss.astype(F32) / F32(n)).astype(F32)
    passQC = (maf >= F32(min_maf)) & (mr <= F32(max_miss))
    if drawn is not None:
        passVR = drawn & (mac >= vr_min_mac)
        passQC &= ~passVR
    else:
        passVR = np.zeros_like(passQC)
    return dict(passQC=passQC, passVR=passVR, maf=maf, alt=alt, altpre=altpre, fill=fill,
                mr=mr, mac=mac, ac=ac, miss=miss, n=n)


def qc(G, idx, drawn=None, **kw):
    ac, miss, n = counts(G, idx)
    return stats(ac, miss, n, drawn, **kw)


def drawn_mask(M, vrdraw_bin):
    import subprocess
    out = subprocess.run([vrdraw_bin, str(M)], check=True, capture_output=True, text=True).stdout
    d = np.zeros(M, bool)
    d[np.array([int(x) for x in out.split()], dtype=np.int64)] = True
    return d
