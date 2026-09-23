#!/usr/bin/env python3
"""Minimal reader for a step-2 .sgs TRAIT file (see sgs_format.hpp).

Returns full-precision fp64 BETA / SE / Tstat / var plus the p-value string,
which is what the text output quantises to 6 printed digits.
"""
import struct, numpy as np

C_BETA,C_SE,C_TSTAT,C_VAR,C_PVAL,C_PVALNA,C_ISSPA = 1,2,3,4,5,6,7
C_N = 22
NAME = {1:"BETA",2:"SE",3:"Tstat",4:"var",5:"p.value",6:"p.value.NA",7:"Is.SPA",
        8:"BETA_c",9:"SE_c",10:"Tstat_c",11:"var_c",12:"pval_c",13:"pvalNA_c",
        14:"AF_case",15:"AF_ctrl",16:"N_case",17:"N_ctrl",
        18:"N_case_hom",19:"N_case_het",20:"N_ctrl_hom",21:"N_ctrl_het",22:"N"}
POD4 = {C_N, 16, 17}     # uint32 columns
POD1 = {C_ISSPA}         # char columns
F64  = {1,2,3,4,8,9,10,11,14,15,18,19,20,21}
PVAL = {5,6,12,13}

class R:
    def __init__(s,b): s.b=b; s.i=0
    def take(s,n): q=s.b[s.i:s.i+n]; s.i+=n; return q
    def u8(s): return s.take(1)[0]
    def u32(s): return struct.unpack('<I', s.take(4))[0]
    def u64(s): return struct.unpack('<Q', s.take(8))[0]
    def st(s): n=s.u32(); return s.take(n).decode()
    def sstr(s):
        n=s.u8()
        if n==255: n=s.u32()
        return s.take(n).decode()

def read_trait(path):
    b = open(path,'rb').read()
    r = R(b)
    assert r.take(8) == b'SAIGESGT', path
    ver = r.u32(); flags = r.u32()
    f32 = bool(flags & 2)
    hdrline = r.st(); name = r.st(); ttype = r.st()
    isCond = r.u8(); isMore = r.u8()
    markerPath = r.st(); textPath = r.st()
    ncol = r.u32(); cols = [r.u8() for _ in range(ncol)]
    out = {NAME[c]: [] for c in cols}
    nrows_tot = 0
    while True:
        magic = struct.unpack('<I', b[r.i:r.i+4])[0]
        if magic != 0x314B4C42: break
        r.u32(); n = r.u32(); bflags = r.u32()
        nrows_tot += n
        if bflags & 1: r.take(n)
        for bit, _ in ((2,'AC'),(4,'AF'),(8,'MISS')):
            if bflags & bit: rd_f64(r, n, f32)
        for c in cols:
            if c in F64:      out[NAME[c]].append(rd_f64(r, n, f32))
            elif c in PVAL:   out[NAME[c]].append(rd_pval(r, n, f32))
            elif c in POD4:   out[NAME[c]].append(rd_pod(r, n, 4, '<u4'))
            elif c in POD1:   out[NAME[c]].append(rd_pod(r, n, 1, 'u1'))
            else: raise RuntimeError("col %d" % c)
    assert struct.unpack('<I', b[r.i:r.i+4])[0] == 0x21444E45, "no END!"
    return {k: np.concatenate(v) for k, v in out.items()}, nrows_tot

def rd_f64(r, n, f32):
    e = r.u8()
    if e == 0: return np.frombuffer(r.take(8*n), '<f8').copy()
    if e == 1: return np.full(n, struct.unpack('<d', r.take(8))[0])
    if e == 3: return np.frombuffer(r.take(4*n), '<f4').astype('f8')
    if e == 4: return np.full(n, struct.unpack('<f', r.take(4))[0], dtype='f8')
    raise RuntimeError("enc %d" % e)

def rd_pod(r, n, w, dt):
    e = r.u8()
    if e == 0: return np.frombuffer(r.take(w*n), dt).copy()
    if e == 1: return np.full(n, np.frombuffer(r.take(w), dt)[0])
    raise RuntimeError("enc %d" % e)

def rd_pval(r, n, f32):
    e = r.u8()
    if e == 2:   d = np.frombuffer(r.take(8*n), '<f8').copy()
    elif e == 5: d = np.frombuffer(r.take(4*n), '<f4').astype('f8')
    else: raise RuntimeError("pval enc %d" % e)
    nexc = r.u32()
    out = d.astype(object)
    for _ in range(nexc):
        i = r.u32(); s = r.sstr(); out[i] = s
    return out
