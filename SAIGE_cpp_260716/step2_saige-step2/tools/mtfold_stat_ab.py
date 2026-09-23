#!/usr/bin/env python3
"""fp64 A/B of the three Z_t's on simulated genotypes, using the real models.

Reproduces scoreTestBatchMT's quantitative tail in numpy and reports
max |d(-log10 p)| of the fitted fold and of the naive fold (trait 1's A for
everyone) against the wide path.  3000 simulated markers, log-uniform MAF.
"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtfold_armaio import load_arma

b = '/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out/m/'
traits = ['y%d' % i for i in range(1, 9)]
X0 = load_arma(b + 'y1/X.arma'); N = X0.shape[0]
XtXinv = np.linalg.inv(X0.T @ X0)
rng = np.random.default_rng(7); M = 3000
maf = 10 ** rng.uniform(np.log10(0.002), np.log10(0.5), M)
G = rng.binomial(2, maf[None, :], size=(N, M)).astype(np.float64)
Gsq = (G * G).sum(0)
Z0 = X0.T @ G
A1 = load_arma(b + 'y1/XVX_inv_XV.arma')

def nlp(chi):
    out = np.empty_like(chi)
    for i, c in enumerate(chi):
        z = math.sqrt(c / 2.0); e = math.erfc(z)
        out[i] = -math.log10(e) if e > 1e-300 else (z*z)/math.log(10) + math.log10(z*math.sqrt(math.pi))
    return out

print("%-4s %12s %12s %10s" % ("", "fit", "naive", "max chi2"))
for t in traits:
    A = load_arma(b + t + '/XVX_inv_XV.arma'); XVX = load_arma(b + t + '/XVX.arma')
    S_a = load_arma(b + t + '/S_a.arma').ravel(); res = load_arma(b + t + '/res.arma').ravel()
    tau0 = json.load(open(b + t + '/nullmodel.json'))['theta'][0]
    K = XtXinv @ (X0.T @ A)
    GR = G.T @ res
    def stat(Z):
        zxz = np.einsum('ij,ij->j', Z, XVX @ Z); saz = S_a @ Z
        gwz = np.einsum('ij,ij->j', Z0, Z)
        return (GR - saz) / tau0, zxz * tau0 + Gsq - 2 * gwz
    So, vo = stat(A.T @ G); Sf, vf = stat(K.T @ Z0); Sn, vn = stat(A1.T @ G)
    lo = nlp(So**2/vo); lf = nlp(Sf**2/vf); ln = nlp(Sn**2/vn)
    print("%-4s %12.3e %12.3e %10.1f" % (t, np.max(np.abs(lf-lo)), np.max(np.abs(ln-lo)),
                                          (So**2/vo).max()))
