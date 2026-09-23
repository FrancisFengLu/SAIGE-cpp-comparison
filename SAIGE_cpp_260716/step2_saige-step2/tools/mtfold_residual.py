#!/usr/bin/env python3
"""Residual of the three ways to get the p x p fold map K_t, on the real models.

A_t = XVX_inv_XV of trait t as step 1 stored it (n x p, = V_t X (X'V_t X)^-1).
  fit             K_t = (X0'X0)^-1 X0' A_t          (what the C++ does)
  analytic tau0   K_t = (1/tau0_json) * XVX_inv_t   (tau0 from nullmodel.json)
  analytic V0     K_t = V[0] * XVX_inv_t            (V as step 1 stored it)
  naive           A_t replaced by A_1 wholesale     (the MARG_SHARECOV experiment)
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtfold_armaio import load_arma

b = '/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out/m/'
traits = ['y%d' % i for i in range(1, 9)]
X0 = load_arma(b + 'y1/X.arma')
print("X byte-identical across traits:",
      all(np.array_equal(X0, load_arma(b + t + '/X.arma')) for t in traits))
XtXinv = np.linalg.inv(X0.T @ X0)
A1 = load_arma(b + 'y1/XVX_inv_XV.arma')
print("%-4s %11s %11s %11s %11s" % ("", "fit", "analytic(tau0)", "analytic(V0)", "naive"))
for t in traits:
    A = load_arma(b + t + '/XVX_inv_XV.arma')
    Xi = load_arma(b + t + '/XVX_inv.arma')
    V0 = load_arma(b + t + '/V.arma').ravel()[0]
    tau0 = json.load(open(b + t + '/nullmodel.json'))['theta'][0]
    s = np.max(np.abs(A))
    print("%-4s %11.3e %11.3e %11.3e %11.3e" % (
        t,
        np.max(np.abs(A - X0 @ (XtXinv @ (X0.T @ A)))) / s,
        np.max(np.abs(A - (1.0 / tau0) * (X0 @ Xi))) / s,
        np.max(np.abs(A - V0 * (X0 @ Xi))) / s,
        np.max(np.abs(A - A1)) / s))
