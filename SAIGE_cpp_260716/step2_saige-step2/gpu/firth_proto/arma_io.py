"""Minimal reader for armadillo arma_binary files (ARMA_MAT_BIN_FN008, column-major doubles)."""
import numpy as np
def load_arma(path):
    with open(path, 'rb') as f:
        hdr = f.readline().decode().strip()
        dims = f.readline().decode().split()
        r, c = int(dims[0]), int(dims[1])
        data = np.frombuffer(f.read(), dtype=np.float64)
    assert hdr.startswith('ARMA_MAT_BIN'), hdr
    return data[:r * c].reshape((c, r)).T.copy() if c > 1 else data[:r].copy()
