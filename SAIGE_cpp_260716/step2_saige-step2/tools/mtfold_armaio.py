import numpy as np

def load_arma(path):
    with open(path,'rb') as f:
        magic = f.readline().decode().strip()
        dims  = f.readline().decode().strip().split()
        r, c = int(dims[0]), int(dims[1])
        buf = f.read()
    a = np.frombuffer(buf, dtype='<f8', count=r*c)
    if c == 1 and 'COL' in magic:
        return a.copy()
    return a.reshape((c, r)).T.copy()   # column-major
