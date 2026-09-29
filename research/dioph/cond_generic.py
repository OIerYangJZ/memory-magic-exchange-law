"""Conditional variant: frames for which the low-height normal count is at most the volume term
max(0, 4eta-2gam) (e.g. badly approximable Pi at all relevant scales).  Float only."""
import numpy as np, fast
def normals_generic(a0, gl, eta, Fstar, dichotomy=True):
    gam = np.minimum(gl, a0)
    return np.maximum(0.0, 4 * eta - 2 * gam)
fast.normals = normals_generic
for al in [2.44, 2.46, 2.48, 2.49]:
    print(al, fast.best_ab(al, a_range=(0, 0.02), b_range=(1.55, 1.61), step=0.01, den=50, deta=0.01), flush=True)
