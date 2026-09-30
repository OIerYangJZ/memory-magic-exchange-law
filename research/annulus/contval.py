"""Continuum value(g) (layers -> 0: eta layers of width 1e-3, offsets/per-sphere at the same y) at
a0 = a, b = (8-2a0)/3 (the kappa=0 limit), T = a0.  Reports max_g value - E1."""
import sys, numpy as np
a0 = float(sys.argv[1]); a = a0; b = (8 - 2*a0)/3; T = a0
sys.argv = ['m', repr(a0), repr(a), repr(b), repr(T)]
import model2 as M
M.Xs = np.linspace(0, 3, 3001)
def value(g, deta=0.001):
    S0, U0 = M.box_S0(g), M.box_U0(g)
    best = M.E1 + M.cost(S0, U0, g, 0.0); run = -np.inf
    for eta in np.arange(deta, 2.0, deta):
        run = max(run, M.low_layer(g, eta, eta))        # continuum: same y
        if run >= best: break
        v = max(M.E1 + M.cost(S0, U0, g, eta), run)
        best = min(best, v)
    return best
worst = max((value(g) - M.E1, g) for g in np.arange(0.305, 2.0, 0.02))
print(f"a0={a0:.4f} E1={M.E1:.5f} slack={a0-M.E1:.5f}  max_g(value-E1)={worst[0]:+.6f} at g={worst[1]:.3f}", flush=True)
