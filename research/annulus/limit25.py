"""Continuum structure of the lemma set at the limit a0 = a = b = T = 8/5 (E1 = 8/5).
For each g: eta0(g) = least eta with HIGH cost = 0, and the max over y <= eta0 of the
continuum LOW bound (layers of width -> 0: offsets and per-sphere at the same y).  Prints the margin."""
import sys, numpy as np
a0 = a = b = T = 1.6
sys.argv = ['m', repr(a0), repr(a), repr(b), repr(T)]
import model2 as M
M.Xs = np.linspace(0, 3, 3001)
res = []
for g in np.arange(0.005, 2.02, 0.01):
    S0, U0 = M.box_S0(g), M.box_U0(g)
    # eta0: bisection on cost(g, eta) == 0
    lo, hi = 0.0, 3.0
    if M.cost(S0, U0, g, 0.0) <= 1e-12: eta0 = 0.0
    else:
        for _ in range(40):
            mid = (lo + hi) / 2
            if M.cost(S0, U0, g, mid) <= 1e-12: hi = mid
            else: lo = mid
        eta0 = hi
    ys = np.linspace(0, eta0, 121)[1:] if eta0 > 0 else []
    lowmax = max([M.low_layer(g, y, y) for y in ys], default=-np.inf)
    res.append((g, eta0, lowmax, M.E1 - lowmax))
    print(f"g={g:.3f} eta0={eta0:.4f} maxLOW={lowmax:.5f} margin={M.E1-lowmax:+.5f}", flush=True)
m = min(res, key=lambda r: r[3]); print("MIN MARGIN", m)
