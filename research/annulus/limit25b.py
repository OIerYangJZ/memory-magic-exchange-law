"""At the limit a0=a=b=T=8/5: for each g and y up to eta0(g)+xi, the margins of the LOW routes
that are NOT the E1<=T identity.  For rank d in {2,3}: m_d = E1 - N_d(X*) where X* = end of the flat part
of route P (route P = E1 exactly for X <= X*); if m_d >= 0 the max-min equals min(N_d(0),E1) <= E1.
Also m01 = E1 - v01, m4 = E1 - v4."""
import sys, numpy as np
a0 = a = b = T = 1.6
sys.argv = ['m', repr(a0), repr(a), repr(b), repr(T)]
import model2 as M
E1 = M.E1
def eta0(g):
    S0, U0 = M.box_S0(g), M.box_U0(g)
    if M.cost(S0, U0, g, 0.0) <= 1e-12: return 0.0
    lo, hi = 0.0, 3.0
    for _ in range(40):
        mid = (lo+hi)/2
        if M.cost(S0, U0, g, mid) <= 1e-12: hi = mid
        else: lo = mid
    return hi
worst = {}
for xi in [0.0, 0.02, 0.05]:
    mins = dict(m01=9, m2=9, m3=9, m4=9)
    arg = {}
    for g in np.arange(0.405, 2.0, 0.01):
        e0 = eta0(g)
        gam = min(g, a0); S0, U0 = M.box_S0(g), M.box_U0(g); Sp = max(S0, U0)
        for y in np.linspace(0.0, e0 + xi, 60)[1:]:
            pm = M.pieces_max(g, y, y)
            m01 = E1 - pm
            m4 = E1 - (max(0, 4*y - 2*gam) + pm)
            # rank 3: route P flat for X <= X3 = y - (Sp+U0+2); N3(X) = (3y-gam-X)_+ + pm
            X3 = max(0.0, y - (Sp + U0 + 2)); N3 = max(0, 3*y - gam - X3) + pm; m3 = E1 - N3
            # rank 2: route P flat for X <= X2 = y - (Sp+1); N2(X) = (2y - max(X, G3p))_+ + pm
            G3p = T - 4 + a0 + gam
            X2 = max(0.0, y - (Sp + 1)); N2 = max(0, 2*y - max(X2, G3p)) + pm; m2 = E1 - N2
            for k, v in (('m01', m01), ('m2', m2), ('m3', m3), ('m4', m4)):
                if v < mins[k]: mins[k] = v; arg[k] = (round(g,3), round(y,3), round(e0,3))
    print(f"xi={xi}: " + "  ".join(f"{k}={mins[k]:+.4f}@{arg[k]}" for k in mins), flush=True)
