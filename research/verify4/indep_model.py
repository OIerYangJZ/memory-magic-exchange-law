"""Independent float re-implementation of the D1-D7 exponent model (written from DIOPH_NOTES.md,
not from cert_dioph.py).  Pointwise in g (no interval endpoints), continuous eta layers and
continuous crossing depth D, brute-force cell optimisation on a grid.  Gives the 'true' optimum
of the lemma set, to compare with the certificate."""
import numpy as np, sys
a0 = float(eval(sys.argv[1])) if len(sys.argv) > 1 else 28/17
a  = float(eval(sys.argv[2])) if len(sys.argv) > 2 else 28/17
b  = float(eval(sys.argv[3])) if len(sys.argv) > 3 else 1.57
T  = float(eval(sys.argv[4])) if len(sys.argv) > 4 else 11183/6800
kappa = max(0.0, 2 - a/2 - 3*b/4); E1 = 2*(a-a0) + b + kappa
XS_grid = np.linspace(-2.6, 1.0, 1801)
def cellcost(S0, U0, g, eta):
    """min over admissible (XL, XS) of N + max(0,P); XL on a grid too"""
    XS = XS_grid[:, None]
    XL = np.linspace(-2.6, 1.0, 1801)[None, :]
    ok = (XS <= XL) & (XL < 1 - g - 1e-9) & (2*XL - 1 <= XS)
    N = np.maximum.reduce([np.zeros_like(XS*XL), S0 - XL + 0*XS, U0 - XS + 0*XL, S0 + U0 - XL - XS, 2*(U0 - XS) + 0*XL])
    P = (XL + XS + np.minimum(XS, 2*XL - (1 - g)) - eta)/3 + 1
    c = np.where(ok, N + np.maximum(0, P), np.inf)
    return c.min()
def cost(S0, U0, g, eta):
    c = cellcost(S0, U0, g, eta)
    if 2*b >= a: c = min(c, max(0.0, (S0 + 2*U0 - eta)/3 + 1))
    return c
def box_S0(g): return min(1 - b, (2 - a - min(g, a))/2, 1 - g)
def box_U0(g): return min(1 - a, 1 - g)
def normals(g, eta):
    gam = min(g, a0); G = max(0.0, gam + T/2 - 2)
    return max(eta, 2*eta - G, 3*eta - gam, 4*eta - 2*gam)
def low_at(g, eta):
    """normals(eta) + max over classes of offsets(eta) + per-sphere(eta)  (continuous layers)"""
    lde = max(1 - g - a0, 1 - 2*a0); ldb = max(1 - g - a, 1 - 2*a)
    UT = min(1 - a0, 1 - g); U0 = box_U0(g)
    best = -np.inf
    Ds = [None] + list(np.linspace(lde, max(lde, 1 - 2*g), 25))
    for D in Ds:
        if D is None:
            offs = max(0, eta + lde + 1); ST = min((1 + lde)/2, 1 - g); S0 = min(1 - b, (1 + ldb)/2, 1 - g)
        else:
            offs = max(0, eta + D + 1)
            ST = min((1 + lde)/2, lde + (1 - D)/2, 1 - g)
            S0 = min(1 - b, (1 + ldb)/2, ldb + (1 - D)/2, 1 - g)
        tube = max(0, (ST + 2*UT - eta)/3 + 1)
        boxes = max(0, ST - (1 - b)) + 2*max(0, UT - (1 - a))
        per = min(tube, boxes + cost(S0, U0, g, eta)) if normals(g, eta) + offs + tube > T - 0.05 else tube
        best = max(best, offs + per)
    return normals(g, eta) + best
def value(g, deta=0.01):
    S0, U0 = box_S0(g), box_U0(g)
    v0 = E1 + cost(S0, U0, g, 0.0)
    best = v0; arg = 0.0
    run = -np.inf
    for eta in np.arange(deta, 1.6, deta):
        run = max(run, low_at(g, eta - deta), low_at(g, eta))
        if run >= best: break
        hi = E1 + cost(S0, U0, g, eta)
        v = max(hi, run)
        if v < best: best, arg = v, eta
    return best, arg, v0
if __name__ == "__main__":
    print(f"a0={a0:.6f} a={a:.6f} b={b:.4f} T={T:.6f} kappa={kappa:.5f} E1={E1:.5f}")
    gs = np.concatenate([np.arange(0.01, 2.02, 0.02)])
    worst = (-1, None)
    for g in gs:
        v, eta, v0 = value(g)
        if v > worst[0]: worst = (v, g, eta)
        if abs(g*50 - round(g*50)) < 1e-9 and round(g*50) % 5 == 0: print(f"g={g:.2f}  best={v:.5f}  eta*={eta:.2f}  (no split {v0:.5f})", flush=True)
    print("WORST", worst, " vs a0", a0)
