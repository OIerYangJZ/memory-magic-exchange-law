"""Float model of the D1-D7 lemma set plus:
  A1 annular fibration (rank-2 dichotomy constant G3 = gam + a0 + T - 4)
  P1 rank-1: O(1) primitive normals
  P2 rank-2 per-box projection along V-perp:   E1 + (S0 + 1 + X - y)_+
  P3 rank-3 per-box projection along m:        E1 + (S0 + U0 + 2 + X - y)_+
Pointwise in g; adversary maximises over the unknown rank d and heights X (and angle theta2 for d=2)."""
import numpy as np, sys
a0 = float(eval(sys.argv[1])); a = float(eval(sys.argv[2])); b = float(eval(sys.argv[3])); T = float(eval(sys.argv[4]))
FLAGS = set(sys.argv[5].split(',')) if len(sys.argv) > 5 else {'A1','P1','P2','P3'}
kappa = max(0.0, 2 - a/2 - 3*b/4); E1 = 2*(a-a0) + b + kappa
XS_grid = np.linspace(-2.6, 1.0, 901)[:, None]; XL_grid = np.linspace(-2.6, 1.0, 901)[None, :]
def cellcost(S0, U0, g, eta):
    XS, XL = XS_grid, XL_grid
    ok = (XS <= XL) & (XL < 1 - g - 1e-9) & (2*XL - 1 <= XS)
    N = np.maximum.reduce([np.zeros_like(XS*XL), S0 - XL + 0*XS, U0 - XS + 0*XL, S0 + U0 - XL - XS, 2*(U0 - XS) + 0*XL])
    P = (XL + XS + np.minimum(XS, 2*XL - (1 - g)) - eta)/3 + 1
    return np.where(ok, N + np.maximum(0, P), np.inf).min()
_cc = {}
def cost(S0, U0, g, eta):
    k = (round(S0,6), round(U0,6), round(g,6), round(eta,6))
    if k not in _cc:
        c = cellcost(S0, U0, g, eta)
        if 2*b >= a: c = min(c, max(0.0, (S0 + 2*U0 - eta)/3 + 1))
        _cc[k] = c
    return _cc[k]
def box_S0(g): return min(1 - b, (2 - a - min(g, a))/2, 1 - g)
def box_U0(g): return min(1 - a, 1 - g)
def pieces_max(g, eta_n, eta_c):
    """max over tangent/crossing classes of offsets(eta_n) + per-sphere(eta_c)"""
    lde = max(1 - g - a0, 1 - 2*a0); ldb = max(1 - g - a, 1 - 2*a)
    UT = min(1 - a0, 1 - g); U0 = box_U0(g)
    best = -np.inf
    for D in [None] + list(np.linspace(lde, max(lde, 1 - 2*g), 13)):
        if D is None:
            offs = max(0, eta_n + lde + 1); ST = min((1 + lde)/2, 1 - g); S0 = min(1 - b, (1 + ldb)/2, 1 - g)
        else:
            offs = max(0, eta_n + D + 1)
            ST = min((1 + lde)/2, lde + (1 - D)/2, 1 - g); S0 = min(1 - b, (1 + ldb)/2, ldb + (1 - D)/2, 1 - g)
        tube = max(0, (ST + 2*UT - eta_c)/3 + 1)
        boxes = max(0, ST - (1 - b)) + 2*max(0, UT - (1 - a))
        per = min(tube, boxes + cost(S0, U0, g, eta_c))
        best = max(best, offs + per)
    return best
Xs = np.linspace(0, 3, 301)
def low_layer(g, yj, yp):
    """bound for the points in boxes whose sphere has class ~g and height in [R'^yp, R'^yj)"""
    gam = min(g, a0); S0, U0 = box_S0(g), box_U0(g)
    pm = pieces_max(g, yj, yp)
    # rank 0/1
    v01 = (0 if 'P1' in FLAGS else yj) + pm
    # rank 4
    v4 = max(0, 4*yj - 2*gam) + pm
    # rank 3: adversary picks X = log H_W
    n3 = np.maximum(0, 3*yj - gam - Xs) + pm
    if 'P3' in FLAGS:
        pp = E1 + np.maximum(0, S0 + U0 + 2 + Xs - yp)
        v3 = np.max(np.minimum(n3, pp))
    else:
        v3 = n3[0]
    # rank 2: adversary picks X = log H_V, th = log sin theta2 (in [-inf,0]); annulus fibration must fail
    if 'A1' in FLAGS:
        G = max(0.0, gam + a0 + T - 4)
    else:
        G = max(0.0, gam + T/2 - 2)
    if 'P2' in FLAGS:
        # route N normals = max(0, 2y - X - (th+gam)_+); failure region; min over X,th of saving handled by grid
        best = -np.inf
        for th in np.linspace(-a0-0.5, 0, 60):
            s2 = max(th, -a0)
            if 'A1' in FLAGS: Xmin = max(0.0, T - (4 - a0 + s2))
            else: Xmin = max(0.0, T - (4 + 2*s2))
            X = Xs[Xs >= Xmin - 1e-12]
            if len(X) == 0: continue
            nN = np.maximum(0, 2*yj - X - max(0, th + gam)) if 'P1' in FLAGS else np.maximum(yj, 2*yj - X - max(0, th+gam))
            rN = nN + pm
            rP = E1 + np.maximum(0, S0 + 1 + X - yp)
            best = max(best, np.max(np.minimum(rN, rP)))
        v2 = best
    else:
        v2 = max(0 if 'P1' in FLAGS else yj, 2*yj - G) + pm
    return max(v01, v2, v3, v4)
def value(g, deta=0.01):
    S0, U0 = box_S0(g), box_U0(g)
    best = E1 + cost(S0, U0, g, 0.0); arg = 0.0; run = -np.inf
    for eta in np.arange(deta, 2.0, deta):
        run = max(run, low_layer(g, eta, eta - deta))
        if run >= best: break
        v = max(E1 + cost(S0, U0, g, eta), run)
        if v < best: best, arg = v, eta
    return best, arg
if __name__ == "__main__":
    gs = np.arange(0.01, 2.02, 0.02)
    worst = (-1, None, None)
    for g in gs:
        v, eta = value(g)
        if v > worst[0]: worst = (v, g, eta)
    print(f"a0={a0:.5f} b={b} T={T:.5f} E1={E1:.4f} flags={sorted(FLAGS)}  WORST {worst[0]:.5f} at g={worst[1]:.2f} eta*={worst[2]:.2f}   margin a0-worst={a0-worst[0]:+.5f}", flush=True)
