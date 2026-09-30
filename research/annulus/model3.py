"""model2 + D1 slicing with kappa > 0, using the constraints that Minkowski's normal l_B imposes:
   transverse tilt <= R'^{kappa-2+a}  -> gamma_eff = max(min(g,a0), 2-a-kappa)
   height <= R'^{ymax},  ymax = kappa - 2 + max(2b, a)
   class g >= gmin = min(2 - a - kappa, b)"""
import numpy as np, sys, importlib
a0 = float(eval(sys.argv[1])); a = float(eval(sys.argv[2])); b = float(eval(sys.argv[3])); T = float(eval(sys.argv[4]))
sys.argv = sys.argv[:5]
import model2 as M
kappa = max(0.0, 2 - a/2 - 3*b/4); E1 = 2*(a-a0) + b + kappa
M.E1 = E1
gmin = min(2 - a - kappa, b); ymax = kappa - 2 + max(2*b, a); gt = 2 - a - kappa
_orig_low = M.low_layer
def low_layer(g, yj, yp):
    # emulate gamma_eff by evaluating the normal counts with g' = max(g, gt) where only gamma enters
    return _low(g, yj, yp)
def _low(g, yj, yp):
    gam = max(min(g, a0), gt); S0, U0 = M.box_S0(g), M.box_U0(g)
    pm = M.pieces_max(g, yj, yp)
    v01 = pm; v4 = max(0, 4*yj - 2*gam) + pm
    n3 = np.maximum(0, 3*yj - gam - M.Xs) + pm
    pp = E1 + np.maximum(0, max(S0,U0) + U0 + 2 + M.Xs - yp); v3 = np.max(np.minimum(n3, pp))
    best = -np.inf
    for th in np.linspace(-a0-0.5, 0, 60):
        s2 = max(th, -a0); Xmin = max(0.0, T - (4 - a0 + s2))
        X = M.Xs[M.Xs >= Xmin - 1e-12]
        if len(X) == 0: continue
        rN = np.maximum(0, 2*yj - X - max(0, th + gam)) + pm
        rP = E1 + np.maximum(0, max(S0,U0) + 1 + X - yp)
        best = max(best, np.max(np.minimum(rN, rP)))
    return max(v01, best, v3, v4)
def value(g, deta=0.01):
    if g < gmin: return -np.inf, None
    S0, U0 = M.box_S0(g), M.box_U0(g)
    best = E1 + M.cost(S0, U0, g, 0.0); arg = 0.0; run = -np.inf
    for eta in np.arange(deta, 2.0, deta):
        run = max(run, _low(g, eta, eta - deta))
        if run >= best: break
        hi = E1 + M.cost(S0, U0, g, eta) if eta < ymax else -np.inf   # no slice above ymax
        v = max(hi, run)
        if v < best: best, arg = v, eta
        if eta >= ymax: break
    return best, arg
if __name__ == "__main__":
    worst = (-1, None, None)
    for g in np.arange(0.01, 2.02, 0.04):
        v, eta = value(g)
        if v > worst[0]: worst = (v, g, eta)
    print(f"a0={a0:.4f} a={a:.3f} b={b:.3f} kappa={kappa:.3f} E1={E1:.4f} gmin={gmin:.3f} ymax={ymax:.3f}  WORST {worst[0]:.4f} at g={worst[1]:.2f} eta={worst[2]}  a0-worst={a0-worst[0]:+.4f}", flush=True)
