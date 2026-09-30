"""Float brute force of the continuum exponent  sup_{g,y} min(HIGH(g,y), LOW(g,y))  (no z3).
Formulas: indep.py (referee's model).  X: exact sup of min(route N, route P) by bisection (N nonincreasing,
P nondecreasing).  D: grid of nD points on [lde, 1-2g] incl. both ends (also records whether the max over D
is always at an endpoint: D = lde, which equals the tangent class, or D = 1-2g).  Cells: 'full' (exact vertex enumeration) or 'cand'.
usage: float_brute.py mu mode [dg] [dy] [nD] [procs]"""
import sys, numpy as np
from multiprocessing import Pool
import indep as I

def _b(*xs): return np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in xs])
O = I.Ops(lambda *xs: np.maximum.reduce(_b(*xs)), lambda *xs: np.minimum.reduce(_b(*xs)),
          lambda c, x, y: np.where(c, x, y), lambda *cs: np.logical_and.reduce(np.broadcast_arrays(*cs)), 1e9)

def maxmin(N, P, lo=0.0, hi=30.0, it=70):
    """sup_{X>=0} min(N(X), P(X)), vectorised; N nonincreasing, P nondecreasing, both continuous."""
    n0, p0 = N(np.full_like(N.shape_hint, lo)), P(np.full_like(N.shape_hint, lo))
    L = np.full_like(n0, lo); H = np.full_like(n0, hi)
    for _ in range(it):
        M = (L + H) / 2
        f = N(M) - P(M)
        L = np.where(f > 0, M, L); H = np.where(f > 0, H, M)
    v = np.minimum(N(H), P(H))
    return np.where(n0 <= p0, n0, v)

def F_at(mu, T, g, y, mode, nD):
    a0 = 1.6 + mu
    M = I.model(O, a0, T, g, y, mode)
    tang = M['piece'](None)
    pm = tang.copy(); dmax_arg = np.zeros_like(y, dtype=bool)
    lde = float(M['lde']); Dmax = 1 - 2 * g
    deep_is_max = True
    if Dmax >= lde:
        best = np.full_like(y, -9.0); bestD = np.zeros_like(y)
        for D in np.linspace(lde, Dmax, nD):
            v = M['piece'](D)
            upd = v > best + 1e-12
            bestD = np.where(upd, D, bestD); best = np.maximum(best, v)
        vdeep = M['piece'](Dmax)
        deep_is_max = bool(np.all(np.maximum(vdeep, tang) >= best - 1e-12))
        pm = np.maximum(pm, best)
    y_ = y
    class Fn:
        pass
    def mk(fun):
        f = lambda X: fun(X); f.shape_hint = y_; return f
    N2 = mk(lambda X: I.routes(O, M, y_, X, pm)[0]); P2 = mk(lambda X: I.routes(O, M, y_, X, pm)[1])
    N3 = mk(lambda X: I.routes(O, M, y_, X, pm)[2]); P3 = mk(lambda X: I.routes(O, M, y_, X, pm)[3])
    V4 = I.routes(O, M, y_, 0 * y_, pm)[4]
    low = np.maximum.reduce([pm, V4, maxmin(N2, P2), maxmin(N3, P3)])
    return np.minimum(M['HIGH'], low), M['HIGH'], low, deep_is_max

def scan_g(args):
    mu, T, g, mode, dy, nD = args
    y = np.arange(dy, 1.3, dy)
    F, H, L, deep = F_at(mu, T, g, y, mode, nD)
    k = int(np.argmax(F))
    # local refinement in y
    yr = np.linspace(max(1e-9, y[k] - 2 * dy), y[k] + 2 * dy, 801)
    Fr, _, _, deep2 = F_at(mu, T, g, yr, mode, nD)
    j = int(np.argmax(Fr))
    return g, float(Fr[j]), float(yr[j]), deep and deep2

if __name__ == '__main__':
    mu = float(sys.argv[1]); mode = sys.argv[2]
    dg = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002
    dy = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0005
    nD = int(sys.argv[5]) if len(sys.argv) > 5 else 41
    procs = int(sys.argv[6]) if len(sys.argv) > 6 else 12
    E1 = (8 - 2 * (1.6 + mu)) / 3
    T = E1 + (10 / 33) * mu          # enters only through G3'
    gs = np.arange(dg, 2.01 + 1e-12, dg)
    with Pool(procs) as p:
        res = p.map(scan_g, [(mu, T, float(g), mode, dy, nD) for g in gs], chunksize=4)
    res.sort(key=lambda r: -r[1])
    g, v, yv, _ = res[0]
    print(f"mu={mu} mode={mode} dg={dg} dy={dy} nD={nD}: sup min(HIGH,LOW) - E1 = {v - E1:+.7f}"
          f"  -> c = {(v - E1) / mu:.5f}  (10/33 = {10/33:.5f})  at g={g:.4f} y={yv:.5f}")
    print("  max over D always at an endpoint (tangent = D=lde, or D=1-2g):", all(r[3] for r in res))
    top = [r for r in res if r[1] > v - 1e-6]
    print(f"  g-range within 1e-6 of sup: [{min(r[0] for r in top):.3f}, {max(r[0] for r in top):.3f}] ({len(top)} grid points)")
    # profile of value(g) - E1 at a few g
    prof = sorted(res, key=lambda r: r[0])
    print("  profile (g, sup_y min - E1):", " ".join(f"{r[0]:.2f}:{r[1]-E1:+.4f}" for r in prof[::max(1, len(prof)//20)]))
