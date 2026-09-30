"""Referee's independent float model of the full assembly with A1 + P1 + P2 + P3 (written from
main.tex App. E/F and research/annulus/NOTES.md, not from cert_proj.py / model2.py).

Differences from the authors' float model:
  * the cell cost (Lemma E4 + D2) is minimised EXACTLY (vertex enumeration of the convex piecewise-
    linear objective on each of the two sagitta regions), not on a grid;
  * everything is pointwise in g (no interval endpoints); layers of width dy (default 0.002) with
    normals/offsets at the top and per-sphere / route P at the bottom; crossing depth D on a grid of
    ND points plus both ends;
  * the rank-2 route N is maximised directly over (X, theta) subject to the A1 failure condition,
    i.e. the formula 2y - max(X, G3') is NOT used (it is re-derived and cross-checked);
  * the max-min over X uses bisection on the crossing (both routes are monotone and continuous).

usage: indep_low.py alpha [dy] [dg] [ND]      (a0 = a = 4/alpha, b = (8-2a0)/3 + 1e-4, T = a0 - 1e-4)
   or: indep_low.py a0 a b T [dy] [dg] [ND]   (explicit)
prints sup_g value(g), where value(g) = min over y* of max(HIGH(y*), LOW layers <= y*), vs T and a0."""
import numpy as np, sys, math

args = sys.argv[1:]
if len(args) >= 4 and all(('/' in s or '.' in s) for s in args[:4]) and float(eval(args[1])) > 1:
    a0, a, b, T = (float(eval(s)) for s in args[:4]); rest = args[4:]
else:
    al = float(eval(args[0])); a0 = 4 / al; a = a0; b = (8 - 2 * a0) / 3 + 1e-4; T = a0 - 1e-4; rest = args[1:]
dy = float(rest[0]) if len(rest) > 0 else 0.002
dg = float(rest[1]) if len(rest) > 1 else 0.0025
NDg = int(rest[2]) if len(rest) > 2 else 24
assert 2 * a + 3 * b > 8 and a >= a0 and a >= b and a0 + a >= 2 * b and 2 * b >= a
E1 = 2 * (a - a0) + b

# ------------------------------------------------------------ exact cell minimisation
def cell_min(S0, U0, g, y):
    """min over cells X_S <= X_L < 1-g, 2X_L - 1 <= X_S of N + max(0,P)  (E4, D2)"""
    R0 = 1 - g
    Ns = np.array([[0, 0, 0], [-1, 0, S0], [0, -1, U0], [-1, -1, S0 + U0], [0, -2, 2 * U0]], float)
    best = math.inf
    for reg in (0, 1):
        # P = 1 + (XL + XS + min(XS, 2XL - R0) - y)/3 ; region 0: XS <= 2XL-R0 ; region 1: XS >= 2XL-R0
        P = np.array([1 / 3, 2 / 3, 1 - y / 3]) if reg == 0 else np.array([1.0, 1 / 3, 1 - (R0 + y) / 3])
        funcs = [Ni + s * P for Ni in Ns for s in (0, 1)]
        lines = []
        for i in range(len(funcs)):
            for j in range(i + 1, len(funcs)):
                d = funcs[i] - funcs[j]
                if abs(d[0]) + abs(d[1]) > 1e-15: lines.append((d[0], d[1], -d[2]))
        cons = [(-1.0, 1.0, 0.0),            # XS - XL <= 0
                (1.0, 0.0, R0 - 1e-9),       # XL <= R0
                (2.0, -1.0, 1.0),            # 2XL - XS <= 1
                ]
        regc = (-2.0, 1.0, -R0) if reg == 0 else (2.0, -1.0, R0)   # reg0: XS - 2XL <= -R0
        allc = cons + [regc]
        lines += [c for c in allc]
        L = np.array(lines)
        i, j = np.triu_indices(len(L), 1)
        A = np.stack([L[i, :2], L[j, :2]], 1)          # (k,2,2)
        rhs = np.stack([L[i, 2], L[j, 2]], 1)
        detA = A[:, 0, 0] * A[:, 1, 1] - A[:, 0, 1] * A[:, 1, 0]
        ok = np.abs(detA) > 1e-12
        A, rhs, detA = A[ok], rhs[ok], detA[ok]
        XL = (rhs[:, 0] * A[:, 1, 1] - rhs[:, 1] * A[:, 0, 1]) / detA
        XS = (A[:, 0, 0] * rhs[:, 1] - A[:, 1, 0] * rhs[:, 0]) / detA
        feas = np.ones_like(XL, bool)
        for (p, q, r) in allc: feas &= p * XL + q * XS <= r + 1e-9
        XL, XS = XL[feas], XS[feas]
        if len(XL) == 0: continue
        N = np.max(Ns[:, 0:1] * XL + Ns[:, 1:2] * XS + Ns[:, 2:3], axis=0)
        Pv = P[0] * XL + P[1] * XS + P[2]
        best = min(best, float(np.min(N + np.maximum(0, Pv))))
    return best

_cache = {}
def cost(S0, U0, g, y):
    k = (round(S0, 9), round(U0, 9), round(g, 9), round(y, 9))
    if k not in _cache:
        patch = max(0.0, 1 + (S0 + 2 * U0 - y) / 3)          # F3(b)/D3, needs 2b >= a (asserted)
        _cache[k] = min(patch, cell_min(S0, U0, g, y))
    return _cache[k]

# ------------------------------------------------------------ per-g geometry
def l_eps(g): return max(1 - 2 * a0, 1 - g - a0)
def l_rho(g): return max(1 - 2 * a, 1 - g - a)
def box_S0(g): return min(1 - b, 0.5 * (1 + l_rho(g)), 1 - g)
def box_U0(g): return min(1 - a, 1 - g)

def pm(g, yn, yc):
    """max over tangent / crossing classes of offsets(yn) + per-sphere(yc)  (F4, F2, F3, E4)"""
    le, lr = l_eps(g), l_rho(g)
    UT = min(1 - a0, 1 - g); U0 = box_U0(g)
    ST0 = min(0.5 * (1 + le), 1 - g); S00 = box_S0(g)
    def per(ST, S0):
        tube = max(0.0, 1 + (ST + 2 * UT - yc) / 3)
        boxes = max(0.0, ST - (1 - b)) + 2 * max(0.0, UT - (1 - a))
        return min(tube, boxes + cost(S0, U0, g, yc))
    best = max(0.0, yn + le + 1) + per(ST0, S00)          # tangent class
    Dmax = 1 - 2 * g
    if Dmax > le:
        for D in np.linspace(le, Dmax, NDg + 1):
            ST = min(ST0, le + 0.5 * (1 - D)); S0 = min(S00, lr + 0.5 * (1 - D))
            best = max(best, max(0.0, yn + D + 1) + per(ST, S0))
    return best

def maxmin(fN, fP, Xhi=6.0):
    """max over X >= 0 of min(fN(X), fP(X)), fN nonincreasing, fP nondecreasing, continuous"""
    if fN(0.0) <= fP(0.0): return fN(0.0)
    lo, hi = 0.0, Xhi
    while fN(hi) > fP(hi): hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if fN(mid) > fP(mid): lo = mid
        else: hi = mid
    return min(fN(lo), fP(lo), fN(hi), fP(hi)) if False else max(min(fN(lo), fP(lo)), min(fN(hi), fP(hi)))

THG = np.concatenate([np.linspace(-4.0, 0.0, 801)])
def routeN2(y, gam, X, pmv):
    """rank 2, route N: max over theta (log sin th2 <= 0) with the A1 bound failing:
       X + 4 - a0 + max(th, -a0) >= T ;  count exponent  (2y - X - max(0, th + gam))_+   (D6(a) + P1)"""
    s2 = np.maximum(THG, -a0)
    okm = X + 4 - a0 + s2 >= T - 1e-12
    if not okm.any(): return pmv            # A1 succeeds for every theta: nothing to count
    th = THG[okm]
    # max over admissible theta: the count decreases in th + gam, so take the smallest admissible th
    val = np.max(np.maximum(0.0, 2 * y - X - np.maximum(0.0, th + gam)))
    # exact endpoint (continuous theta): smallest admissible theta
    c = T + a0 - 4 - X
    thmin = -np.inf if c < -a0 else c
    val = max(val, max(0.0, 2 * y - X - max(0.0, thmin + gam)))
    return val + pmv

def low(g, yn, yc, detail=False):
    gam = min(g, a0)
    S0, U0 = box_S0(g), box_U0(g); eb = max(S0, U0)
    p = pm(g, yn, yc)
    v01 = p                                                              # P1
    v4 = max(0.0, 4 * yn - 2 * gam) + p
    v3 = maxmin(lambda X: max(0.0, 3 * yn - gam - X) + p,                # route N, rank 3 (H_W kept)
                lambda X: E1 + max(0.0, eb + U0 + 2 + X - yc))           # route P3
    v2 = maxmin(lambda X: routeN2(yn, gam, X, p),                        # route N, rank 2 (A1 + D6a + P1)
                lambda X: E1 + max(0.0, eb + 1 + X - yc))                # route P2
    if detail: return dict(pm=p, v01=v01, v2=v2, v3=v3, v4=v4)
    return max(v01, v2, v3, v4)

def value(g, zero=False):
    S0 = (1 - b) if zero else box_S0(g); U0 = box_U0(g)
    geff = max(g, dg) if zero else g
    best = E1 + cost(S0, U0, geff, 0.0); arg = 0.0
    if zero or g <= 0: return best, arg
    run = -math.inf; y = 0.0
    while y < 3.0:
        y += dy
        run = max(run, low(g, y, y - dy))
        if run >= best: break
        hi = E1 + cost(S0, U0, g, y)
        v = max(hi, run)
        if v < best: best, arg = v, y
    return best, arg

if __name__ == "__main__":
    print(f"a0={a0:.6f} a={a:.6f} b={b:.6f} T={T:.6f} E1={E1:.6f}  (alpha=4/a0={4/a0:.5f})  dy={dy} dg={dg} ND={NDg}")
    z, _ = value(0.0, zero=True)
    worst = (z, 'zero', 0)
    gs = np.arange(dg, 2.01 + 1e-12, dg)
    for g in gs:
        v, y = value(float(g))
        if v > worst[0]: worst = (v, float(g), y)
    print(f"sup_g value = {worst[0]:.6f} at g={worst[1]}, y*={worst[2]:.3f};  T - sup = {T - worst[0]:+.6f};  a0 - sup = {a0 - worst[0]:+.6f}")
