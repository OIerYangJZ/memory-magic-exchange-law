"""Continuum version of the lemma set of research/annulus (cert_proj.py) near the 5/2 limit.
All functions are written over a generic number type so that the same code builds z3 terms
(verify_z3.py) and evaluates floats (this file's __main__).  Layers, g-intervals and crossing classes are
infinitesimal (their widths -> 0 slowly with R', costing R'^{o(1)}).  Prover choices (cells, eta*) are
explicit; adversary choices (y, X, D, rank, piece) are free variables.

Parameters: a0 = 8/5 + mu, a = a0, b = (8 - 2 a0)/3 (+ o(1)), E1 = b, target T."""
FR = lambda p, q=1: p / q

class Num:
    """backend: float by default; verify_z3 swaps in z3 versions of mx/mn/ite"""
    mx = staticmethod(lambda *xs: max(xs))
    mn = staticmethod(lambda *xs: min(xs))
    ite = staticmethod(lambda c, x, y: x if c else y)
    AND = staticmethod(lambda *cs: all(cs))
    BIG = 1e9

def params(mu, T):
    a0 = FR(8, 5) + mu; a = a0; b = (8 - 2 * a0) / 3
    return dict(mu=mu, a0=a0, a=a, b=b, E1=b, T=T)

def patch_box(P, g):
    mx, mn = Num.mx, Num.mn
    a, b = P['a'], P['b']
    S0 = mn(1 - b, (2 - a - mn(g, a)) / 2, 1 - g)
    U0 = mn(1 - a, 1 - g)
    return S0, U0

def cellval(S0, U0, R0, eta, XL, XS):
    mx, mn, ite, AND = Num.mx, Num.mn, Num.ite, Num.AND
    ok = AND(XS <= XL, XL <= R0, 2 * XL - 1 <= XS)
    N = mx(0, S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))
    Pv = (XL + XS + mn(XS, 2 * XL - R0) - eta) / 3 + 1
    return ite(ok, N + mx(0, Pv), Num.BIG)

def cost_upper(S0, U0, R0, eta):
    """min over the patch fibration and a finite list of explicit candidate cells"""
    mx, mn = Num.mx, Num.mn
    vals = [mx(0, (S0 + 2 * U0 - eta) / 3 + 1)]
    XS = U0
    for XL in (S0, (XS + R0) / 2, XS, mn(S0, (XS + R0) / 2)):
        vals.append(cellval(S0, U0, R0, eta, XL, XS))
    return mn(*vals)

def high(P, g, eta):
    S0, U0 = patch_box(P, g)
    return P['E1'] + cost_upper(S0, U0, 1 - g, eta)

def piece_total(P, g, y, D, tangent):
    """offsets(y) + per-sphere(y) for the tangent class (tangent=True) or crossing depth D"""
    mx, mn = Num.mx, Num.mn
    a0, a, b = P['a0'], P['a'], P['b']
    lde = mx(1 - g - a0, 1 - 2 * a0); ldb = mx(1 - g - a, 1 - 2 * a)
    UT = mn(1 - a0, 1 - g); U0b = mn(1 - a, 1 - g)
    if tangent:
        offc = lde + 1
        ST = mn((1 + lde) / 2, 1 - g)
        S0b = mn(1 - b, (1 + ldb) / 2, 1 - g)
    else:
        offc = D + 1
        ST = mn((1 + lde) / 2, lde + (1 - D) / 2, 1 - g)
        S0b = mn(1 - b, (1 + ldb) / 2, ldb + (1 - D) / 2, 1 - g)
    offs = mx(0, y + offc)
    tube = mx(0, (ST + 2 * UT - y) / 3 + 1)
    boxes = mx(0, ST - (1 - b)) + 2 * mx(0, UT - (1 - a))
    per = mn(tube, boxes + cost_upper(S0b, U0b, 1 - g, y))
    return offs + per

def routes(P, g, y, X, pm):
    """(rank-2 route N, route P), (rank-3 route N, route P), v01, v4 -- as functions of adversary X"""
    mx, mn = Num.mx, Num.mn
    a0, T, E1 = P['a0'], P['T'], P['E1']
    gam = mn(g, a0)
    S0, U0 = patch_box(P, g); Sp = mx(S0, U0)
    G3p = T - 4 + a0 + gam
    N2 = mx(0, 2 * y - mx(X, G3p)) + pm; P2 = E1 + mx(0, Sp + 1 + X - y)
    N3 = mx(0, 3 * y - gam - X) + pm; P3 = E1 + mx(0, Sp + U0 + 2 + X - y)
    v4 = mx(0, 4 * y - 2 * gam) + pm
    return (N2, P2), (N3, P3), pm, v4

if __name__ == "__main__":
    import sys, numpy as np
    mu = float(sys.argv[1]); c = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    E1 = (8 - 2 * (1.6 + mu)) / 3; T = E1 + c * mu
    P = params(mu, T)
    Ds = np.linspace(0, 1, 41)
    Xs = np.linspace(0, 3, 601)
    def low(g, y):
        lde = max(1 - g - P['a0'], 1 - 2 * P['a0']); Dmax = 1 - 2 * g
        tots = [piece_total(P, g, y, None, True)]
        if Dmax > lde: tots += [piece_total(P, g, y, lde + f * (Dmax - lde), False) for f in Ds]
        pm = max(tots)
        best = pm
        for X in Xs:
            (N2, P2), (N3, P3), v01, v4 = routes(P, g, y, X, pm)
            best = max(best, min(N2, P2), min(N3, P3), v4)
        return best
    worst = (-9, None)
    for g in np.arange(0.01, 2.01, 0.01):
        # optimal eta by scan (float), for comparison with the explicit formula
        best = (9, None)
        run = -9
        for eta in np.arange(0, 2.0, 0.002):
            if eta > 0: run = max(run, low(g, eta))
            v = max(high(P, g, eta), run)
            if v < best[0]: best = (v, eta)
            if run > best[0]: break
        if best[0] - T > worst[0]: worst = (best[0] - T, g)
        if abs(g * 20 - round(g * 20)) < 1e-9:
            print(f"g={g:.2f} value-E1={best[0]-E1:+.5f} eta*={best[1]:.4f} formula g-2/5+mu/11={g-0.4+mu/11:.4f}", flush=True)
    print(f"mu={mu} T-E1={T-E1:.5f} (a0-E1={1.6+mu-E1:.5f})  max_g(value - T)={worst[0]:+.6f} at g={worst[1]:.2f}")
