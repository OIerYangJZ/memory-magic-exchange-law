import numpy as np
def patch(a, b, g, zero):
    S0 = (1 - b) if zero else min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    return S0, min(1 - a, 1 - g)
XS_GRID = np.arange(0.0, -4.0, -0.005)
def cellcost(a, b, gl, gr, zero=False, DEL=1e-7):
    S0, U0 = patch(a, b, gl, zero); R0 = 1 - gr
    XS = XS_GRID[XS_GRID <= min(U0, R0 - DEL)]
    c1 = np.minimum((XS + R0) / 2, (R0 - XS - 3) / 3 - DEL)
    c2 = np.where(-3 - 2 * XS - DEL > (XS + R0) / 2, -3 - 2 * XS - DEL, -np.inf)
    XL = np.maximum(np.minimum.reduce([np.maximum(c1, c2), np.full_like(XS, R0 - DEL), (1 + XS) / 2, np.full_like(XS, S0)]), XS)
    ok = (XS <= XL) & (XL < R0) & (2 * XL - 1 <= XS) & (XL + XS + np.minimum(XS, 2 * XL - R0) + 3 < 0)
    if not ok.any(): return np.inf
    c = np.maximum.reduce([np.zeros_like(XS), S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS)])
    return c[ok].min()
def level2(a, b, gl, gr, zero, fib):
    c = cellcost(a, b, gl, gr, zero)
    if fib: c = min(c, max(0.0, patch(a, b, gl, zero)[0] + 1))
    return c
def E(a0, a, b, fib, den=200, top=2.01):
    w = level2(a, b, 0.0, 1 / den, True, fib)
    gs = np.arange(0, top + 1e-12, 1 / den)
    for gl, gr in zip(gs[:-1], gs[1:]): w = max(w, level2(a, b, gl, gr, False, fib))
    return 2 * (a - a0) + b + w, w
def best(fib):
    res = None
    for alpha in np.arange(2.22, 2.40, 0.005):
        a0 = 4 / alpha; hit = None
        for a in np.arange(a0, a0 + 0.08, 0.01):
            for b in np.arange(max((8 - 2 * a) / 3 + 1e-3, 0), min(a, (a0 + a) / 2), 0.01):
                e, w = E(a0, a, b, fib)
                if e < a0: hit = (round(a, 3), round(b, 3), round(e, 4), round(w, 4)); break
            if hit: break
        if not hit: break
        res = (round(alpha, 3), hit)
    return res
print("cells only       :", best(False), flush=True)
print("cells+fibration  :", best(True), flush=True)
