"""Re-optimise the Result-E scheme with the extra whole-sphere / region-fibration option:
   #(X' cap Y cap region of sigma1-diameter D) <= R'^{o(1)} (1 + D * r2),  r2 <= 2R',
so the level-2 cost for a sphere of exponent g is min(N_cells(g), max(0, S0(g) + 1)),
where S0 = log e_s (the patch diameter is ~ e_s).  Float search (exact certificate afterwards)."""
import numpy as np, sys
sys.path.insert(0, '../push2'); sys.path.insert(0, '../cert')
from fractions import Fraction as F

def patch(a, b, g, zero):
    S0 = (1 - b) if zero else min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    return S0, min(1 - a, 1 - g)

def ok(gr, XL, XS):
    R0 = 1 - gr
    return XS <= XL < R0 and 2 * XL - 1 <= XS and XL + XS + min(XS, 2 * XL - R0) + 3 < 0

def cellcost(a, b, gl, gr, zero=False, steps=2500, DEL=1e-7):
    S0, U0 = patch(a, b, gl, zero); R0 = 1 - gr
    best = None; top = min(U0, R0 - DEL)
    for k in range(steps):
        XS = top - k / 1000
        if XS < -4: break
        c1 = min((XS + R0) / 2, (R0 - XS - 3) / 3 - DEL)
        cands = [c1]
        if -3 - 2 * XS - DEL > (XS + R0) / 2: cands.append(-3 - 2 * XS - DEL)
        XL = max(min(max(cands), R0 - DEL, (1 + XS) / 2, S0), XS)
        if not ok(gr, XL, XS): continue
        c = max(0, S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))
        if best is None or c < best: best = c
        if c == 0: break
    return best

def level2(a, b, gl, gr, zero=False, use_fib=True):
    c = cellcost(a, b, gl, gr, zero)
    if use_fib:
        S0, _ = patch(a, b, gl, zero)
        c = min(c, max(0.0, S0 + 1))   # region fibration: 1 + e_s * r2, r2 <= 2R'
    return c

def E(a0, a, b, use_fib, n=804, den=400):
    if not (a >= a0 and a >= b and 2 * a + 3 * b > 8 and a0 + a >= 2 * b): return np.inf
    w = level2(a, b, 0.0, 1 / den, True, use_fib)
    gs = [k / den for k in range(n + 1)]
    for gl, gr in zip(gs, gs[1:]):
        w = max(w, level2(a, b, gl, gr, False, use_fib))
    return 2 * (a - a0) + b + w, w

def best_alpha(use_fib):
    best = (0, None)
    for alpha in np.arange(2.20, 2.40, 0.005):
        a0 = 4 / alpha; found = None
        for a in np.arange(a0, a0 + 0.1, 0.005):
            bmin = max((8 - 2 * a) / 3 + 1e-4, 0)
            for b in np.arange(bmin, min(a, (a0 + a) / 2) + 1e-9, 0.005):
                e, w = E(a0, a, b, use_fib)
                if e < a0: found = (round(a, 4), round(b, 4), round(e, 4), round(w, 4)); break
            if found: break
        if found: best = (round(alpha, 3), found)
        else: break
    return best

print("cells only        :", best_alpha(False))
print("cells + fibration :", best_alpha(True))
