"""Rigorous exponent model for Result E (all quantities are exponents of R', R' ~ 2^{t/4}).

Parameters: a0 = 4L/t (tube radius eps = R'^{-a0}); box transverse rho = R'^{-a}, along lambda = R'^{-b}.
Constraints (level 1): a >= a0, a >= b >= 0, 2a + 3b > 8.     #boxes: 2(a - a0) + b.
Sphere Y of radius r = R'^{R0}, R0 = 1 - g (g = 0 means r >= R'/2).
Patch extents:  S0 = log e_s,  U0 = log e_u
   g = 0 : S0 = 1 - b
   g > 0 : S0 = min(1 - b, (2 - a - min(g, a))/2, 1 - g)      [sublevel-set lemma, r <= R'/2]
   U0 = min(1 - a, 1 - g)
Cells: tube-frame sub-boxes, along x_L = R'^XL, transverse x_S = R'^XS (both transverse dirs), with
   XS <= XL <= R0 - 0 (2 x_L <= r),  XL - b <= XS,  2XL - 1 <= XS,  1 - a - b <= XS  (radial containment)
Coplanarity (level 2):  XL + XS + min(XS, 2XL - R0) + 3 < 0
Cell count (Crofton / grid-face count):
   N = max(0, S0 - XL, U0 - XS, S0 + U0 - XL - XS, U0 + min(U0, max(S0, (R0 + S0)/2)) - 2XS)
Count exponent for the box: level-2 cost N (level 3 = circles: R'^{o(1)}).
"""
import numpy as np

def patch(a, b, g):
    if g == 0:
        S0 = 1 - b
    else:
        S0 = min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    U0 = min(1 - a, 1 - g)
    return S0, U0

def cell_cost(S0, U0, g, XL, XS):
    R0 = 1 - g
    return max(0.0, S0 - XL, U0 - XS, S0 + U0 - XL - XS,
               U0 + min(U0, max(S0, (R0 + S0) / 2)) - 2 * XS)

def feasible_cell(a, b, g, XL, XS, margin=0.0):
    R0 = 1 - g
    if not (XS <= XL <= R0 - 1e-12): return False
    if XL - b > XS or 2 * XL - 1 > XS or 1 - a - b > XS: return False
    return XL + XS + min(XS, 2 * XL - R0) + 3 < -margin

def level2(a, b, g, h=0.005, margin=0.0):
    S0, U0 = patch(a, b, g)
    best = np.inf
    for XL in np.arange(min(S0, 1 - g) , -3.0, -h):
        # for fixed XL the feasible XS form an interval [lo, XL] in which cost decreases in XS,
        # so take the largest feasible XS
        for XS in np.arange(min(XL, U0), -4.0, -h):
            if feasible_cell(a, b, g, XL, XS, margin):
                best = min(best, cell_cost(S0, U0, g, XL, XS))
                break
    return best

GS = np.concatenate([[0.0], np.linspace(1e-6, 2.5, 251)])
def box_exponent(a, b):
    return max(level2(a, b, g) for g in GS)

def optimise(a0, h=0.01):
    best = (np.inf, None)
    for a in np.arange(a0, a0 + 1.5, h):
        for b in np.arange(0, min(a, 3.0) + 1e-9, h):
            if 2 * a + 3 * b <= 8: continue
            E1 = 2 * (a - a0) + b
            if E1 >= best[0]: break
            e = E1 + box_exponent(a, b)
            if e < best[0]:
                best = (e, (round(a, 3), round(b, 3)))
    return best

if __name__ == "__main__":
    import sys
    for t in [float(x) for x in sys.argv[1:]]:
        e, arg = optimise(4 / t)
        print(f"t={t:.2f}L exponent={e:.4f} (a,b)={arg} bits={e*t/4:.4f}L ratio t/bits={t/max(e*t/4,1e-9):.3f}", flush=True)
