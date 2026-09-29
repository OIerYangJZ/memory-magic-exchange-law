"""Exact rational certificate for Theorem thm:elem (elementary tube bound) at a0 = 9/5.

Checks, in exact rational arithmetic, that with a = 91/50, b = 73/50 (so 2a+3b > 8, a >= b,
a0 + a >= 2b) every sphere-radius exponent g in [0, 2.01] admits a cell (XL, XS) satisfying the
hypotheses of Lemma lem:E3 (XS <= XL < 1-g, 2XL-1 <= XS, coplanarity (E3)), with cell-count
exponent N of Lemma lem:E4 at most 0.30501.  Constraints are checked at the right end of each
g-interval (they tighten as g grows) and N at the left end (N is nonincreasing in g); g = 0
(r > R'/2) is checked separately on [0, 1/400]; g > 2.01 is trivial (the whole sphere is coplanar).
Output: E(9/5) = 2(a-a0) + b + max cost <= 1.80501, i.e. kappa = E/4 < 5/11 (rate 11/5).
The extra constraint 1-2a <= XS is not needed by the lemma and only restricts the search."""
from fractions import Fraction as F

import sys
A0 = F(sys.argv[1]) if len(sys.argv) > 1 else F(9, 5)
a, b = F(182, 100), F(146, 100)
assert a >= A0 and a >= b >= 0 and 2 * a + 3 * b > 8 and A0 + a >= 2 * b   # hypotheses of Lemma lem:E1
E1 = 2 * (a - A0) + b

def patch(g, zero=False):
    if zero:
        S0 = 1 - b
    else:
        S0 = min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    return S0, min(1 - a, 1 - g)

def N(S0, U0, g, XL, XS):
    R0 = 1 - g
    return max(F(0), S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))

def ok(g, XL, XS):            # feasibility + coplanarity at a given g (strict)
    R0 = 1 - g
    return (XS <= XL < R0 and 2 * XL - 1 <= XS and 1 - 2 * a <= XS
            and XL + XS + min(XS, 2 * XL - R0) + 3 < 0)

def best_cell(gl, gr, zero=False, step=F(1, 200)):
    """search a rational cell valid on [gl, gr]; return (cost at gl, XL, XS)"""
    S0, U0 = patch(gl, zero)
    best = None
    XL = min(S0, 1 - gr) - F(1, 10**6)
    while XL > -3:
        XS = min(XL, U0)
        while XS > -4:
            if ok(gr, XL, XS):
                c = N(S0, U0, gl, XL, XS)
                if best is None or c < best[0]:
                    best = (c, XL, XS)
                break
            XS -= step
        XL -= step
    return best

TARGET = A0 - E1                      # level-2 budget
worst = F(0)
c0 = best_cell(F(0), F(1, 400), zero=True)     # r > R'/2: g in [0, 1/400] once R' >= 2^400
worst = max(worst, c0[0]); print("r > R'/2 (g in [0,1/400]):", float(c0[0]), "cell", float(c0[1]), float(c0[2]))
grid = [F(k, 400) for k in range(0, 805)]          # g in (0, 2.01], step 1/400
for gl, gr in zip(grid, grid[1:]):
    c = best_cell(gl, gr)
    assert c is not None, (gl, gr)
    worst = max(worst, c[0])
print("E1 = 2(a-a0)+b =", float(E1), "  max level-2 cost =", float(worst))
E = E1 + worst
print("a0 =", A0, " E(a0) <=", E, "=", float(E), "  kappa = E/4 =", float(E / 4), "  5/11 =", float(F(5, 11)))
print("CERTIFIED kappa < 5/11 (rate 11/5):", E / 4 < F(5, 11))
