"""Fast exact certificate engine for the lemma set of App. app:elem (same hypotheses as
research/cert/certificate.py), with an analytic best-XL-for-given-XS search.
Usage:  python fastcert.py a0 a b target      (all rationals, e.g. 179/100)
Checks: a >= max(a0,b), a0 + a >= 2b, 2a+3b > 8, and E(a0) = 2(a-a0)+b + max_g min N < target."""
import sys
from fractions import Fraction as F

def patch(a, b, g, zero):
    S0 = (1 - b) if zero else min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    return S0, min(1 - a, 1 - g)

def N(S0, U0, XL, XS):
    return max(F(0), S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))

def ok(gr, XL, XS):
    R0 = 1 - gr
    return (XS <= XL < R0 and 2 * XL - 1 <= XS
            and XL + XS + min(XS, 2 * XL - R0) + 3 < 0)

DEL = F(1, 10**7)
def best_cell(a, b, gl, gr, zero=False, steps=4000):
    S0, U0 = patch(a, b, gl, zero)
    R0 = 1 - gr
    best = None
    top = min(U0, R0 - DEL)
    for k in range(steps):
        XS = top - F(k, 1000)
        if XS < -4: break
        c1 = min((XS + R0) / 2, (R0 - XS - 3) / 3 - DEL)
        cands = [c1]
        if -3 - 2 * XS - DEL > (XS + R0) / 2:
            cands.append(-3 - 2 * XS - DEL)
        XL = min(max(cands), R0 - DEL, (1 + XS) / 2, S0)   # XL > S0 is pointless
        XL = max(XL, XS)
        if not ok(gr, XL, XS):
            continue
        c = N(S0, U0, XL, XS)
        if best is None or c < best[0]:
            best = (c, XL, XS)
        if c == 0: break
    return best

def certify(a0, a, b, target, nint=804, den=400):
    assert a >= a0 and a >= b >= 0 and 2 * a + 3 * b > 8 and a0 + a >= 2 * b
    E1 = 2 * (a - a0) + b
    worst, where = F(0), None
    c0 = best_cell(a, b, F(0), F(1, den), zero=True)
    assert c0 is not None
    worst, where = c0[0], ('r>R/2', c0)
    grid = [F(k, den) for k in range(0, nint + 1)]
    for gl, gr in zip(grid, grid[1:]):
        c = best_cell(a, b, gl, gr)
        assert c is not None, (gl, gr)
        if c[0] > worst:
            worst, where = c[0], (float(gl), c)
    E = E1 + worst
    return E, E1, worst, where, E < target

if __name__ == "__main__":
    a0, a, b, target = (F(x) for x in sys.argv[1:5])
    E, E1, worst, where, okk = certify(a0, a, b, target)
    print(f"a0={a0} a={a} b={b}: E1={float(E1):.5f} worst={float(worst):.5f} at {where[0]} "
          f"E={float(E):.6f} target={float(target):.6f} CERTIFIED={okk}")
