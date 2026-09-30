"""Exact (z3, linear real arithmetic) verification of the continuum inequalities near the 5/2 limit.

Claim: for all mu in (0, mu0], with a0 = a = 8/5 + mu, b = E1 = (8 - 2 a0)/3 and T = E1 + c*mu,
there is no (g, y) with 0 < g <= 201/100, y > 0 such that
      HIGH(g, y) > T   and   LOW(g, y) > T.
Since HIGH is nonincreasing in eta, this is equivalent to: for every g the least eta* with HIGH <= T
has LOW(g, y) <= T for all y <= eta*  (the HIGH/LOW split of App. F closes at every class g).
Adversary choices (rank, X, crossing depth D, piece) are existential; prover choices are explicit
(the minimum over the patch fibration and the candidate cells of cont_model.cost_upper), so the negated
claim is a quantifier-free LRA formula.  UNSAT = claim proved (exactly, over Q)."""
import sys, time, z3
import cont_model as CM

def zmax(*xs):
    xs = [x if isinstance(x, z3.ExprRef) else z3.RealVal(x) for x in xs]
    r = xs[0]
    for x in xs[1:]: r = z3.If(x > r, x, r)
    return r
def zmin(*xs):
    xs = [x if isinstance(x, z3.ExprRef) else z3.RealVal(x) for x in xs]
    r = xs[0]
    for x in xs[1:]: r = z3.If(x < r, x, r)
    return r
CM.Num.mx = staticmethod(zmax); CM.Num.mn = staticmethod(zmin)
CM.Num.ite = staticmethod(lambda c, x, y: z3.If(c, x, y))
CM.Num.AND = staticmethod(lambda *cs: z3.And(*cs))
CM.Num.BIG = z3.RealVal(10**6)
CM.FR = lambda p, q=1: z3.Q(p, q)

def build(mu0, c):
    mu, g, y, X, D = z3.Reals('mu g y X D')
    a0 = z3.Q(8, 5) + mu
    E1 = (8 - 2 * a0) / 3
    T = E1 + c * mu
    P = dict(mu=mu, a0=a0, a=a0, b=E1, E1=E1, T=T)
    gam = zmin(g, a0)
    lde = zmax(1 - g - a0, 1 - 2 * a0)
    dom = [mu > 0, mu <= mu0, g > 0, g <= z3.Q(201, 100), y > 0, y <= 3, X >= 0, X <= 6]
    high_fail = CM.high(P, g, y) > T
    cases = []
    for tangent in (True, False):
        tot = CM.piece_total(P, g, y, D, tangent)
        (N2, P2), (N3, P3), _, _ = CM.routes(P, g, y, X, tot)
        v4 = zmax(0, 4 * y - 2 * gam) + tot
        fail = z3.Or(tot > T, v4 > T, z3.And(N3 > T, P3 > T), z3.And(N2 > T, P2 > T))
        if not tangent:
            fail = z3.And(D >= lde, D <= 1 - 2 * g, fail)
        cases.append(fail)
    return z3.And(*dom, high_fail, z3.Or(*cases)), (mu, g, y, X, D)

if __name__ == "__main__":
    from fractions import Fraction as F
    mu0 = F(sys.argv[1]) if len(sys.argv) > 1 else F(1, 300)
    c = F(sys.argv[2]) if len(sys.argv) > 2 else F(1, 2)
    t = time.time()
    phi, vs = build(z3.Q(mu0.numerator, mu0.denominator), z3.Q(c.numerator, c.denominator))
    s = z3.Solver(); s.add(phi)
    r = s.check()
    print(f"mu0={mu0} c={c}: {r}  ({time.time()-t:.1f}s)", flush=True)
    if r == z3.sat:
        m = s.model(); print({str(v): m.eval(v) for v in vs})
    elif r == z3.unsat:
        print("VERIFIED: for all mu in (0, mu0], HIGH/LOW split closes at every class g with T = E1 + c*mu < a0")
