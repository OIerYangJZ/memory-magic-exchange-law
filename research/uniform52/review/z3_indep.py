"""Referee's independent z3 encoding of  not C(a0, T)  (formulas from indep.py, not cont_model.py).

  not C :  exists g in (0, 201/100], y > 0, X >= 0, D, and a piece (tangent | crossing D in [lde, 1-2g])
           such that HIGH(g,y) > T and [ pm > T  or  V4 > T  or (N2 > T and P2 > T) or (N3 > T and P3 > T) ].

No upper bounds on y or X (the authors' y <= 3, X <= 6 are dropped).
usage:
  z3_indep.py check  MODE mu0 c          (mu ranges over (0, mu0])
  z3_indep.py fixed  MODE mu c           (mu fixed)
  z3_indep.py opt    MODE mu             (exact sup of min(HIGH, LOW) at fixed mu via z3 Optimize)
MODE in {cand, full}.
"""
import sys, time
from fractions import Fraction as F
import z3
import indep as I

def _c(x): return x if isinstance(x, z3.ExprRef) else z3.RealVal(x)
def zmax(*xs):
    xs = [_c(x) for x in xs]; r = xs[0]
    for x in xs[1:]: r = z3.If(x > r, x, r)
    return r
def zmin(*xs):
    xs = [_c(x) for x in xs]; r = xs[0]
    for x in xs[1:]: r = z3.If(x < r, x, r)
    return r
O = I.Ops(zmax, zmin, z3.If, z3.And, z3.RealVal(10 ** 6))

def branches(mu, T, mode, g, y, X, D):
    a0 = z3.Q(8, 5) + mu
    M = I.model(O, a0, T, g, y, mode)
    out = []
    for kind in ('tangent', 'cross'):
        pm = M['piece'](None if kind == 'tangent' else D)
        N2, P2, N3, P3, V4 = I.routes(O, M, y, X, pm)
        guard = z3.BoolVal(True) if kind == 'tangent' else z3.And(D >= M['lde'], D <= 1 - 2 * g)
        out.append((kind, guard, pm, V4, (N2, P2), (N3, P3)))
    return M, out

def notC(mu, T, mode):
    g, y, X, D = z3.Reals('g y X D')
    M, br = branches(mu, T, mode, g, y, X, D)
    fails = [z3.And(gd, z3.Or(pm > T, V4 > T, z3.And(N2 > T, P2 > T), z3.And(N3 > T, P3 > T)))
             for (_, gd, pm, V4, (N2, P2), (N3, P3)) in br]
    dom = [g > 0, g <= z3.Q(201, 100), y > 0, X >= 0]
    return z3.And(*dom, M['HIGH'] > T, z3.Or(*fails)), (g, y, X, D)

def Qf(f): f = F(f); return z3.Q(f.numerator, f.denominator)

def check(mu_lo_open, mu_hi, c, mode, fixed=False):
    mu = z3.Real('mu')
    E1 = (8 - 2 * (z3.Q(8, 5) + mu)) / 3
    T = E1 + Qf(c) * mu
    phi, vs = notC(mu, T, mode)
    s = z3.Solver()
    s.add(phi, mu == Qf(mu_hi) if fixed else z3.And(mu > 0, mu <= Qf(mu_hi)))
    r = s.check()
    m = s.model() if r == z3.sat else None
    return r, m, (mu,) + vs

def opt(mu_val, mode):
    """sup over (g,y,X,D,branch) of min(HIGH, branch-LOW) = the continuum exponent (exact)."""
    mu = Qf(mu_val)
    E1 = (8 - 2 * (z3.Q(8, 5) + mu)) / 3
    t = z3.Real('t')
    T = t   # G3' uses T; for the exponent we set T = t (self-consistent threshold)
    g, y, X, D = z3.Reals('g y X D')
    M, br = branches(mu, T, mode, g, y, X, D)
    alts = [z3.And(gd, z3.Or(pm >= t, V4 >= t, z3.And(N2 >= t, P2 >= t), z3.And(N3 >= t, P3 >= t)))
            for (_, gd, pm, V4, (N2, P2), (N3, P3)) in br]
    o = z3.Optimize()
    o.add(g > 0, g <= z3.Q(201, 100), y > 0, X >= 0, M['HIGH'] >= t, z3.Or(*alts), t <= 3)
    h = o.maximize(t)
    r = o.check()
    return r, o.model() if r == z3.sat else None, t, (g, y, X, D), z3.simplify(E1)

if __name__ == '__main__':
    cmd, mode = sys.argv[1], sys.argv[2]
    t0 = time.time()
    if cmd in ('check', 'fixed'):
        mu0, c = F(sys.argv[3]), F(sys.argv[4])
        r, m, vs = check(0, mu0, c, mode, fixed=(cmd == 'fixed'))
        print(f"{cmd} mode={mode} mu{'=' if cmd == 'fixed' else '<='}{mu0} c={c}: {r} ({time.time() - t0:.1f}s)")
        if m is not None:
            print({str(v): m.eval(v) for v in vs})
    elif cmd == 'opt':
        mu = F(sys.argv[3])
        r, m, t, vs, E1 = opt(mu, mode)
        if m is not None:
            tv = m.eval(t); tvF = F(str(tv)) if '/' in str(tv) or str(tv).lstrip('-').isdigit() else None
            E1F = F(8 - 2 * (F(8, 5) + mu), 3)
            print(f"opt mode={mode} mu={mu}: sup = {tv}  (sup-E1)/mu = {float((tvF - E1F) / mu) if tvF is not None else '?'}"
                  f"  [10/33 = {10/33:.6f}]  at", {str(v): str(m.eval(v)) for v in vs}, f"({time.time() - t0:.1f}s)")
        else:
            print('opt', r)
