"""Exact value at FIXED class g (z3 Optimize), per LOW branch: tests NOTES section 3's claim that the value
is E1 + (10/33)mu for every g in [2/5, 8/5 - 7mu/3] and that it comes from the rank-3 routes."""
import z3
from fractions import Fraction as F
import indep as I, z3_indep as Z
def opt(mu, gval, keep):
    mu = Z.Qf(mu); t = z3.Real('t'); y, X, D = z3.Reals('y X D'); g = Z.Qf(gval)
    M, br = Z.branches(mu, t, 'cand', g, y, X, D)
    alts = []
    for (_, gd, pm, V4, (N2, P2), (N3, P3)) in br:
        opts = dict(pm=pm >= t, v4=V4 >= t, r2=z3.And(N2 >= t, P2 >= t), r3=z3.And(N3 >= t, P3 >= t))
        alts.append(z3.And(gd, z3.Or(*[opts[k] for k in keep])))
    o = z3.Optimize(); o.add(y > 0, X >= 0, M['HIGH'] >= t, z3.Or(*alts), t <= 3)
    o.maximize(t); o.check(); m = o.model()
    E1 = F(8 - 2 * (F(8, 5) + F(mu.as_fraction())), 3)
    return float((F(str(m.eval(t))) - E1) / F(mu.as_fraction()))
for mu in ('1/10', '1/300'):
    for gv in ('1/2', '3/5', '4/5', '1', '6/5', '7/5', '3/2'):
        r = {k: opt(mu, gv, [k]) for k in ('r2', 'r3', 'v4', 'pm')}
        print(f"mu={mu:5s} g={gv:4s}: (value-E1)/mu per branch  " + "  ".join(f"{k}={v:+.5f}" for k, v in r.items())
              + f"   max={max(r.values()):+.5f}  (10/33={10/33:.5f})", flush=True)
