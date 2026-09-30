"""Exact continuum exponent restricted to one LOW branch at a time (z3 Optimize, referee's model)."""
import sys, z3
from fractions import Fraction as F
import indep as I, z3_indep as Z
def opt(mu, mode, keep):
    mu = Z.Qf(mu); t = z3.Real('t'); g, y, X, D = z3.Reals('g y X D')
    M, br = Z.branches(mu, t, mode, g, y, X, D)
    alts = []
    for (_, gd, pm, V4, (N2, P2), (N3, P3)) in br:
        opts = dict(pm=pm >= t, v4=V4 >= t, r2=z3.And(N2 >= t, P2 >= t), r3=z3.And(N3 >= t, P3 >= t))
        alts.append(z3.And(gd, z3.Or(*[opts[k] for k in keep])))
    o = z3.Optimize(); o.add(g > 0, g <= z3.Q(201, 100), y > 0, X >= 0, M['HIGH'] >= t, z3.Or(*alts), t <= 3)
    o.maximize(t); r = o.check(); m = o.model()
    E1 = F(8 - 2 * (F(8, 5) + F(mu.as_fraction())), 3)
    tv = F(str(m.eval(t)))
    return float((tv - E1) / F(mu.as_fraction())), {str(v): round(float(F(str(m.eval(v)))), 5) for v in (g, y, X)}
for mu in ('1/300', '1/10'):
    for keep in (['pm'], ['v4'], ['r2'], ['r3'], ['pm', 'v4', 'r2', 'r3']):
        c, at = opt(mu, 'cand', keep)
        print(f"mu={mu} branches={'+'.join(keep):14s} (sup-E1)/mu = {c:+.6f} at {at}", flush=True)
