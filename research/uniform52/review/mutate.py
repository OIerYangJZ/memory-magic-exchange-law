"""Mutation tests of the authors' encoding (verify_z3.build on cont_model).  Each mutation perturbs one
formula by +-eps in favour of the adversary (expect unsat -> sat at the tight c = 10/33) or of the prover
(expect sat -> unsat at c = 303/1000 if the formula is binding).  Also reports which disjunct the sat model
uses.  Run from research/uniform52/review with ../ on sys.path."""
import sys, os, z3
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import verify_z3 as V          # patches cont_model.Num to z3
import cont_model as CM

orig = dict(high=CM.high, piece_total=CM.piece_total, routes=CM.routes, cost_upper=CM.cost_upper,
            patch_box=CM.patch_box)
EPS = z3.Q(1, 10 ** 5)

def run(label, c, expect):
    phi, vs = V.build(z3.Q(1, 300), c)
    s = z3.Solver(); s.add(phi); r = s.check()
    flag = 'ok ' if str(r) == expect else 'XX '
    print(f"{flag}{label:60s} c={c}: {r} (expected {expect})", flush=True)
    return r

def reset():
    for k, v in orig.items(): setattr(CM, k, v)

def mut_routes(which, delta):
    def R(P, g, y, X, pm):
        (N2, P2), (N3, P3), v01, v4 = orig['routes'](P, g, y, X, pm)
        d = dict(N2=N2, P2=P2, N3=N3, P3=P3); d[which] = d[which] + delta
        return (d['N2'], d['P2']), (d['N3'], d['P3']), v01, v4
    return R

C_TIGHT, C_BELOW = z3.Q(10, 33), z3.Q(303, 1000)
reset(); run('baseline', C_TIGHT, 'unsat'); run('baseline', C_BELOW, 'sat')

# which disjunct does the sat model at c = 303/1000 use?
phi, (mu, g, y, X, D) = V.build(z3.Q(1, 300), C_BELOW)
s = z3.Solver(); s.add(phi); s.check(); m = s.model()
a0 = z3.Q(8, 5) + mu; E1 = (8 - 2 * a0) / 3; T = E1 + C_BELOW * mu
P = dict(mu=mu, a0=a0, a=a0, b=E1, E1=E1, T=T)
for tangent in (True, False):
    tot = CM.piece_total(P, g, y, D, tangent)
    (N2, P2), (N3, P3), _, v4 = CM.routes(P, g, y, X, tot)
    ev = lambda e: float(m.eval(e).as_fraction()) if hasattr(m.eval(e), 'as_fraction') else m.eval(e)
    print(f"   sat model piece={'tangent' if tangent else 'crossing'}: tot-T={ev(tot - T):+.6f} v4-T={ev(v4 - T):+.6f} "
          f"N3-T={ev(N3 - T):+.6f} P3-T={ev(P3 - T):+.6f} N2-T={ev(N2 - T):+.6f} P2-T={ev(P2 - T):+.6f} "
          f"HIGH-T={ev(CM.high(P, g, y) - T):+.6f}")

# adversary-favouring mutations (expect sat at the tight threshold if binding)
for w in ('N3', 'P3', 'N2', 'P2'):
    CM.routes = mut_routes(w, EPS); run(f'{w} + 1e-5', C_TIGHT, 'sat' if w in ('N3', 'P3') else 'unsat'); reset()
CM.high = lambda P, g, eta: orig['high'](P, g, eta) + EPS; run('HIGH + 1e-5', C_TIGHT, 'sat'); reset()
CM.piece_total = lambda P, g, y, D, t: orig['piece_total'](P, g, y, D, t) + EPS; run('pm + 1e-5', C_TIGHT, 'sat'); reset()
# prover-favouring mutations (expect unsat at c = 303/1000 if binding)
for w in ('N3', 'P3'):
    CM.routes = mut_routes(w, -z3.Q(1, 10 ** 3)); run(f'{w} - 1e-3', C_BELOW, 'unsat'); reset()
CM.high = lambda P, g, eta: orig['high'](P, g, eta) - z3.Q(1, 10 ** 3); run('HIGH - 1e-3', C_BELOW, 'unsat'); reset()
# structural mutations
def only_crossing_D_top(P, g, y, D, tangent):
    return orig['piece_total'](P, g, y, D, tangent)
CM.cost_upper = lambda S0, U0, R0, eta: orig['cost_upper'](S0, U0, R0, eta) + EPS
run('cell cost + 1e-5 (HIGH and per-sphere box bound)', C_TIGHT, 'sat'); reset()
# a formula error that should be caught: route P3 constant 2 -> 2.01 (adversary)
def R_bad(P, g, y, X, pm):
    (N2, P2), (N3, P3), v01, v4 = orig['routes'](P, g, y, X, pm)
    return (N2, P2), (N3, P3 + z3.Q(1, 100)), v01, v4
CM.routes = R_bad; run('P3 constant 2 -> 2.01', C_TIGHT, 'sat'); reset()
# a harmless-looking error in G3' (T -> T - 1/10): rank 2 should stay non-binding?
def R_g3(P, g, y, X, pm):
    Q = dict(P); Q['T'] = P['T'] - z3.Q(1, 10)
    return orig['routes'](Q, g, y, X, pm)
CM.routes = R_g3; run("G3' with T - 1/10 (rank-2 route N larger; informational)", C_TIGHT, 'unsat'); reset()
# Python-level min/max on z3 terms must fail loudly, not silently
try:
    max(z3.Real('u'), 0); print('XX Python max on z3 term did NOT raise')
except z3.Z3Exception as e:
    print('ok Python max/min on a z3 term raises (cannot be silently evaluated):', str(e)[:60])
