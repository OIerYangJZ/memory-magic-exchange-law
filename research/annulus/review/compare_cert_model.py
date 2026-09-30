"""Spot check: the certificate's exact LOW-layer value (loop body of cert_proj.py, copied verbatim)
must be >= the referee's pointwise float model (indep_low.py) at the same (g, y) with dy = 1/100,
because the certificate evaluates every quantity at the weaker endpoint and uses cell witnesses.
A certificate value BELOW the model would indicate a bug (an optimistic endpoint or a missing term)."""
import sys, os
sys.dont_write_bytecode = True   # do not write caches outside review/
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'dioph'))
from cert_dioph import Cost, patch_box, low_pieces
from cert_proj import maxmin, crossings
a0, a, b, T = Fr(400, 249), Fr(400, 249), Fr(1192747, 747000), Fr(319751, 199200)
sys.argv = ['x', '400/249', '400/249', '1192747/747000', '319751/199200', '0.01', '0.0025', '12']
import indep_low as Mdl
den, eden, nD = 400, 100, 12
E1 = 2 * (a - a0) + b
cost = Cost(True)

def cert_low(i, m):
    gl, gr = Fr(i, den), Fr(i + 1, den)
    S0, U0 = patch_box(a, b, gl); R0 = 1 - gr
    gam = min(gl, a0)
    pieces = low_pieces(a0, a, b, gl, gr, nD)
    UT = min(1 - a0, 1 - gl); U0b = min(1 - a, 1 - gl)
    Sp = max(S0, U0)
    G3p = T - 4 + a0 + gam
    yj, yp = Fr(m, eden), Fr(m - 1, eden)
    pm = Fr(-10)
    for (offc, S0b, ST) in pieces:
        offs = max(Fr(0), yj + offc)
        tube = max(Fr(0), (ST + 2 * UT - yp) / 3 + 1)
        per = tube
        boxes = max(Fr(0), ST - (1 - b)) + 2 * max(Fr(0), UT - (1 - a))
        cv, _ = cost(S0b, U0b, R0, yp, need=None)
        if cv is not None: per = min(per, boxes + cv)
        pm = max(pm, offs + per)
    v01 = pm
    v4 = max(Fr(0), 4 * yj - 2 * gam) + pm
    A3 = 3 * yj - gam; B3 = Sp + U0 + 2 - yp
    f3 = lambda X: max(Fr(0), A3 - X) + pm
    g3 = lambda X: E1 + max(Fr(0), B3 + X)
    c3 = [Fr(0), A3, -B3] + crossings([(A3 + pm, -1), (pm, 0)], [(E1 + B3, 1), (E1, 0)])
    v3 = max(maxmin(f3, g3, c3), min(pm, g3(max(A3, Fr(0)) + 10)))
    A2 = 2 * yj; B2 = Sp + 1 - yp
    f2 = lambda X: max(Fr(0), A2 - max(X, G3p)) + pm
    g2 = lambda X: E1 + max(Fr(0), B2 + X)
    c2 = [Fr(0), A2, -B2, max(G3p, Fr(0))] + crossings([(A2 - G3p + pm, 0), (A2 + pm, -1), (pm, 0)], [(E1 + B2, 1), (E1, 0)])
    v2 = maxmin(f2, g2, c2)
    hv, _ = cost(S0, U0, R0, yj, need=None)
    return dict(pm=pm, v2=v2, v3=v3, v4=v4, low=max(v01, v2, v3, v4), high=E1 + hv)

worst = 0; cnt = 0; worstH = 0
for i in list(range(200, 805, 15)):
    for m in range(40, 130, 6):
        c = cert_low(i, m)
        g = i / den; y = m / eden
        d = Mdl.low(g, y, y - 0.01, detail=True)
        mlow = max(d['v01'], d['v2'], d['v3'], d['v4'])
        diff = float(c['low']) - mlow
        S0, U0 = Mdl.box_S0(g), Mdl.box_U0(g)
        mhigh = Mdl.E1 + Mdl.cost(S0, U0, (i + 1) / den, y)       # model HIGH at the sagitta endpoint
        worst = min(worst, diff); worstH = min(worstH, float(c['high']) - mhigh); cnt += 1
print(f"{cnt} (interval, layer) pairs: min over pairs of [cert LOW - model LOW] = {worst:+.2e};  "
      f"min of [cert HIGH - model HIGH] = {worstH:+.2e}   (should be >= -1e-9)")
