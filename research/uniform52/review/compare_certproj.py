"""Faithfulness check: the per-layer LOW and HIGH of the reviewed certificate (cert_proj.py lines 63-100,
copied verbatim, with gl = gr = g, yp = yj = y, and pieces = tangent + one crossing class D_l = D_r = D)
against the referee's continuum model (indep.py, exact cells, exact max over X) at random points.
Expected: cert >= model (cert cells are rationalised float witnesses) and equal up to ~1e-6."""
import sys, os, random, numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'dioph'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'annulus'))
from cert_dioph import Cost, patch_box
from cert_proj import maxmin, crossings
import indep as I, float_brute as FB
random.seed(11)
def cert_layer(a0, a, b, T, g, y, D, cost):
    E1 = 2 * (a - a0) + b
    S0, U0 = patch_box(a, b, g); R0 = 1 - g
    gam = min(g, a0)
    lde = max(1 - g - a0, 1 - 2 * a0); ldb = max(1 - g - a, 1 - 2 * a)
    pieces = [(lde + 1, min(1 - b, (1 + ldb) / 2, 1 - g), min((1 + lde) / 2, 1 - g))]
    if D is not None:
        pieces.append((D + 1, min(1 - b, (1 + ldb) / 2, ldb + (1 - D) / 2, 1 - g), min((1 + lde) / 2, lde + (1 - D) / 2, 1 - g)))
    UT = min(1 - a0, 1 - g); U0b = min(1 - a, 1 - g)
    Sp = max(S0, U0); G3p = T - 4 + a0 + gam
    yj = yp = y
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
    return max(v01, v2, v3, v4), E1 + hv
def model_layer(mu, T, g, y, D):
    a0 = 1.6 + mu
    ya = np.array([y])
    M = I.model(FB.O, a0, T, g, ya, 'full')
    pm = M['piece'](None)
    if D is not None: pm = np.maximum(pm, M['piece'](D))
    mk = lambda k: (lambda X: I.routes(FB.O, M, ya, X, pm)[k])
    fs = []
    for k in range(4):
        f = mk(k); f.shape_hint = ya; fs.append(f)
    V4 = I.routes(FB.O, M, ya, 0 * ya, pm)[4]
    low = max(float(pm[0]), float(V4[0]), float(FB.maxmin(fs[0], fs[1])[0]), float(FB.maxmin(fs[2], fs[3])[0]))
    return low, float(M['HIGH'][0])
worst_below = 0.0; worst_above = 0.0; n = 0; worst_at = None
for it in range(1500):
    mu = Fr(random.choice([1, 3, 10, 30, 100]), 1000)
    a0 = Fr(8, 5) + mu; a = a0; b = (8 - 2 * a0) / 3; T = b + mu * Fr(random.randint(0, 60), 100)
    g = Fr(random.randint(1, 2010), 1000); y = Fr(random.randint(1, 1300), 1000)
    lde = max(1 - g - a0, 1 - 2 * a0); Dmax = 1 - 2 * g
    D = lde + (Dmax - lde) * Fr(random.randint(0, 1000), 1000) if Dmax >= lde else None
    cost = Cost(True)
    lc, hc = cert_layer(a0, a, b, T, g, y, D, cost)
    lm, hm = model_layer(float(mu), float(T), float(g), float(y), None if D is None else float(D))
    for c_, m_ in ((lc, lm), (hc, hm)):
        d = float(c_) - m_
        if d < worst_below: worst_below, worst_at = d, (float(mu), float(g), float(y), D)
        worst_above = max(worst_above, d); n += 1
print(f"{n} comparisons (LOW and HIGH): min(cert - model) = {worst_below:.2e} at {worst_at}, max(cert - model) = {worst_above:.2e}")
