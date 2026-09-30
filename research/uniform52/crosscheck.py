"""Cross-check cont_model (continuum, candidate cells) against the reviewed exact certificate code
(research/dioph/cert_dioph.py pieces + research/annulus/cert_proj.py assembly) at random points.
With gl = gr = g, yp = yj = y and D_l = D_r = D the certificate's per-layer quantities must equal the
continuum ones, except that cont_model's cost_upper (patch + candidate cells) must be >= the certificate's
float-searched cell cost (it is an upper bound, so the continuum model is at most as strong)."""
import sys, os, random
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dioph'))
import cert_dioph as CD
import cont_model as CM
random.seed(7)
mu = F(1, 300); a0 = F(8, 5) + mu; a = a0; b = (8 - 2 * a0) / 3; T = b + mu / 2
P = dict(mu=float(mu), a0=float(a0), a=float(a), b=float(b), E1=float(b), T=float(T))
cost = CD.Cost(True)
bad = 0; worst_gap = 0.0
for it in range(3000):
    g = F(random.randint(1, 2010), 1000); y = F(random.randint(1, 1300), 1000)
    # --- certificate side (exact) ---
    S0, U0 = CD.patch_box(a, b, g); R0 = 1 - g
    ch, _ = cost(S0, U0, R0, y)
    high_cert = float(b + ch)
    lde = max(1 - g - a0, 1 - 2 * a0); ldb = max(1 - g - a, 1 - 2 * a)
    UT = min(1 - a0, 1 - g); U0b = min(1 - a, 1 - g)
    Dmax = 1 - 2 * g
    pieces = [('t', None)]
    if Dmax > lde: pieces.append(('c', lde + (Dmax - lde) * F(random.randint(0, 1000), 1000)))
    for kind, D in pieces:
        if kind == 't':
            offc = lde + 1; ST = min((1 + lde) / 2, 1 - g); S0b = min(1 - b, (1 + ldb) / 2, 1 - g)
        else:
            offc = D + 1; ST = min((1 + lde) / 2, lde + (1 - D) / 2, 1 - g)
            S0b = min(1 - b, (1 + ldb) / 2, ldb + (1 - D) / 2, 1 - g)
        offs = max(F(0), y + offc); tube = max(F(0), (ST + 2 * UT - y) / 3 + 1)
        boxes = max(F(0), ST - (1 - b)) + 2 * max(F(0), UT - (1 - a))
        cv, _ = cost(S0b, U0b, R0, y)
        tot_cert = float(offs + min(tube, boxes + cv))
        tot_cont = CM.piece_total(P, float(g), float(y), None if D is None else float(D), kind == 't')
        # continuum must be >= certificate (only the cell costs may differ, in the safe direction)
        if tot_cont < tot_cert - 1e-9:
            bad += 1
            cvw = cost(S0b, U0b, R0, y)
            print('PIECE violation', kind, 'g', g, 'y', y, 'D', D, 'cert', tot_cert, 'cont', tot_cont)
            print('   S0b', float(S0b), 'U0b', float(U0b), 'R0', float(R0), 'cert cell cost', float(cvw[0]), cvw[1],
                  'cont cost', CM.cost_upper(float(S0b), float(U0b), float(R0), float(y)),
                  'tube', float(tube), 'boxes', float(boxes), 'offs', float(offs))
        worst_gap = max(worst_gap, tot_cont - tot_cert)
    hc = CM.high(P, float(g), float(y))
    if hc < high_cert - 1e-9: bad += 1
print(f"checked 3000 points: violations of cont >= cert: {bad};  max(cont - cert) = {worst_gap:.4f}")

