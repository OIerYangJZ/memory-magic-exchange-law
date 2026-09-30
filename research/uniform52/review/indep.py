"""Referee's independent continuum model (written from main.tex App. E/F and ../annulus/NOTES.md,
NOT from cont_model.py).  Backend-agnostic: `Ops` supplies max/min/ite/and so the same formulas give
z3 terms (exact) or numpy arrays (float brute force).

Prover cell cost, two modes:
  'cand' : patch fibration + the four candidate cells of the authors (X_S = U0; X_L in
           {S0, (U0+R0)/2, U0, min(S0,(U0+R0)/2)}), re-implemented here;
  'full' : patch fibration + the EXACT minimum over all admissible cells
           (X_S <= X_L <= R0, 2X_L - 1 <= X_S), by enumerating the vertices of the line arrangement
           on which N + (P)_+ is piecewise linear.  The line normals are constant, so each vertex is an
           affine function of (S0, U0, R0, y) and the exact min stays in QF-LRA.
"""
import itertools

class Ops:
    def __init__(self, mx, mn, ite, AND, big):
        self.mx, self.mn, self.ite, self.AND, self.big = mx, mn, ite, AND, big

# arrangement lines a1*XL + a2*XS = c(S0,U0,R0,y)
def _lines(S0, U0, R0, y):
    return [
        (1, 0, S0),                 # XL = S0
        (0, 1, U0),                 # XS = U0
        (1, 1, S0 + U0),            # XL + XS = S0 + U0
        (-1, 1, U0 - S0),           # S0 - XL = U0 - XS
        (-1, 2, 2 * U0 - S0),       # S0 - XL = 2(U0 - XS)
        (1, 2, y - 3),              # Pa = 0   (Pa = 1 + (XL + 2XS - y)/3)
        (3, 1, R0 + y - 3),         # Pb = 0   (Pb = 1 + (3XL + XS - R0 - y)/3)
        (-2, 1, -R0),               # Pa = Pb  (XS = 2XL - R0)
        (-1, 1, 0 * S0),            # XS = XL
        (1, 0, R0),                 # XL = R0
        (-2, 1, -1 + 0 * S0),       # XS = 2XL - 1
    ]

def cell_value(o, S0, U0, R0, y, XL, XS):
    N = o.mx(0 * S0, S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))
    P = 1 + (XL + XS + o.mn(XS, 2 * XL - R0) - y) / 3
    return N + o.mx(0 * S0, P)

def admissible(o, R0, XL, XS):
    return o.AND(XS <= XL, XL <= R0, 2 * XL - 1 <= XS)

def cost(o, S0, U0, R0, y, mode):
    vals = [o.mx(0 * S0, 1 + (S0 + 2 * U0 - y) / 3)]         # patch fibration (D3)
    if mode == 'cand':
        XS = U0
        cands = [(S0, XS), ((XS + R0) / 2, XS), (XS, XS), (o.mn(S0, (XS + R0) / 2), XS)]
    else:
        cands = []
        for (a1, a2, c), (b1, b2, d) in itertools.combinations(_lines(S0, U0, R0, y), 2):
            det = a1 * b2 - a2 * b1
            if det == 0: continue
            cands.append(((c * b2 - a2 * d) / det, (a1 * d - c * b1) / det))
    for XL, XS in cands:
        vals.append(o.ite(admissible(o, R0, XL, XS), cell_value(o, S0, U0, R0, y, XL, XS), o.big))
    return o.mn(*vals)

def model(o, a0, T, g, y, mode):
    """returns dict of the continuum quantities at class g, height y (b = E1 = (8-2a0)/3, a = a0)."""
    a = a0; b = (8 - 2 * a0) / 3; E1 = b
    # Lemma E2 patch of a box (r <= R'/2)
    S0 = o.mn(1 - b, (1 - a + o.mx(1 - g, 1 - a)) / 2, 1 - g)
    U0 = o.mn(1 - a, 1 - g)
    R0 = 1 - g
    HIGH = E1 + cost(o, S0, U0, R0, y, mode)
    # F4: delta' and delta'_B exponents
    lde = o.mx(1 - 2 * a0, 1 - g - a0)
    ldb = o.mx(1 - 2 * a, 1 - g - a)
    UT = o.mn(1 - a0, 1 - g)
    U0b = o.mn(1 - a, 1 - g)
    def piece(D):
        if D is None:   # tangent class
            off = o.mx(0 * g, y + lde + 1)
            ST = o.mn((1 + lde) / 2, 1 - g)
            S0b = o.mn(1 - b, (1 + ldb) / 2, 1 - g)
        else:           # crossing class, Delta ~ R'^D
            off = o.mx(0 * g, y + D + 1)
            ST = o.mn((1 + lde) / 2, lde + (1 - D) / 2, 1 - g)
            S0b = o.mn(1 - b, (1 + ldb) / 2, ldb + (1 - D) / 2, 1 - g)
        tube = o.mx(0 * g, 1 + (ST + 2 * UT - y) / 3)
        boxes = o.mx(0 * g, ST - (1 - b)) + 2 * o.mx(0 * g, UT - (1 - a))
        box = boxes + cost(o, S0b, U0b, R0, y, mode)
        return off + o.mn(tube, box)
    gam = o.mn(g, a0)
    ebar = o.mx(S0, U0)
    G3p = T - 4 + a0 + gam
    return dict(E1=E1, S0=S0, U0=U0, HIGH=HIGH, lde=lde, piece=piece, gam=gam, ebar=ebar, G3p=G3p)

def routes(o, M, y, X, pm):
    E1, gam, ebar, U0, G3p = M['E1'], M['gam'], M['ebar'], M['U0'], M['G3p']
    N2 = o.mx(0 * y, 2 * y - o.mx(X, G3p)) + pm          # F5(a) sector count + A1 dichotomy + P1
    P2 = E1 + o.mx(0 * y, ebar + 1 + X - y)               # P2 per box, R'^{E1} boxes
    N3 = o.mx(0 * y, 3 * y - gam - X) + pm                # rank 3 with H_W kept
    P3 = E1 + o.mx(0 * y, ebar + U0 + 2 + X - y)          # P3 per box
    V4 = o.mx(0 * y, 4 * y - 2 * gam) + pm
    return N2, P2, N3, P3, V4
