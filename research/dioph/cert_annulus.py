"""Exact rational certificate for the Diophantine-dichotomy bound  #X_t <= R'^{T+o(1)},  T < a0.

usage:  python cert_dioph.py a0 a b T [den] [eta_den] [nD]
        e.g.  python cert_dioph.py 5/3 5/3 31/20 41/25        (alpha = 4/a0 = 12/5)

Every inequality below is checked in exact rational arithmetic (fractions.Fraction).  Floats are
used only to *find* the witnesses (cell exponents, eta layers); each witness is then rationalised
and re-verified exactly.  The lemmas behind each term are in DIOPH_NOTES.md (numbers D1-D7).

Structure (exponents in log base R'):
  E1   = 2(a-a0) + b + kappa,  kappa = max(0, 2 - a/2 - 3b/4)                     [D1]
  for each g-interval [gl,gr] (grid 1/den, up to 201/100; beyond, the sphere is coplanar):
    either  E1 + cost(gl,gr,eta=0) <= T,
    or an eta* = m/eta_den with
       HIGH  E1 + cost(gl,gr,eta*) <= T                                             [D2,D3]
       LOW   for every layer j=1..m:  normals(eta_j) + max_pieces [offs + per] <= T  [D4-D7]
  zero case (r > R'/2, g < 1/den):  E1 + cost with S0 = 1-b, eta=0 <= T.
"""
from fractions import Fraction as Fr
import sys, math

def fr(s): return Fr(s)

# ---------------------------------------------------------------- exact level-2 cost pieces
def N_cells(S0, U0, XL, XS):
    return max(0, S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS))

def P_cell(XL, XS, R0, eta):
    return (XL + XS + min(XS, 2 * XL - R0) - eta) / 3 + 1

def cell_ok(XL, XS, R0):
    return XS <= XL and XL < R0 and 2 * XL - 1 <= XS

def cell_value(S0, U0, R0, eta, XL, XS):
    """exact cost of the witness cell, or None if inadmissible"""
    if not cell_ok(XL, XS, R0): return None
    return N_cells(S0, U0, XL, XS) + max(0, P_cell(XL, XS, R0, eta))

def patch_dir(S0, U0, eta):
    return max(Fr(0), (S0 + 2 * U0 - eta) / 3 + 1)

# ---------------------------------------------------------------- float witness search
def float_cell_argmin(S0, U0, R0, eta, step=0.002, span=2.2, DEL=1e-6):
    import numpy as np
    S0, U0, R0, eta = float(S0), float(U0), float(R0), float(eta)
    top = min(U0, R0 - DEL)
    XS = top - np.arange(0, span, step)
    best = (math.inf, None)
    for XL in (S0 - U0 + XS, (XS + R0) / 2, S0 + 0 * XS, XS, (1 + XS) / 2, R0 - DEL + 0 * XS,
               eta - 3 - 2 * XS, (eta - 3 - XS + R0) / 3):
        XL = np.maximum(np.minimum(np.minimum(XL, S0), np.minimum(R0 - DEL, (1 + XS) / 2)), XS)
        ok = (XL < R0) & (2 * XL - 1 <= XS + 1e-12)
        N = np.maximum.reduce([0 * XS, S0 - XL, U0 - XS, S0 + U0 - XL - XS, 2 * (U0 - XS)])
        P = (XL + XS + np.minimum(XS, 2 * XL - R0) - eta) / 3 + 1
        c = np.where(ok, N + np.maximum(0, P), np.inf)
        k = int(np.argmin(c))
        if c[k] < best[0]: best = (float(c[k]), (float(XL[k]), float(XS[k])))
    return best

def rational_cell(S0, U0, R0, eta, fw, D=100000):
    """rationalise a float witness and nudge it into the admissible region"""
    XL = Fr(fw[0]).limit_denominator(D); XS = Fr(fw[1]).limit_denominator(D)
    XS = min(XS, U0)
    if XL >= R0: XL = R0 - Fr(1, D)
    if 2 * XL - 1 > XS: XL = (1 + XS) / 2
    if XL < XS: XL = XS
    if XL > S0 and S0 >= XS: XL = S0
    return XL, XS

class Cost:
    """exact upper bound for the level-2 cost; returns (value, witness description)"""
    def __init__(self, allow_patch):
        self.allow_patch = allow_patch
        self.cache = {}
    def __call__(self, S0, U0, R0, eta, need=None):
        key = (S0, U0, R0, eta)
        if key in self.cache: return self.cache[key]
        best = (None, None)
        if self.allow_patch:
            v = patch_dir(S0, U0, eta); best = (v, 'patch')
            if need is not None and v <= need:
                self.cache[key] = best; return best
        fv, fw = float_cell_argmin(S0, U0, R0, eta)
        if fw is not None:
            XL, XS = rational_cell(S0, U0, R0, eta, fw)
            v = cell_value(S0, U0, R0, eta, XL, XS)
            if v is not None and (best[0] is None or v < best[0]): best = (v, ('cell', XL, XS))
        self.cache[key] = best
        return best

# ---------------------------------------------------------------- model pieces (exact)
def patch_box(a, b, g, zero=False):
    S0 = (1 - b) if zero else min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
    return S0, min(1 - a, 1 - g)

def normals(a0, gl, eta, Fstar):
    gam = min(gl, a0)
    G = max(Fr(0), gam + a0 + Fstar - 4)   # annular fibration (Lemma A1), replaces gam + Fstar/2 - 2
    return max(eta, 2 * eta - G, 3 * eta - gam, 4 * eta - 2 * gam)

def low_pieces(a0, a, b, gl, gr, nD):
    """list of (offset exponent without eta, per-box S0, tube ST) for the tangent class and
    crossing classes; offsets use the right end, per-sphere quantities the left end."""
    lde_l = max(1 - gl - a0, 1 - 2 * a0); lde_r = max(1 - gr - a0, 1 - 2 * a0)
    ldb_l = max(1 - gl - a, 1 - 2 * a)
    Dmax = 1 - 2 * gl
    out = []
    ST = min((1 + lde_l) / 2, 1 - gl)
    S0 = min(1 - b, (1 + ldb_l) / 2, 1 - gl)
    out.append((lde_l + 1, S0, ST))
    if Dmax > lde_r:
        for i in range(nD):
            Dl = lde_r + (Dmax - lde_r) * i / nD; Dr = lde_r + (Dmax - lde_r) * (i + 1) / nD
            ST = min((1 + lde_l) / 2, lde_l + (1 - Dl) / 2, 1 - gl)
            S0 = min(1 - b, (1 + ldb_l) / 2, ldb_l + (1 - Dl) / 2, 1 - gl)
            out.append((Dr + 1, S0, ST))
    return out

def main():
    a0, a, b, T = map(fr, sys.argv[1:5])
    den = int(sys.argv[5]) if len(sys.argv) > 5 else 200
    eden = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    nD = int(sys.argv[7]) if len(sys.argv) > 7 else 8
    Fstar = T
    assert a >= a0 and a >= b >= 0 and a0 + a >= 2 * b, "box constraints (D1/E1)"
    allow_patch = 2 * b >= a
    kappa = max(Fr(0), 2 - a / 2 - 3 * b / 4)
    E1 = 2 * (a - a0) + b + kappa
    cost = Cost(allow_patch)
    print(f"a0={a0} a={a} b={b} T={T} (={float(T):.6f}, a0={float(a0):.6f})  kappa={kappa}  E1={E1}={float(E1):.6f}")
    assert T < a0 and E1 <= T, 'need E1 <= T (also covers g > 201/100)'
    # zero case
    S0, U0 = patch_box(a, b, Fr(0), zero=True)
    v, w = cost(S0, U0, 1 - Fr(1, den), Fr(0), need=T - E1)
    assert v is not None and E1 + v <= T, ("zero case", v)
    worst = E1 + v; worst_at = 'zero'
    nint = math.ceil(Fr(201, 100) * den)
    eta_max_used = Fr(0); layers_total = 0
    for i in range(nint):
        gl = Fr(i, den); gr = Fr(i + 1, den)
        S0, U0 = patch_box(a, b, gl); R0 = 1 - gr
        v0, _ = cost(S0, U0, R0, Fr(0), need=T - E1)
        if v0 is not None and E1 + v0 <= T:
            if E1 + v0 > worst: worst, worst_at = E1 + v0, (float(gl), 0)
            continue
        assert gl > 0, f"interval {i} needs the split but gl=0"
        pieces = low_pieces(a0, a, b, gl, gr, nD)
        UT = min(1 - a0, 1 - gl); U0b = min(1 - a, 1 - gl)
        low_run = Fr(-10); ok = False
        for m in range(1, 4 * eden):
            eta = Fr(m, eden); eta_prev = Fr(m - 1, eden)
            # layer m: heights in [R'^eta_prev, R'^eta)
            nm = normals(a0, gl, eta, Fstar)
            worst_piece = Fr(-10)
            for (offc, S0b, ST) in pieces:
                offs = max(Fr(0), eta + offc)
                tube = max(Fr(0), (ST + 2 * UT - eta_prev) / 3 + 1)
                per = tube
                if nm + offs + per > T:
                    boxes = max(Fr(0), ST - (1 - b)) + 2 * max(Fr(0), UT - (1 - a))
                    cv, _ = cost(S0b, U0b, R0, eta_prev, need=T - nm - offs - boxes)
                    if cv is not None: per = min(per, boxes + cv)
                worst_piece = max(worst_piece, offs + per)
            low_run = max(low_run, nm + worst_piece)
            if low_run > T: break
            hv, _ = cost(S0, U0, R0, eta, need=T - E1)
            if hv is not None and E1 + hv <= T:
                ok = True
                val = max(E1 + hv, low_run)
                if val > worst: worst, worst_at = val, (float(gl), float(eta))
                eta_max_used = max(eta_max_used, eta); layers_total += m
                break
        assert ok, f"FAILED on interval [{gl},{gr}]"
    print(f"all {nint} g-intervals certified; worst exponent {worst} = {float(worst):.6f} at {worst_at}")
    print(f"max eta* = {eta_max_used}, total layers checked = {layers_total}")
    print(f"CERTIFIED: #X_t <= R'^(T+o(1)) with T = {T} < a0 = {a0};  rate alpha >= 4/a0 = {4 / a0} = {float(4 / a0):.6f}")

if __name__ == "__main__":
    main()
