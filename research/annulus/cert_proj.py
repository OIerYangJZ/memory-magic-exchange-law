"""Exact rational certificate for  #X_t <= R'^{T+o(1)},  T < a0,  with the lemma set of
research/dioph (D1-D7, App. F of the paper) plus four new ingredients (NOTES.md in this directory):

  A1  annular fibration: #X_t <= R'^{o(1)}(1 + H_V R'^4 eps (sin th2 + eps))       [replaces F5(b)]
  P1  rank <= 1 (and the one-line case of the sector count): O(1) primitive normals
  P2  rank 2, per box: projection along V-perp,   #(X_t cap B) <= R'^{o(1)}(1 + e_s R' H_V / h(n_B))
  P3  rank 3, per box: projection along m,        #(X_t cap B) <= R'^{o(1)}(1 + e_s e_u R'^2 H_W / h(n_B))

usage:  python cert_proj.py a0 a b T [den] [eta_den] [nD]
Every inequality is exact over Q; floats only propose cell witnesses, which are rationalised and
re-checked.  kappa = 0 is required (Lemma E1 as proved in the paper: 2a + 3b > 8).
"""
from fractions import Fraction as Fr
import sys, math, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'dioph'))
from cert_dioph import Cost, patch_box, low_pieces

def maxmin(f, g, cands):
    """max over X >= 0 of min(f(X), g(X)) for f nonincreasing, g nondecreasing, both piecewise
    linear with breakpoints/crossings among cands (and the limit X -> infinity handled by caller)."""
    best = None
    for X in cands:
        if X < 0: continue
        v = min(f(X), g(X))
        if best is None or v > best: best = v
    return best

def crossings(f_pieces, g_pieces):
    """candidate crossing points of affine pieces  c1 - X  and  c2 + X  etc.  f_pieces/g_pieces:
    lists of (const, slope)."""
    out = []
    for (cf, sf) in f_pieces:
        for (cg, sg) in g_pieces:
            if sf != sg: out.append((cg - cf) / (sf - sg))
    return out

def main():
    a0, a, b, T = map(Fr, sys.argv[1:5])
    den = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    eden = int(sys.argv[6]) if len(sys.argv) > 6 else 100
    nD = int(sys.argv[7]) if len(sys.argv) > 7 else 12
    assert a >= a0 and a >= b >= 0 and a0 + a >= 2 * b, "box constraints (E1)"
    assert 2 * a + 3 * b > 8, "kappa = 0 required (Lemma E1)"
    assert 2 * b >= a, "patch fibration (D3) and E2 geometry use 2b >= a"
    E1 = 2 * (a - a0) + b
    assert T < a0 and E1 <= T
    cost = Cost(True)
    print(f"a0={a0} a={a} b={b} T={T} (={float(T):.6f}, a0={float(a0):.6f})  E1={E1}={float(E1):.6f}")
    S0z, U0z = patch_box(a, b, Fr(0), zero=True)
    v, _ = cost(S0z, U0z, 1 - Fr(1, den), Fr(0), need=T - E1)
    assert v is not None and E1 + v <= T, "zero case"
    worst, worst_at = E1 + v, 'zero'
    nint = math.ceil(Fr(201, 100) * den)
    eta_max = Fr(0); nlayers = 0
    for i in range(nint):
        gl, gr = Fr(i, den), Fr(i + 1, den)
        S0, U0 = patch_box(a, b, gl); R0 = 1 - gr
        v0, _ = cost(S0, U0, R0, Fr(0), need=T - E1)
        if v0 is not None and E1 + v0 <= T:
            if E1 + v0 > worst: worst, worst_at = E1 + v0, (float(gl), 0)
            continue
        assert gl > 0
        gam = min(gl, a0)
        pieces = low_pieces(a0, a, b, gl, gr, nD)
        UT = min(1 - a0, 1 - gl); U0b = min(1 - a, 1 - gl)
        Sp = max(S0, U0)                      # projected patch diameter exponent
        G3p = T - 4 + a0 + gam                # annular-fibration threshold (A1)
        run = Fr(-10); ok = False
        for m in range(1, 4 * eden):
            yj, yp = Fr(m, eden), Fr(m - 1, eden)
            # offsets + per-sphere, maximised over tangent/crossing classes
            pm = Fr(-10)
            for (offc, S0b, ST) in pieces:
                offs = max(Fr(0), yj + offc)
                tube = max(Fr(0), (ST + 2 * UT - yp) / 3 + 1)
                per = tube
                boxes = max(Fr(0), ST - (1 - b)) + 2 * max(Fr(0), UT - (1 - a))
                cv, _ = cost(S0b, U0b, R0, yp, need=None)
                if cv is not None: per = min(per, boxes + cv)
                pm = max(pm, offs + per)
            # rank 0/1 (P1)
            v01 = pm
            # rank 4
            v4 = max(Fr(0), 4 * yj - 2 * gam) + pm
            # rank 3: route N  (3y - gam - X)_+ + pm   vs  route P3  E1 + (Sp + U0 + 2 + X - yp)_+
            A3 = 3 * yj - gam; B3 = Sp + U0 + 2 - yp
            f3 = lambda X: max(Fr(0), A3 - X) + pm
            g3 = lambda X: E1 + max(Fr(0), B3 + X)
            c3 = [Fr(0), A3, -B3] + crossings([(A3 + pm, -1), (pm, 0)], [(E1 + B3, 1), (E1, 0)])
            v3 = max(maxmin(f3, g3, c3), min(pm, g3(max(A3, Fr(0)) + 10)))
            # rank 2: route N  (2y - max(X, G3p))_+ + pm   vs  route P2  E1 + (Sp + 1 + X - yp)_+
            A2 = 2 * yj; B2 = Sp + 1 - yp
            f2 = lambda X: max(Fr(0), A2 - max(X, G3p)) + pm
            g2 = lambda X: E1 + max(Fr(0), B2 + X)
            c2 = [Fr(0), A2, -B2, max(G3p, Fr(0))] + crossings([(A2 - G3p + pm, 0), (A2 + pm, -1), (pm, 0)], [(E1 + B2, 1), (E1, 0)])
            v2 = maxmin(f2, g2, c2)
            low = max(v01, v2, v3, v4)
            run = max(run, low)
            if run > T: break
            hv, _ = cost(S0, U0, R0, yj, need=T - E1)
            if hv is not None and E1 + hv <= T:
                ok = True
                val = max(E1 + hv, run)
                if val > worst: worst, worst_at = val, (float(gl), float(yj))
                eta_max = max(eta_max, yj); nlayers += m
                break
        assert ok, f"FAILED on interval [{gl},{gr}]  (run={float(run):.5f})"
    print(f"all {nint} g-intervals certified; worst exponent {worst} = {float(worst):.6f} at {worst_at}")
    print(f"max eta* = {eta_max}, layers checked = {nlayers}")
    print(f"CERTIFIED: #X_t <= R'^(T+o(1)), T = {T} < a0 = {a0};  alpha >= 4/a0 = {4 / a0} = {float(4 / a0):.6f}")

if __name__ == "__main__":
    main()
