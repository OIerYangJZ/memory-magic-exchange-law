"""Unit tests of the max-min over X in cert_proj.py (lines copied verbatim from its main loop) against
brute force on a fine X grid, and of the closed form  max_theta routeN2 = (2y - max(X, G3'))_+  against a
direct maximisation over theta subject to the A1 failure condition."""
import sys, os, random
sys.dont_write_bytecode = True   # do not write caches outside review/
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'dioph'))
from cert_proj import maxmin, crossings

rng = random.Random(1)
def R(lo, hi): return Fr(rng.randint(int(lo * 1000), int(hi * 1000)), 1000)
Xg = np.linspace(0, 8, 160001)
bad3 = bad2 = 0; under = 0
for trial in range(3000):
    yj = R(0.5, 1.4); yp = yj - Fr(1, 100); gam = R(0.5, 1.61); pm = R(0, 0.5); E1 = R(1.55, 1.62)
    Sp = R(-0.7, -0.4); U0 = R(-0.7, -0.4); T = R(1.59, 1.63); a0 = R(1.6, 1.65)
    G3p = T - 4 + a0 + gam
    # ---- verbatim from cert_proj.py
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
    # ---- brute force (float)
    F3 = np.maximum(0, float(A3) - Xg) + float(pm); Gg3 = float(E1) + np.maximum(0, float(B3) + Xg)
    b3 = np.max(np.minimum(F3, Gg3))
    F2 = np.maximum(0, float(A2) - np.maximum(Xg, float(G3p))) + float(pm); Gg2 = float(E1) + np.maximum(0, float(B2) + Xg)
    b2 = np.max(np.minimum(F2, Gg2))
    if abs(float(v3) - b3) > 1e-4: bad3 += 1
    if abs(float(v2) - b2) > 1e-4: bad2 += 1
    if float(v3) < b3 - 1e-4 or float(v2) < b2 - 1e-4: under += 1
print(f"cert_proj max-min vs brute force: rank-3 mismatches {bad3}/3000, rank-2 mismatches {bad2}/3000, "
      f"cases where the certificate UNDER-estimates: {under}")

# closed form of the rank-2 route-N adversary
th = np.linspace(-6, 0, 600001)
bad = 0
for trial in range(2000):
    y = float(R(0.5, 1.4)); gam = float(R(0.3, 1.61)); T = float(R(1.59, 1.63)); a0 = float(R(1.6, 1.65))
    X = float(R(0, 2.5))
    s2 = np.maximum(th, -a0)
    fail = X + 4 - a0 + s2 > T                   # A1 bound exceeds R'^T
    direct = np.max(np.maximum(0, 2 * y - X - np.maximum(0, th[fail] + gam))) if fail.any() else 0.0
    closed = max(0.0, 2 * y - max(X, T - 4 + a0 + gam))
    if direct > closed + 2e-5: bad += 1
print(f"rank-2 route N: direct max over theta exceeds (2y - max(X,G3'))_+ in {bad}/2000 cases")
