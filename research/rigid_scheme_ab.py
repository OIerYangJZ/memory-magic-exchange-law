"""Exponent bookkeeping for the rigid determinant scheme on S^3 (all exponents in units of log R',
R' = 2^{tau/4}; tube radius eps = R'^{-a0}, a0 = 4L/tau).

Level 1 (affine d=1 on S^3): box with transverse eps-exponent a' >= a0, along-exponent b':
   5 points affinely dependent if 2a' + 3b' > 8   -> points on a 2-sphere Y = S^3 cap H.
Level 2 on Y (radius r_Y = R'^{1-g}), sub-caps of Y of along/transverse exponents p, q
   (absolute sizes R'^{1-p}, R'^{1-q}); p >= max(b', g), q >= max(a', g):
   (a) high degree, difference coordinates, weighted triangle:  (p-g)(q-g) > (2-g)^2/2
   (b) plane method (4 points coplanar), p <= q:               3p + q > 6 + g
   -> points on a curve of Y; level 3 (circles: divisor bound; others assumed O(R'^o(1))).
Count exponent = [2(a'-a0) + b'] + max_g min_{p,q} [(p - max(b',g)) + (q - max(a',g))].
"""
import numpy as np, sys
USE_A = sys.argv[1] == "A" if len(sys.argv) > 1 else True

def level2(ap, bp, g, grid):
    best = np.inf
    ps = grid[grid >= max(bp, g)]
    qs = grid[grid >= max(ap, g)]
    for p in ps:
        # smallest feasible q for this p
        cands = []
        # (a)
        if USE_A and p - g > 1e-12:
            qa = g + (2 - g) ** 2 / 2 / (p - g)
            cands.append(max(qa, max(ap, g)))
        # (b) with the orientation p<=q or q<=p (symmetric version: 3*min+max > 6+g)
        qb1 = 6 + g - 3 * p          # if p <= q
        if qb1 >= p:
            cands.append(max(qb1, max(ap, g)))
        qb2 = (6 + g - p) / 3        # if q <= p
        qb2 = max(qb2, max(ap, g))
        if qb2 <= p:
            cands.append(qb2)
        if g > 2:                    # tiny sphere: everything coplanar
            cands.append(max(ap, g))
        for q in cands:
            best = min(best, (p - max(bp, g)) + (q - max(ap, g)))
    return best

def count_exponent(a0):
    grid = np.linspace(0, 6, 241)
    best = (np.inf, None)
    for ap in grid[grid >= a0][:60]:
        bp_min = max((8 - 2 * ap) / 3, 0)
        for bp in grid[(grid >= bp_min - 1e-12)][:40]:
            E1 = 2 * (ap - a0) + bp
            if E1 >= best[0]:
                continue
            E2 = max(level2(ap, bp, g, grid) for g in np.linspace(0, 3.2, 33))
            if E1 + E2 < best[0]:
                best = (E1 + E2, (round(ap, 3), round(bp, 3), round(E2, 3)))
    return best

L = 1.0
rows = []
for tau in np.linspace(1.0, 3.2, 23):
    a0 = 4 * L / tau
    E, arg = count_exponent(a0)
    bits = min(E * tau / 4, tau / 2)          # also the spectral bound tau/2 (bits in units of L)
    rows.append((tau, bits))
    print(f"tau={tau:.2f}L  a0={a0:.3f}  exponent={E:.3f} {arg}  bits<= {E*tau/4:.3f}L   spectral {tau/2:.3f}L")
# rate: min over tau of tau / min(bits, 1) restricted to the envelope with b <= 1 (bits per coordinate <= L)
rate = min(t / min(b, 1.0) for t, b in rows if b > 1e-9)
print("rate lower bound from this profile (per bit):", round(rate, 3))
