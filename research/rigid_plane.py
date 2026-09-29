"""Plane-only rigid scheme (hyperplane -> plane -> circle -> divisor bound), refined extents.
Exponents in units of log R' (R' = 2^{tau/4}); tube radius R'^{-a0}, a0 = 4L/tau.
Level 1: boxes (transverse R'^{-a'}, along R'^{-b'}), need 2a'+3b' > 8.
Level 2: sphere Y of radius R'^{1-g}; Y cap B has extents
   e_s = min(lambda R', r_Y, sqrt(r_Y rho R'))  -> exponent p0 = max(b', g, (g+a')/2)
   e_u = min(rho R', r_Y)                        -> exponent q0 = max(a', g)
   regions of Y with exponents p >= p0, q >= q0 (x_s = R'^{1-p} >= x_u = R'^{1-q} i.e. p <= q, or swap):
   4 points coplanar if vol_sigma1 * vol_sigma2 < c, vol_sigma1 <= x_long * x_short * min(x_short, x_long^2/r_Y),
   vol_sigma2 <= 8 R'^3.
"""
import numpy as np, sys

def feasible(pl, ps, g):
    # pl = exponent of long side (smaller exponent), ps = short side; sizes R'^{1-pl} >= R'^{1-ps}
    long_, short_ = 1 - pl, 1 - ps
    sag = min(short_, 2 * long_ - (1 - g))
    return long_ + short_ + sag + 3 < 0

def level2(ap, bp, g, grid):
    p0 = max(bp, g, (g + ap) / 2) if g <= ap else max(g, bp)
    q0 = max(ap, g)
    best = np.inf
    for p in grid[grid >= p0]:
        for q in grid[grid >= q0]:
            pl, ps = min(p, q), max(p, q)
            if feasible(pl, ps, g):
                best = min(best, (p - p0) + (q - q0))
                break
    if g > 2:           # tiny sphere: whole sphere coplanar
        best = min(best, 0.0) if feasible(g, g, g) else best
    return best

def count_exponent(a0, grid):
    best = (np.inf, None)
    for ap in grid[grid >= a0][:50]:
        for bp in grid[grid >= max((8 - 2 * ap) / 3, 0) + 1e-9][:40]:
            E1 = 2 * (ap - a0) + bp
            if E1 >= best[0]:
                continue
            E2 = max(level2(ap, bp, g, grid) for g in np.linspace(0, 3.5, 71))
            if E1 + E2 < best[0]:
                best = (E1 + E2, (round(ap, 3), round(bp, 3), round(E2, 3)))
    return best

grid = np.round(np.arange(0, 6.0001, 0.02), 4)
rows = []
for tau in [1.2, 1.5, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8]:
    E, arg = count_exponent(4 / tau, grid)
    rows.append((tau, E * tau / 4))
    print(f"tau={tau:.2f}L  exponent={E:.3f} {arg}  bits<={E*tau/4:.3f}L  (spectral {tau/2:.3f}L)")
bits = dict(rows)
ts = sorted(bits)
# tau* where bits reaches 1 (linear interpolation)
for t1, t2 in zip(ts, ts[1:]):
    if bits[t1] < 1 <= bits[t2]:
        tstar = t1 + (1 - bits[t1]) * (t2 - t1) / (bits[t2] - bits[t1])
        print("bits reach L at tau* ~", round(tstar, 3), "L")
print("min over grid of tau/min(bits,1):", round(min(t / min(b, 1) for t, b in rows), 3))
