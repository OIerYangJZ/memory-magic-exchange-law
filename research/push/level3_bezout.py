"""Option: at level 2 use high-degree determinants on the (rigid) sphere Y (condition A),
which puts the points of a region on a curve Gamma = Y cap {G=0} of degree <= 2 d2;
then level 3: cut Gamma into pieces of diameter z = R'^{1-p3}; 4 points of a piece are
coplanar once z^2 min(z, z^2/r_Y) R'^3 < c, and a plane meets Gamma in <= 2 deg points
unless the circle Y cap plane lies in {G=0} (then divisor bound).  Curve length inside a
region of sides x_L >= x_S is <= C(x_L + x_S) (Crofton), so #pieces << 1 + x_L/z.
All exponents in units of log R'."""
import numpy as np, sys
sys.path.insert(0, '.')
src = open('rigid_plane_copy.py').read().split('grid = np.round')[0]
exec(src)

def p3_min(g):
    # smallest p3 with 2(1-p3) + min(1-p3, 2(1-p3)-(1-g)) + 3 < 0, and p3 >= g
    for p3 in np.arange(g, 6, 0.005):
        lo = 1 - p3
        if 2 * lo + min(lo, 2 * lo - (1 - g)) + 3 < 0:
            return p3
    return np.inf

P3 = {}
def level2(ap, bp, g, grid):
    p0 = max(bp, g, (g + ap) / 2) if g <= ap else max(g, bp)
    q0 = max(ap, g)
    if g not in P3:
        P3[g] = p3_min(g)
    best = np.inf
    for p in grid[grid >= p0]:
        for q in grid[grid >= q0]:
            pl, ps = min(p, q), max(p, q)
            c = None
            if feasible(pl, ps, g):
                c = (p - p0) + (q - q0)
            elif g < 2 and (p - g) * (q - g) > (2 - g) ** 2 / 2:
                c = (p - p0) + (q - q0) + max(0.0, P3[g] - pl)
            if c is not None:
                best = min(best, c)
                if feasible(pl, ps, g):
                    break
    return best

grid = np.round(np.arange(0, 6.0001, 0.02), 4)
for tau in (2.1, 2.2, 2.25, 2.3):
    E, arg = count_exponent(4 / tau, grid)
    print(f"tau={tau}L bits<={E*tau/4:.4f}L {arg}")
