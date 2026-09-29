import numpy as np
exec(open('rigid_plane.py').read().split('grid = np.round')[0])
def level2(ap, bp, g, grid):
    p0 = max(bp, g, (g + ap) / 2) if g <= ap else max(g, bp)
    q0 = max(ap, g)
    best = np.inf
    for p in grid[grid >= p0]:
        for q in grid[grid >= q0]:
            pl, ps = min(p, q), max(p, q)
            okB = feasible(pl, ps, g)
            okA = (p - g) * (q - g) > (2 - g) ** 2 / 2 if g < 2 else True
            if okA or okB:
                best = min(best, (p - p0) + (q - q0)); break
    return best
grid = np.round(np.arange(0, 6.0001, 0.02), 4)
for tau in (2.2, 2.3, 2.4, 2.5, 2.6):
    E, arg = count_exponent(4 / tau, grid)
    print(f"tau={tau}L bits<={E*tau/4:.3f}L {arg}")
