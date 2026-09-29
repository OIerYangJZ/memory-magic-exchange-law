"""Finer optimisation of the plane-only scheme near tau*, plus an explicit certificate:
for the chosen (a', b') every sphere-radius exponent g on a fine grid gets explicit (p, q)
and all inequalities are re-checked with a margin."""
import numpy as np
exec(open('rigid_plane.py').read().split('grid = np.round')[0])

grid = np.round(np.arange(0, 6.0001, 0.01), 4)
best = None
for tau in np.arange(2.16, 2.26, 0.01):
    E, arg = count_exponent(4 / tau, grid)
    b = E * tau / 4
    print(f"tau={tau:.2f}L bits<={b:.4f}L {arg}")
    if b >= 1 and best is None:
        best = tau
print("first tau on grid with bits >= L:", best)

# certificate at tau = 2.20 L (bits must be < L there)
tau = 2.20; a0 = 4 / tau
E, (ap, bp, E2) = count_exponent(a0, grid)
margin = 1e-6
assert 2 * ap + 3 * bp > 8 + margin, "level 1"
worst = 0
for g in np.linspace(0, 4, 401):
    p0 = max(bp, g, (g + ap) / 2) if g <= ap else max(g, bp)
    q0 = max(ap, g)
    found = None
    for p in grid[grid >= p0]:
        for q in grid[grid >= q0]:
            pl, ps = min(p, q), max(p, q)
            lo, sh = 1 - pl, 1 - ps
            sag = min(sh, 2 * lo - (1 - g))
            if lo + sh + sag + 3 < -margin:
                c = (p - p0) + (q - q0)
                if found is None or c < found[0]:
                    found = (c, p, q)
                break
    worst = max(worst, found[0])
E1 = 2 * (ap - a0) + bp
print(f"certificate tau={tau}L: a'={ap}, b'={bp}, E1={E1:.3f}, max_g E2={worst:.3f}, "
      f"bits <= {(E1 + worst) * tau / 4:.4f} L  (< 1 required)")
