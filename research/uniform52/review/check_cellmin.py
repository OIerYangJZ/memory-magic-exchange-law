"""Validate indep.cost(mode='full') (vertex enumeration) against a dense (XL, XS) grid, and check
cand >= full (candidate cells are an upper bound).  Random (S0, U0, R0, y) in the ranges that occur."""
import numpy as np, float_brute as FB, indep as I
O = FB.O
rng = np.random.default_rng(1)
XL = np.linspace(-3.5, 1.0, 1801)[:, None]; XS = np.linspace(-3.5, 1.0, 1801)[None, :]
worst_grid_below = 0.0; worst_cand_below = 0.0; gap_max = 0.0
for it in range(300):
    g = rng.uniform(0.001, 2.01); mu = rng.uniform(0, 0.1); a = 1.6 + mu; b = (8 - 2 * a) / 3
    R0 = 1 - g; U0 = min(1 - a, 1 - g)
    S0 = min(1 - b, (1 - a + max(1 - g, 1 - a)) / 2, 1 - g) - rng.uniform(0, 0.8) * (rng.random() < 0.5)
    y = rng.uniform(0, 1.3)
    full = float(I.cost(O, S0, U0, R0, y, 'full')); cand = float(I.cost(O, S0, U0, R0, y, 'cand'))
    ok = (XS <= XL) & (XL <= R0) & (2 * XL - 1 <= XS)
    N = np.maximum.reduce([0 * XL * XS, S0 - XL + 0 * XS, U0 - XS + 0 * XL, S0 + U0 - XL - XS, 2 * (U0 - XS) + 0 * XL])
    P = 1 + (XL + XS + np.minimum(XS, 2 * XL - R0) - y) / 3
    c = np.where(ok, N + np.maximum(0, P), np.inf)
    grid = min(float(c.min()), max(0.0, 1 + (S0 + 2 * U0 - y) / 3))
    worst_grid_below = max(worst_grid_below, full - grid)       # grid finds lower than 'full' => enumeration bug
    worst_cand_below = max(worst_cand_below, full - cand)       # cand below full => not an upper bound
    gap_max = max(gap_max, cand - full)
print(f"300 random cases: max(full - densegrid) = {worst_grid_below:.2e} (must be <= ~0), "
      f"max(full - cand) = {worst_cand_below:.2e} (must be <= 0), max(cand - full) = {gap_max:.4f}")
