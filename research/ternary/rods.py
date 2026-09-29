"""Max number of sphere points in a thin rod: for every pair (i,j) at distance <= ell, count points
within delta of the segment [p_i, p_j] extended to length ell.  Report the max, and whether any rod with
>= 3 points exists (3 points in a rod of width delta << spacing would already be a 'near-collinear'
cluster)."""
import pickle, math, numpy as np, sys
S = math.sqrt(2)
N, R = pickle.load(open(sys.argv[1], "rb"))
P = np.array([[a + b * S for (a, b) in x] for x in R]); n = len(P)
r = math.sqrt(N[0] + N[1] * S); sp = r * math.sqrt(4 * math.pi / n)
D = np.linalg.norm(P[:, None] - P[None], axis=2)
for lf, df in [(1.0, 0.1), (1.0, 0.05), (2.0, 0.05), (2.0, 0.02), (4.0, 0.02)]:
    ell, dl = lf * sp, df * sp
    best = 0
    for i in range(n):
        for j in np.where((D[i] <= ell) & (D[i] > 0))[0]:
            u = (P[j] - P[i]) / D[i, j]
            X = P - P[i]; t = X @ u
            perp = np.linalg.norm(X - np.outer(t, u), axis=1)
            cnt = int(np.sum((perp <= dl) & (t >= -ell) & (t <= ell)))
            best = max(best, cnt)
    print(f"  rod length {2*lf:.0f}x spacing, width {df}x spacing: max points = {best}")
