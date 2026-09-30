"""Referee check of P3 step 3 / P2 step 4 (float, Monte Carlo): for a four-box Q with sides
(e_u, e_s, 2e_u, 2e_u) in a random frame of R^4, a random 3-space W1, and a random 2-plane S in W1
meeting P_W1(Q), the section S cap P_W1(Q) has area <= C max(e_s,e_u) e_u and diameter <= C max(e_s,e_u);
for a random line in a random 2-plane V1, the section of P_V1(Q) has length <= C max(e_s,e_u)."""
import numpy as np
rng = np.random.default_rng(2)

def hull_area(P):
    P = np.unique(P, axis=0); P = P[np.lexsort((P[:, 1], P[:, 0]))]
    if len(P) < 3: return 0.0
    def cross(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in P[::-1]:
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    H = np.array(lo[:-1] + up[:-1])
    x, y = H[:, 0], H[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

worstA = worstD = worstL = 0
for trial in range(400):
    eu = 1.0; es = 10 ** rng.uniform(-1, 3); eb = max(es, eu)
    F, _ = np.linalg.qr(rng.standard_normal((4, 4)))
    sides = np.array([eu, es, 2 * eu, 2 * eu])
    n = 400000
    Q = (rng.uniform(-0.5, 0.5, (n, 4)) * sides) @ F.T
    # P3: W1 random 3-space, plane S in W1 through a random point of P(Q)
    W, _ = np.linalg.qr(rng.standard_normal((4, 3)))
    Y = Q @ W                                    # coordinates in W1 (3D)
    nv = rng.standard_normal(3); nv /= np.linalg.norm(nv)
    c = Y[rng.integers(n)] @ nv
    tau = 0.01 * eu
    sel = np.abs(Y @ nv - c) < tau
    if sel.sum() > 50:
        B = np.linalg.svd(np.eye(3) - np.outer(nv, nv))[0][:, :2]
        P2 = Y[sel] @ B
        area = hull_area(P2); diam = np.max(np.linalg.norm(P2 - P2.mean(0), axis=1)) * 2
        worstA = max(worstA, area / (eb * eu)); worstD = max(worstD, diam / eb)
    # P2: V1 random 2-plane, a random line in it through a random point of P(Q)
    V, _ = np.linalg.qr(rng.standard_normal((4, 2)))
    Z = Q @ V
    d = rng.standard_normal(2); d /= np.linalg.norm(d); m = np.array([-d[1], d[0]])
    c = Z[rng.integers(n)] @ m
    sel = np.abs(Z @ m - c) < tau
    if sel.sum() > 10:
        t = Z[sel] @ d
        worstL = max(worstL, (t.max() - t.min()) / eb)
print(f"P3: section area / (e_bar e_u) <= {worstA:.3f};  section diameter / e_bar <= {worstD:.3f}")
print(f"P2: section length / e_bar <= {worstL:.3f}")
