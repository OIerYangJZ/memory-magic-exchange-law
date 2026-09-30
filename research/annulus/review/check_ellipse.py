"""Referee check of A1 steps 2-3 (float): for random K-planes V (sigma1-image V1) and random planes Pi,
 (i)  P_{V1perp}|_Pi has singular values sin(theta1) <= sin(theta2) (principal angles of V1 vs Pi);
 (ii) every tube point projects into the 3*delta neighbourhood of the ellipse E = P_{V1perp}(R' C);
 (iii) the arc covering: arcs of length <= l0 and turning <= Theta0 (l0*Theta0 = delta1) lie within
      delta1 of their chords; N <= C(1 + (a2/delta1)^(1/2)) and sum of rectangle lengths <= C(a2+delta1);
 (iv) the 'thickened ellipse' area is ~ a2*delta1 (vs a2^2 for the filled box of F5(b)).
We also generate planes V1 *close* to Pi (theta small), the relevant regime."""
import numpy as np, sys, random
from klat import rand_vec, kperp, sig1
import mpmath as mp
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
prng = random.Random(5)

def orth(M):
    q, _ = np.linalg.qr(M); return q
def principal(A, B):
    """principal angles between column spaces of orthonormal A, B"""
    s = np.linalg.svd(A.T @ B, compute_uv=False)
    return np.arccos(np.clip(s, -1, 1))          # ascending angles

def ellipse_dist(y, P, Pi, R):
    """distance from y to the curve s -> P (R C(s)), dense sampling + local refine"""
    s = np.linspace(0, 2 * np.pi, 4001)
    C = R * (np.outer(np.cos(s), Pi[:, 0]) + np.outer(np.sin(s), Pi[:, 1])) @ P.T
    d = np.linalg.norm(C - y, axis=1); k = np.argmin(d)
    s2 = np.linspace(s[k] - 0.002, s[k] + 0.002, 401)
    C2 = R * (np.outer(np.cos(s2), Pi[:, 0]) + np.outer(np.sin(s2), Pi[:, 1])) @ P.T
    return np.min(np.linalg.norm(C2 - y, axis=1))

worst_ratio = 0; worst_sv = 0; covstats = []
for trial in range(40):
    # V1 = sigma1 of a random K-plane; half the trials: Pi chosen near V1 (small principal angles)
    v1, v2 = rand_vec(prng, 3), rand_vec(prng, 3)
    V1 = orth(np.array([[float(c) for c in sig1(v1)], [float(c) for c in sig1(v2)]]).T)
    if trial % 2 == 0:
        Pi = orth(rng.standard_normal((4, 2)))
    else:
        tilt = 10 ** rng.uniform(-3, -0.5)
        Pi = orth(V1 + tilt * rng.standard_normal((4, 2)))
    th = principal(Pi, V1)                        # th[0] <= th[1]
    Pperp = np.eye(4) - V1 @ V1.T
    sv = np.linalg.svd(Pperp @ Pi, compute_uv=False)[::-1]      # ascending
    worst_sv = max(worst_sv, np.max(np.abs(sv - np.sin(th))))
    R = 1.0; eps = 10 ** rng.uniform(-4, -2); delta = eps * R; d1 = 3 * delta
    Piperp = orth(np.linalg.svd(np.eye(4) - Pi @ Pi.T)[0][:, :2])
    for _ in range(60):
        s = rng.uniform(0, 2 * np.pi); phi = rng.uniform(0, 2 * np.arcsin(eps / 2))
        n = Piperp @ rng.standard_normal(2); n /= np.linalg.norm(n)
        x = R * (np.cos(phi) * (np.cos(s) * Pi[:, 0] + np.sin(s) * Pi[:, 1]) + np.sin(phi) * n)
        d = ellipse_dist(Pperp @ x, Pperp, Pi, R)
        worst_ratio = max(worst_ratio, d / delta)
    # (iii) covering of the ellipse with semi-axes A1 = R sin th1 <= a2 = R sin th2
    A1, a2 = R * np.sin(th[0]), R * np.sin(th[1])
    if a2 <= d1 or A1 <= d1:
        covstats.append((a2 / d1, 1, (2 * a2 + 2 * d1) / (a2 + d1), 0.0)); continue
    t = np.linspace(0, 2 * np.pi, 200001)
    pts = np.stack([A1 * np.cos(t), a2 * np.sin(t)], 1)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1); arc = np.concatenate([[0], np.cumsum(seg)])
    tang = np.unwrap(np.arctan2(np.gradient(pts[:, 1]), np.gradient(pts[:, 0])))
    perim = arc[-1]; l0 = np.sqrt(d1 * perim / (2 * np.pi)); Th0 = d1 / l0
    cuts = [0]; i0 = 0
    for i in range(1, len(t)):
        if arc[i] - arc[i0] > l0 or abs(tang[i] - tang[i0]) > Th0:
            cuts.append(i - 1); i0 = i - 1
    cuts.append(len(t) - 1)
    maxdev = 0; Lsum = 0
    for a, b in zip(cuts[:-1], cuts[1:]):
        P0, P1 = pts[a], pts[b]; ch = P1 - P0; L = np.linalg.norm(ch)
        if L == 0: continue
        nrm = np.array([-ch[1], ch[0]]) / L
        dev = np.max(np.abs((pts[a:b + 1] - P0) @ nrm)); maxdev = max(maxdev, dev / d1)
        Lsum += L + 4 * d1
    N = len(cuts) - 1
    covstats.append((a2 / d1, N / (1 + np.sqrt(a2 / d1)), Lsum / (a2 + d1), maxdev))
cs = np.array(covstats)
print(f"(i)  max |singular values of P_V1perp|Pi - sin(principal angles)| = {worst_sv:.2e}")
print(f"(ii) max dist(P_V1perp x, E)/delta over tube points = {worst_ratio:.4f}   (A1 uses delta1 = 3 delta)")
print(f"(iii) a2/delta1 range [{cs[:,0].min():.2f}, {cs[:,0].max():.1f}];  N/(1+sqrt(a2/delta1)) <= {cs[:,1].max():.3f};"
      f"  sum l'_k/(a2+delta1) <= {cs[:,2].max():.3f};  max arc-to-chord deviation/delta1 = {cs[:,3].max():.4f} (must be <= 1)")
