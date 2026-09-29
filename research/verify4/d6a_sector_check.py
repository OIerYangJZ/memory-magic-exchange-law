"""Numerical test of the D6(a) sector bound  #(K cap V cap O*) <= C max(1, h, beta h^2/H_V),
beta = min(1, psi/sin theta2), for a K-plane V and a 2-plane Pi at prescribed principal angles to V1."""
import numpy as np, math, itertools, random
exec(open('d6_dual_check.py').read().split("random.seed(3)")[0])
rng = np.random.default_rng(7)
def orth(A):
    q, _ = np.linalg.qr(A); return q
def enum_ball(Bas, Rad):
    """all integer combos c with |c @ Bas| <= Rad (Bas rows, full rank 4 in R^8 subspace)"""
    G = Bas @ Bas.T; Ginv = np.linalg.inv(G)
    bnd = [int(math.floor(Rad * math.sqrt(Ginv[i, i]))) for i in range(len(Bas))]
    rngs = [range(-k, k + 1) for k in bnd]
    C = np.array(list(itertools.product(*rngs)))
    X = C @ Bas
    return X[np.linalg.norm(X, axis=1) <= Rad + 1e-9]
random.seed(5)
for trial in range(6):
    v1 = sum(random.randint(-2, 2) * Ostar[i] for i in range(8))
    v2 = sum(random.randint(-2, 2) * Ostar[i] for i in range(8))
    Vp = kernel(list(Ostar), [v1, v2]); V = kernel(list(Ostar), list(Vp))
    # LLL-ish: just use V as is; H_V
    HV = covol(V) / 8
    # V1 = sigma1(V) (2-plane in R^4): orthonormal basis f1,f2
    F = orth(V[:, :4].T)[:, :2]
    comp = orth(np.hstack([F, rng.standard_normal((4, 2))]))[:, 2:]
    for (t1, t2) in [(0.0, 0.3), (0.0, 0.05), (0.01, 0.2), (0.0, 1.2)]:
        # Pi spanned by cos t_i f_i + sin t_i g_i  (g_i orthonormal in V1^perp)
        Pi = np.stack([math.cos(t1) * F[:, 0] + math.sin(t1) * comp[:, 0],
                       math.cos(t2) * F[:, 1] + math.sin(t2) * comp[:, 1]], axis=1)
        PPi = Pi @ Pi.T; PPerp = np.eye(4) - PPi
        for hval, psi in [(400.0, 0.02), (1500.0, 0.01)]:
            pts = enum_ball(V, math.sqrt(2.2 * hval))
            s1 = pts[:, :4]; s2v = pts[:, 4:]
            inK = (np.linalg.norm(s1 @ PPi, axis=1) <= math.sqrt(hval)) & \
                  (np.linalg.norm(s1 @ PPerp, axis=1) <= psi * math.sqrt(hval)) & \
                  (np.linalg.norm(s2v, axis=1) <= math.sqrt(hval))
            cnt = int(inK.sum())
            th = np.arccos(np.clip(np.linalg.svd(F.T @ Pi)[1], -1, 1))
            beta = min(1.0, psi / max(math.sin(th.max()), 1e-12))
            pred = max(1, hval, beta * hval ** 2 / HV)
            print(f"H_V={HV:9.1f} theta=({th.min():.3f},{th.max():.3f}) h={hval:6.0f} psi={psi}: count={cnt:6d}  "
                  f"max(1,h,beta h^2/H_V)={pred:10.1f}  ratio={cnt/pred:.3f}   [beta h^2/H_V={beta*hval**2/HV:.1f}]")
