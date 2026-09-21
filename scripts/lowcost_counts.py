"""lowcost_counts.py -- exhaustive Matsumoto-Amano enumeration to a chosen T-count,
(a) counts of words of cost <= tau0 within eps of a grid rotation on several frames,
including an axis-aligned frame in the sense of Proposition (cluster);
(b) exact optimal synthesis cost tau_min(d) on the tuned grids at L = 5, 6.

Words are (T|e)(HT|SHT)^* C, C one of the 24 Cliffords modulo phase.

Called from the package root, like the other scripts:  ../.venv/bin/python scripts/lowcost_counts.py
(~5 min).  Writes data/lowcost_counts.json (task A, Table 2) and data/taumin_exact.json (task B,
Fig. 2(a)); the latter is independent of the older data/taumin.pkl, which comes from enum2.py.
"""
import numpy as np, json, sys, time

H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def su2(M):
    """normalize batch (N,2,2) to det 1"""
    det = M[..., 0, 0] * M[..., 1, 1] - M[..., 0, 1] * M[..., 1, 0]
    return M / np.sqrt(det)[..., None, None]

def proj_key(M):
    """hashable key modulo phase"""
    M = su2(M)
    # fix phase sign: make first nonzero entry have positive real part (or positive imag if real=0)
    a = M[0, 0]
    ref = a if abs(a) > 1e-9 else M[1, 0]
    ph = ref / abs(ref)
    M = M / ph
    return tuple(np.round(M.flatten(), 6))

def cliffords():
    seen = {proj_key(I2): I2}
    frontier = [I2]
    while frontier:
        new = []
        for M in frontier:
            for g in (H, S):
                N = M @ g
                k = proj_key(N)
                if k not in seen:
                    seen[k] = N; new.append(N)
        frontier = new
    assert len(seen) == 24, len(seen)
    return su2(np.array(list(seen.values())))

C24 = cliffords()
HT = H @ T
SHT = S @ H @ T

def prefixes(t):
    """all (T|e)(HT|SHT)^* prefixes of T-count exactly t, shape (3*2^(t-1),2,2)"""
    if t == 0:
        return I2[None]
    P = np.array([T, HT, SHT])
    for _ in range(t - 1):
        P = np.concatenate([P @ HT, P @ SHT])
    return P

def words_of_level(t):
    P = prefixes(t)
    W = (P[:, None] @ C24[None]).reshape(-1, 2, 2)
    return su2(W)

def near_grid(W, G1, G2, Q, eps):
    """for words W (already SU(2)): nearest grid index d and dproj on coset G1 Rz(2 pi d/Q) G2"""
    Wp = np.conj(G1.T)[None] @ W @ np.conj(G2.T)[None]
    a = Wp[:, 0, 0]
    d0 = np.round(-np.angle(a) * Q / np.pi).astype(np.int64) % Q
    re = np.abs(np.real(np.exp(1j * np.pi * d0 / Q) * a))
    dist = np.sqrt(np.maximum(0.0, 2 - 2 * re))
    return d0, dist

def haar_su2(rng):
    z = rng.normal(size=4); z /= np.linalg.norm(z)
    a = z[0] + 1j * z[1]; b = z[2] + 1j * z[3]
    return np.array([[a, -np.conj(b)], [b, np.conj(a)]])

def axis(M):
    M = su2(M[None])[0]
    s = np.sqrt(max(1e-30, 1 - np.real(M[0, 0]) ** 2))
    n = np.array([-np.imag(M[1, 0]), np.real(M[1, 0]), -np.imag(M[0, 0])]) / s
    return n / np.linalg.norm(n)

def axis_angle(n1, n2):
    c = abs(float(np.dot(n1, n2)))
    return np.arccos(min(1.0, c))

def random_word(t, rng):
    P = prefixes(t)
    return su2((P[rng.integers(len(P))] @ C24[rng.integers(24)])[None])[0]

def tuned_Q(eps):
    return int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))

# ---------------------------------------------------------------- task A
def task_A(TMAX=16, Ls=(3, 4, 5, 6, 7, 8)):
    rng = np.random.default_rng(20260910)
    frames = {
        'identity': (I2, I2),
        'Clifford pair': (C24[7], C24[13]),
        'arith T-count 5': (random_word(5, rng), random_word(5, rng)),
        'arith T-count 10': (random_word(10, rng), random_word(10, rng)),
        'arith T-count 15': (random_word(15, rng), random_word(15, rng)),
        'Haar pair A': (haar_su2(rng), haar_su2(rng)),
        'Haar pair B': (haar_su2(rng), haar_su2(rng)),
    }
    # axis-aligned frame: word G1 of cost <= TMAX whose image of the z-axis, G1 z G1^-1,
    # is closest to the rotation axis of V = HT (so that the powers (HT)^j lie near the coset G1 Rz G1^-1)
    sig = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)
    V = su2(HT[None])[0]
    nV = np.array([np.imag(np.trace(sig[i] @ V)) for i in range(3)]); nV /= np.linalg.norm(nV)
    best = (9.0, None, None)
    for t in range(0, TMAX + 1):
        W = words_of_level(t)
        Z = W @ sig[2][None] @ np.conj(np.transpose(W, (0, 2, 1)))
        n = np.stack([np.real(np.trace(sig[i][None] @ Z, axis1=1, axis2=2)) for i in range(3)], 1) / 2
        ang = np.arccos(np.minimum(1.0, np.abs(n @ nV)))
        i = int(np.argmin(ang))
        if ang[i] < best[0]:
            best = (float(ang[i]), W[i].copy(), t)
    eta, G1, tG1 = best
    frames['axis-aligned (HT), cost %d, angle %.2e' % (tG1, eta)] = (G1, np.conj(G1.T))
    print('axis alignment: T-count', tG1, 'axis angle', eta, file=sys.stderr)

    epss = [2.0 ** (-L) for L in Ls]
    Qs = [tuned_Q(e) for e in epss]
    tau0 = [int(np.floor(2 * L - np.log2(24))) for L in Ls]
    counts = {f: np.zeros((TMAX + 1, len(Ls)), dtype=np.int64) for f in frames}
    for t in range(TMAX + 1):
        W = words_of_level(t)
        for f, (G1, G2) in frames.items():
            for k, (e, Q) in enumerate(zip(epss, Qs)):
                _, dist = near_grid(W, G1, G2, Q, e)
                counts[f][t, k] = int(np.sum(dist <= e))
        print('level', t, 'done', file=sys.stderr)
    out = {'Ls': Ls, 'Qs': Qs, 'tau0': tau0, 'eta': eta, 'tG1': tG1, 'TMAX': TMAX,
           'frames': {f: c.tolist() for f, c in counts.items()},
           'cumulative_at_tau0': {f: [int(c[:tau0[k] + 1, k].sum()) for k in range(len(Ls))] for f, c in counts.items()},
           'reference_c0_c1tau0': [8 + 2 * x for x in tau0]}
    return out

# ---------------------------------------------------------------- task B
def task_B(L, TMAX):
    eps = 2.0 ** (-L); Q = tuned_Q(eps)
    taumin = np.full(Q, -1, dtype=np.int64); taumin[0] = 0
    for t in range(0, TMAX + 1):
        P = prefixes(t)
        for c in range(24):
            W = su2(P @ C24[c])
            d0, dist = near_grid(W, I2, I2, Q, eps)
            hit = d0[dist <= eps]
            for d in np.unique(hit):
                if taumin[d] < 0:
                    taumin[d] = t
        if np.all(taumin >= 0):
            break
    return Q, taumin

if __name__ == '__main__':
    t0 = time.time()
    if '--B-only' not in sys.argv:
        A = task_A()
        json.dump(A, open('data/lowcost_counts.json', 'w'), indent=1)
        print('task A done', time.time() - t0, file=sys.stderr)
    res = {}
    for L, TMAX in (() if '--A-only' in sys.argv else ((5, 16), (6, 18))):
        Q, tm = task_B(L, TMAX)
        res[L] = {'Q': Q, 'taumin': tm.tolist()}
        print('L', L, 'Q', Q, 'unreached', int(np.sum(tm < 0)), 'mean over Z_Q', float(tm.mean()),
              'mean d!=0', float(tm[1:].mean()), 'max', int(tm.max()), time.time() - t0, file=sys.stderr)
    if res:
        json.dump(res, open('data/taumin_exact.json', 'w'))
