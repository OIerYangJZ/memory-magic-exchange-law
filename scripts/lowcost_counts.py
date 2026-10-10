"""lowcost_counts.py -- exhaustive Matsumoto-Amano enumeration to a chosen T-count,
(a) counts of words of cost <= tau0 within eps of a grid rotation on several frames,
including an axis-aligned frame in the sense of Proposition (cluster);
(b) exact optimal synthesis cost tau_min(d) on the tuned grids at L = 5, 6.

Words are (T|e)(HT|SHT)^* C, C one of the 24 Cliffords modulo phase.

Run from the repository root:  .venv/bin/python scripts/lowcost_counts.py   (about a minute).
Writes data/lowcost_counts.json (task A, Table S2 of the Supplemental Material) and
data/taumin_exact.json (task B, Fig. 3(a)); the latter is independent of the older data/taumin.pkl,
which comes from enum2.py.

Ties.  On arithmetic frames a word can lie at distance exactly eps from a grid rotation, and then
floating point alone decides the count.  Every word whose float distance is within TIE of eps is
therefore rebuilt from its gate sequence and re-evaluated with mpmath at 50 digits; a distance that
agrees with eps to MP_TIE is an exact tie and counts as within eps, as "d_proj <= eps" requires.
This makes the output independent of the numpy/BLAS version.
"""
import numpy as np, json, sys, time
import mpmath as mp

mp.mp.dps = 50
TIE = 1e-9                   # float distances this close to eps are re-evaluated in high precision
MP_TIE = mp.mpf('1e-35')     # high-precision distances this close to eps are exact ties

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
    """the 24 Cliffords modulo phase, with the H/S sequence that produced each"""
    seen = {proj_key(I2): (I2, '')}
    frontier = [(I2, '')]
    while frontier:
        new = []
        for M, q in frontier:
            for g, name in ((H, 'H'), (S, 'S')):
                N = M @ g
                k = proj_key(N)
                if k not in seen:
                    seen[k] = (N, q + name); new.append((N, q + name))
        frontier = new
    assert len(seen) == 24, len(seen)
    return su2(np.array([v[0] for v in seen.values()])), [v[1] for v in seen.values()]

C24, C24_SEQ = cliffords()
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

def prefix_seq(t, p):
    """gate sequence of prefixes(t)[p]: undo the concatenation order of prefixes()"""
    if t == 0:
        return ''
    length, tail = 3 * 2 ** (t - 1), []
    for _ in range(t - 1):
        half = length // 2
        if p < half:
            tail.append('HT')
        else:
            tail.append('SHT'); p -= half
        length = half
    return ['T', 'HT', 'SHT'][p] + ''.join(reversed(tail))

def word_seq(t, i):
    """gate sequence of words_of_level(t)[i]"""
    return prefix_seq(t, i // 24) + C24_SEQ[i % 24]

_MPG = {'H': mp.matrix([[1, 1], [1, -1]]) / mp.sqrt(2), 'S': mp.matrix([[1, 0], [0, 1j]]),
        'T': mp.matrix([[1, 0], [0, mp.expjpi(mp.mpf(1) / 4)]])}

def mp_su2(M):
    return M / mp.sqrt(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0])

def mp_word(seq):
    M = mp.eye(2)
    for g in seq:
        M = M * _MPG[g]
    return mp_su2(M)

def mp_float(M):
    """an SU(2) frame given only in floating point (Haar frames): its entries are exact by definition"""
    return mp.matrix([[mp.mpc(complex(M[i, j])) for j in range(2)] for i in range(2)])

def mp_near_grid(Wm, G1m, G2m, Q):
    """high-precision version of near_grid for one word"""
    a = (G1m.H * Wm * G2m.H)[0, 0]
    d0 = int(mp.nint(-mp.arg(a) * Q / mp.pi)) % Q
    re = abs(mp.re(mp.expjpi(mp.mpf(d0) / Q) * a))
    return d0, mp.sqrt(max(mp.mpf(0), 2 - 2 * re))

def count_within(W, t, G1, G2, G1m, G2m, Q, eps):
    """grid indices d of the words W = words_of_level(t) (or a slice of it, see task_B) within eps,
    with the near-tie words resolved in high precision; returns (indices d, number of near-ties)"""
    d0, dist = near_grid(W, G1, G2, Q, eps)
    hits = list(d0[dist < eps - TIE])
    border = np.nonzero(np.abs(dist - eps) <= TIE)[0]
    for i in border:
        d, dm = mp_near_grid(mp_word(t(i)) if callable(t) else mp_word(word_seq(t, i)), G1m, G2m, Q)
        if dm <= eps + MP_TIE:
            hits.append(d)
    return np.array(hits, dtype=np.int64), len(border)

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
    """a random word of T-count t and its gate sequence (same draws from rng as before)"""
    P = prefixes(t)
    p = rng.integers(len(P)); c = rng.integers(24)
    return su2((P[p] @ C24[c])[None])[0], prefix_seq(t, int(p)) + C24_SEQ[int(c)]

def tuned_Q(eps):
    return int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))

# ---------------------------------------------------------------- task A
def task_A(TMAX=16, Ls=(3, 4, 5, 6, 7, 8)):
    rng = np.random.default_rng(20260910)
    frames, mpframes = {}, {}
    frames['identity'] = (I2, I2); mpframes['identity'] = (mp.eye(2), mp.eye(2))
    frames['Clifford pair'] = (C24[7], C24[13])
    mpframes['Clifford pair'] = (mp_word(C24_SEQ[7]), mp_word(C24_SEQ[13]))
    for tt in (5, 10, 15):
        (g1, q1), (g2, q2) = random_word(tt, rng), random_word(tt, rng)
        frames['arith T-count %d' % tt] = (g1, g2); mpframes['arith T-count %d' % tt] = (mp_word(q1), mp_word(q2))
    for nm in ('Haar pair A', 'Haar pair B'):
        g1, g2 = haar_su2(rng), haar_su2(rng)
        frames[nm] = (g1, g2); mpframes[nm] = (mp_float(g1), mp_float(g2))
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
            best = (float(ang[i]), W[i].copy(), t, i)
    eta, G1, tG1, iG1 = best
    name = 'axis-aligned (HT), cost %d, angle %.2e' % (tG1, eta)
    frames[name] = (G1, np.conj(G1.T))
    G1m = mp_word(word_seq(tG1, iG1)); mpframes[name] = (G1m, G1m.H)
    print('axis alignment: T-count', tG1, 'axis angle', eta, file=sys.stderr)

    epss = [2.0 ** (-L) for L in Ls]
    Qs = [tuned_Q(e) for e in epss]
    tau0 = [int(np.floor(2 * L - np.log2(24))) for L in Ls]
    counts = {f: np.zeros((TMAX + 1, len(Ls)), dtype=np.int64) for f in frames}
    ties = 0
    for t in range(TMAX + 1):
        W = words_of_level(t)
        for f, (G1, G2) in frames.items():
            for k, (e, Q) in enumerate(zip(epss, Qs)):
                hits, nb = count_within(W, t, G1, G2, *mpframes[f], Q, e)
                counts[f][t, k] = len(hits); ties += nb
        print('level', t, 'done', file=sys.stderr)
    print('near-tie words re-evaluated in high precision:', ties, file=sys.stderr)
    out = {'Ls': Ls, 'Qs': Qs, 'tau0': tau0, 'eta': round(eta, 12), 'tG1': tG1, 'TMAX': TMAX,
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
            hit, _ = count_within(W, lambda p, c=c, t=t: prefix_seq(t, p) + C24_SEQ[c],
                                  I2, I2, mp.eye(2), mp.eye(2), Q, eps)
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
