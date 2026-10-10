#!/usr/bin/env python3
"""
conj1_scan.py -- Matsumoto-Amano enumeration of single-qubit Clifford+T words and a
systematic scan of the constant c of Conjecture 1 (eq. (13) of the paper) over frames:
identity, Clifford, Haar-random, exact axes of cheap words (HT, SHT, ...),
perturbed axes, and angular offsets of the grid.

For each frame (G1,G2), accuracy eps (Q = Q_eps) and cost tau it records
   N_grid(tau) = #{W : t_min(W) <= tau, exists d: dproj(W, G1 Rz(2 pi d/Q) G2) <= eps}
   N_tube(tau) = #{W : t_min(W) <= tau, exists theta: dproj(W, G1 Rz(theta) G2) <= eps}
and the constant that (11) would need with (c0,c1) = (8,2):
   c_req(tau) = (N_grid(tau) - 8 - 2 tau) / (2^tau eps^2).

Usage: python3 conj1_scan.py [TMAX]      (default TMAX = 14; 16 takes a few minutes)
"""
import sys, time, itertools
import numpy as np

TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
rng = np.random.default_rng(914)

# ---------------------------------------------------------------- gates
H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)
T = np.diag([1, np.exp(1j * np.pi / 4)]).astype(complex)
I2 = np.eye(2, dtype=complex)
HT, SHT = H @ T, S @ H @ T

def proj_key(U):
    """hashable key of U modulo phase (for dedup checks)."""
    # normalize phase: make the largest-modulus entry real positive
    f = U.flatten()
    idx = next(i for i in range(4) if abs(f[i]) > 1e-9)
    V = f / (f[idx] / abs(f[idx]))
    return tuple(np.round(V, 7).view(float))

# Clifford group mod phase by closure
cliff = {proj_key(I2): I2}
frontier = [I2]
while frontier:
    new = []
    for C in frontier:
        for g in (H, S):
            U = C @ g
            k = proj_key(U)
            if k not in cliff:
                cliff[k] = U; new.append(U)
    frontier = new
CL = np.array(list(cliff.values()))
assert len(CL) == 24, len(CL)

# ---------------------------------------------------------------- shells
def shells(tmax):
    """yield (t, array of shape (n_t,2,2)) for t = 0..tmax, MA normal form (T|e)(HT|SHT)^* C."""
    yield 0, CL.copy()
    core = np.array([I2])            # products of syllables
    prev_core = None
    for t in range(1, tmax + 1):
        prev_core = core
        core = np.concatenate([core @ HT, core @ SHT])            # 2^t
        a = np.einsum('nij,cjk->ncik', core, CL).reshape(-1, 2, 2)              # prefix e
        b = np.einsum('ij,njk,ckl->ncil', T, prev_core, CL).reshape(-1, 2, 2)   # prefix T
        yield t, np.concatenate([a, b])

t0 = time.time()
SH = dict(shells(TMAX))
print(f"enumerated shells t<=TMAX={TMAX}: sizes {[len(SH[t]) for t in range(TMAX+1)]}  ({time.time()-t0:.1f}s)")
assert all(len(SH[t]) == 36 * 2 ** t for t in range(1, TMAX + 1))
# uniqueness check at small t
for t in range(0, 7):
    keys = {proj_key(U) for U in SH[t]}
    assert len(keys) == len(SH[t]), (t, len(keys))
# unitarity
err = max(np.max(np.abs(np.einsum('nij,nkj->nik', SH[t], SH[t].conj()) - I2)) for t in SH)
print(f"MA normal form unique for t<=6; max unitarity error {err:.2e}")

# ---------------------------------------------------------------- frames
def rz(th): return np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])

def su2(U):
    return U / np.sqrt(np.linalg.det(U))

def eigframe(V):
    """D with columns = eigenvectors of V (det 1): coset D Rz D^dag = rotations about V's axis."""
    w, D = np.linalg.eig(su2(V))
    # order eigenvalues so that D Rz(theta) D^dag matches sign convention; irrelevant for counting
    D = D / np.sqrt(np.linalg.det(D))
    return D

def axis_of(V):
    V = su2(V)
    n = np.array([np.imag(V[0, 1] + V[1, 0]) / 2, np.real(V[1, 0] - V[0, 1]) / 2, np.imag(V[0, 0] - V[1, 1]) / 2])
    # V = cos(phi/2) I - i sin(phi/2) n.sigma  => -i sin n_x = V01 ... (sign conventions aside, only |n| direction matters)
    nrm = np.linalg.norm(n)
    return n / nrm if nrm > 1e-12 else None

def haar():
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    q = q * (np.diag(r) / abs(np.diag(r)))
    return su2(q)

def rot(axis, phi):
    X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1, -1]).astype(complex)
    n = axis / np.linalg.norm(axis)
    return np.cos(phi / 2) * I2 - 1j * np.sin(phi / 2) * (n[0] * X + n[1] * Y + n[2] * Z)

frames = {}
frames['identity'] = (I2, I2)
frames['Clifford pair'] = (CL[5], CL[17])
for i in range(3):
    frames[f'Haar {i}'] = (haar(), haar())
D_HT = eigframe(HT)
frames['axis HT (exact)'] = (D_HT, D_HT.conj().T)
frames['axis SHT (exact)'] = (eigframe(SHT), eigframe(SHT).conj().T)
def axis_key(n):
    """key of the rotation axis +-n, independent of floating-point noise: rounded, -0.0 -> 0.0, and the
    lexicographically larger of +n and -n (choosing the sign by argmax|n_i| is noise-dependent when two
    components have equal modulus)"""
    r = np.round(n, 6) + 0.0
    return max(tuple(r.tolist()), tuple((-r + 0.0).tolist()))

# exact axes of all words with t<=3, distinct axes (up to sign), sampled
axes = {}
for t in range(1, 4):
    for U in SH[t]:
        n = axis_of(U)
        if n is None: continue
        axes.setdefault(axis_key(n), (t, U))
ax_list = sorted(axes.items(), key=lambda kv: kv[0])   # sorted, so the sample does not depend on enumeration order
print(f"distinct rotation axes of words with t<=3: {len(ax_list)}")
sel = [ax_list[i] for i in rng.choice(len(ax_list), size=min(12, len(ax_list)), replace=False)]
for k, (t, U) in sel:
    D = eigframe(U); frames[f'axis of t={t} word {np.round(k,3)}'] = (D, D.conj().T)
# perturbed HT axis
for ang in (0.01, 0.03, 0.1):
    Dp = rot(rng.normal(size=3), ang) @ D_HT
    frames[f'axis HT perturbed {ang} rad'] = (Dp, Dp.conj().T)
# transposed HT frame (G1 = D^dag, G2 = D): the coset of the inverse-axis words
frames['axis HT (G1=D^dag,G2=D)'] = (D_HT.conj().T, D_HT)

# ---------------------------------------------------------------- counting
def Qeps(eps): return int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))

def counts(G1, G2, eps, Q, offset=0.0):
    """per-shell grid and tube hit counts."""
    thr = 2 - eps ** 2                       # |tr| >= 2 - eps^2  <=>  dproj <= eps
    th = 2 * np.pi * (np.arange(Q) + offset) / Q
    e1, e2 = np.exp(1j * th / 2), np.exp(-1j * th / 2)
    g, tb = [], []
    for t in range(TMAX + 1):
        W = SH[t]
        M = np.einsum('ij,njk,kl->nil', G1.conj().T, W, G2.conj().T)
        m0, m1 = M[:, 0, 0], M[:, 1, 1]
        tube = (np.abs(m0) + np.abs(m1)) >= thr
        best = np.zeros(len(W))
        for c in range(0, len(W), 200000):
            v = np.abs(m0[c:c+200000, None] * e1[None, :] + m1[c:c+200000, None] * e2[None, :])
            best[c:c+200000] = v.max(axis=1)
        g.append(int(np.sum(best >= thr))); tb.append(int(np.sum(tube)))
    return np.array(g), np.array(tb)

# ---- 1. reproduce the editor's instance: HT axis, eps=1/32, Q=25, tau<=11
eps = 1 / 32; Q = Qeps(eps); assert Q == 25
g, tb = counts(*frames['axis HT (exact)'], eps, Q)
print(f"\n[editor's instance] HT-axis frame, eps=1/32, Q={Q}: grid hits by exact T-count 0..{TMAX}: {g.tolist()}")
print(f"  cumulative at tau=11: {g[:12].sum()}   (editor: 137; allowance 8+2*11+48*2^11/1024 = {8+22+48*2**11/1024:.0f})")
print(f"  tube hits by shell: {tb.tolist()}; Haar tube prediction 36*2^t*eps^2: {[round(36*2**t*eps**2,1) for t in range(TMAX+1)]}")
print(f"  grid window angular coverage at Q={Q}: {Q*8*np.arcsin(eps/2)/(2*np.pi):.3f} of the coset")
# transposed frame too
g2, _ = counts(*frames['axis HT (G1=D^dag,G2=D)'], eps, Q)
print(f"  same with (G1,G2)=(D^dag,D): {g2.tolist()}, cumulative at 11: {g2[:12].sum()}")

# ---- 2. systematic scan
print("\n[scan] required c = (N_grid(<=tau) - 8 - 2 tau)/(2^tau eps^2), max over tau with 2^tau eps^2 >= 1 and tau <= TMAX")
print(f"{'frame':44s} " + " ".join(f"eps=2^-{k:<2d}" for k in (3, 4, 5, 6)) + "   | tube c_req (same cells)")
rows = []
worst = (0, None)
for name, (G1, G2) in frames.items():
    line = f"{name:44s} "; tline = ""
    for k in (3, 4, 5, 6):
        eps = 2.0 ** -k; Q = Qeps(eps)
        g, tb = counts(G1, G2, eps, Q)
        cg, ct = 0.0, 0.0; arg = None
        for tau in range(TMAX + 1):
            if 2 ** tau * eps ** 2 < 1: continue
            Ng, Nt = g[:tau+1].sum(), tb[:tau+1].sum()
            r = (Ng - 8 - 2 * tau) / (2 ** tau * eps ** 2)
            if r > cg: cg, arg = r, tau
            ct = max(ct, (Nt - 8 - 2 * tau) / (2 ** tau * eps ** 2))
        line += f"{cg:8.1f}@{arg if arg is not None else '-':<3}"; tline += f"{ct:7.1f}"
        rows.append((name, k, cg, arg, ct))
        if cg > worst[0]: worst = (cg, (name, k, arg))
    print(line + "   | " + tline)
print(f"\nmax required c over the scan: {worst[0]:.1f} at {worst[1]}  (the paper assumes c=256; Haar grid value 24; Haar tube value 72)")

# ---- 3. angular offsets on the HT-axis frame
print("\n[offsets] HT-axis frame, grid shifted by a fraction of a step (equivalent to G1 -> G1 Rz(offset)):")
for k in (4, 5, 6):
    eps = 2.0 ** -k; Q = Qeps(eps)
    out = []
    for off in (0.0, 0.25, 0.5, 0.75):
        g, _ = counts(*frames['axis HT (exact)'], eps, Q, offset=off)
        cg = max(((g[:tau+1].sum() - 8 - 2*tau) / (2**tau*eps**2)) for tau in range(TMAX+1) if 2**tau*eps**2 >= 1)
        out.append(f"off={off}: {cg:6.1f}")
    print(f"  eps=2^-{k}, Q={Q}: " + "  ".join(out))
print(f"\ntotal time {time.time()-t0:.0f}s")

# ---- 4. structure of the axis excess: V-orbits inside the tube
print("\n[orbit check] HT-axis frame, eps=1/32: fraction of tube words of shell t that equal HT*W (mod phase) for a tube word W of shell t-1")
G1, G2 = frames['axis HT (exact)']; eps = 1/32; thr = 2 - eps**2
tube_sets = {}
for t in range(TMAX + 1):
    W = SH[t]; M = np.einsum('ij,njk,kl->nil', G1.conj().T, W, G2.conj().T)
    mask = (np.abs(M[:, 0, 0]) + np.abs(M[:, 1, 1])) >= thr
    tube_sets[t] = [U for U in W[mask]]
for t in range(7, TMAX + 1):
    keys_t = {proj_key(U) for U in tube_sets[t]}
    lifted = sum(1 for U in tube_sets[t-1] if proj_key(HT @ U) in keys_t)
    print(f"  t={t:2d}: tube shell {len(tube_sets[t]):4d}, of which HT*(shell t-1 tube) = {lifted:4d}, new = {len(tube_sets[t])-lifted:4d}, Haar shell 36*2^t*eps^2 = {36*2**t*eps**2:.0f}")

# ---- 5. Theorem 8 numbers as a function of c (eps=1e-10, S=delta=0), c0=8, c1=2
print("\n[Theorem 8 parametric] eps=1e-10, S=delta=0, (c0,c1)=(8,2):")
L = np.log2(1e10); Q = int(np.floor(np.pi/(2*np.arcsin(2e-10)))); k = np.log2(Q-1)
def fixpoint(f):
    lo, hi = 0.0, 1e4
    for _ in range(200):
        mid = (lo+hi)/2
        if mid - f(mid) < 0: lo = mid
        else: hi = mid
    return hi
print(f"{'c':>6} {'tau0':>5} {'N0':>5} {'a0':>6} {'b':>7} {'kappa':>6} {'Thm5 T/coord':>13} {'A>/m (14)':>10}")
for c in (24, 48, 53.5, 96, 192):
    tau0 = int(np.floor(2*L - np.log2(c))); N0 = 8 + 2*tau0 + 1; a0 = 1 + np.log2(N0)
    b = 2*L - np.log2(c) - 1; kappa = 1 + b/k
    Tt = fixpoint(lambda T: kappa*(k - a0) - 3*np.log2(1+T))
    print(f"{c:6.1f} {tau0:5d} {N0:5d} {a0:6.2f} {b:7.2f} {kappa:6.3f} {Tt:13.2f} {(k-a0)/k:10.3f}")
