#!/usr/bin/env python3
"""
mixing_bounds.py -- checks for Sec. 10.1 ("The law under mixing") of M2.

Part A: randomized checks of the three elementary inequalities used in
        Lemmas 10-11 (fidelity concentration, window decoding, tube radius).
Part B: evaluation of the constants of Theorem 8 (rate one and rate two under
        mixing) at several accuracies, S = 0, delta = 0, and the accuracy at
        which each bound stops being vacuous.

No external data; runs in a few seconds with numpy only.
"""
import numpy as np

rng = np.random.default_rng(20260912)

# ---------------------------------------------------------------- helpers
def haar_su2(n):
    """n Haar-random SU(2) matrices."""
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    ph = d / np.abs(d)
    q = q * ph[:, None, :]
    det = np.linalg.det(q)
    q = q / np.sqrt(det)[:, None, None]
    return q

def rz(theta):
    return np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)])

def rot(axis, phi):
    """SU(2) rotation by phi about unit vector axis."""
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    Z = np.array([[1, 0], [0, -1]], complex)
    n = axis / np.linalg.norm(axis)
    return np.cos(phi / 2) * np.eye(2) - 1j * np.sin(phi / 2) * (n[0] * X + n[1] * Y + n[2] * Z)

def fid(U, V):
    """process fidelity |tr(V^dag U)/2|^2 of two unitaries."""
    return abs(np.trace(V.conj().T @ U) / 2) ** 2

def dproj(U, V):
    return np.sqrt(max(0.0, 2 - abs(np.trace(V.conj().T @ U))))

def choi(kraus_list, probs):
    """Choi matrix (normalized, trace 1) of the mixed-unitary channel."""
    phi = np.zeros(4, complex)
    phi[0] = phi[3] = 1 / np.sqrt(2)
    J = np.zeros((4, 4), complex)
    for p, U in zip(probs, kraus_list):
        v = np.kron(U, np.eye(2)) @ phi
        J += p * np.outer(v, v.conj())
    return J

# ---------------------------------------------------------------- Part A
print("=== Part A: randomized checks of the elementary inequalities ===")

# A1. Window inequality (Lemma 11):
#     F(C, V R_psi) <= cos^2(psi/2) + sum p s_i + |sin psi| sum p sqrt(s_i),
#     s_i = sin^2(phi_i/2) = 1 - F(X_i, V).
worst = -np.inf
for trial in range(20000):
    k = rng.integers(1, 6)
    V = haar_su2(1)[0]
    # words at controlled small angles from V (so the bound is exercised
    # near equality as well as far from it)
    Xs = []
    for _ in range(k):
        ax = rng.normal(size=3)
        phi = abs(rng.normal()) * rng.choice([0.01, 0.1, 0.5, 2.0])
        phi = min(phi, np.pi)
        Xs.append(V @ rot(ax, phi))
    p = rng.dirichlet(np.ones(k))
    psi = rng.uniform(-np.pi, np.pi)
    Vp = V @ rz(psi)
    lhs = sum(pi * fid(X, Vp) for pi, X in zip(p, Xs))
    s = np.array([1 - fid(X, V) for X in Xs])
    rhs = np.cos(psi / 2) ** 2 + p @ s + abs(np.sin(psi)) * (p @ np.sqrt(s))
    worst = max(worst, lhs - rhs)
print(f"A1 window inequality: max(lhs - rhs) over 20000 random mixtures = {worst:.3e}  (must be <= 0)")

# A2. Tube radius: |tr(V^dag X)/2|^2 >= 1 - eta^2  =>  dproj(X, V) <= sqrt(2) eta
worst = -np.inf
for trial in range(20000):
    V = haar_su2(1)[0]
    ax = rng.normal(size=3)
    phi = rng.uniform(0, np.pi)
    X = V @ rot(ax, phi)
    eta = np.sqrt(max(1e-300, 1 - fid(X, V)))
    worst = max(worst, dproj(X, V) / (np.sqrt(2) * eta))
print(f"A2 tube radius: max dproj/(sqrt2*eta) = {worst:.6f}  (must be <= 1)")

# A3. Fidelity vs Choi trace distance: 2(1-F_e) <= || J(E) - Phi_V ||_1
worst = -np.inf
for trial in range(5000):
    k = rng.integers(1, 6)
    V = haar_su2(1)[0]
    Xs = [V @ rot(rng.normal(size=3), abs(rng.normal()) * 0.3) for _ in range(k)]
    p = rng.dirichlet(np.ones(k))
    JE = choi(Xs, p)
    JV = choi([V], [1.0])
    Fe = np.real(np.trace(JE @ JV))
    tr1 = np.sum(np.abs(np.linalg.eigvalsh(JE - JV)))
    worst = max(worst, 2 * (1 - Fe) - tr1)
print(f"A3 2(1-F_e) - ||J(E)-Phi_V||_1 : max = {worst:.3e}  (must be <= 0)")

# A4. Markov/concentration check on a synthetic mixing pipeline:
#     mixture of two words at unitary accuracy ~sqrt(eps) whose channel is
#     eps-close; verify that the qualifying mass >= 1 - eps/eta^2.
print("A4 concentration on Campbell-type two-word mixtures:")
for eps in [1e-4, 1e-6, 1e-8]:
    viol = 0
    for trial in range(2000):
        V = haar_su2(1)[0]
        a = np.sqrt(eps)
        ax = rng.normal(size=3)
        X1 = V @ rot(ax, a); X2 = V @ rot(ax, -a)      # opposite-side errors
        p = np.array([0.5, 0.5])
        JE = choi([X1, X2], p); JV = choi([V], [1.0])
        Fe = np.real(np.trace(JE @ JV))
        eps_ch = 1 - Fe                                  # <= eps by construction
        eta2 = 4 * eps_ch                                # gamma = 1/4
        qual = sum(pi for pi, X in zip(p, [X1, X2]) if 1 - fid(X, V) <= eta2)
        if qual < 1 - eps_ch / eta2 - 1e-12:
            viol += 1
    print(f"   eps={eps:.0e}: violations of the Markov mass bound = {viol}/2000")


def fixpoint(f, hi=1e6):
    """smallest T >= 0 with T >= f(T), for f decreasing (bisection; a plain
    iteration can oscillate when |f'| > 1)."""
    if f(0.0) <= 0: return 0.0
    lo, h = 0.0, 1.0
    while h < hi and h - f(h) < 0: h *= 2
    for _ in range(200):
        mid = (lo + h) / 2
        if mid - f(mid) < 0: lo = mid
        else: h = mid
    return h

# ---------------------------------------------------------------- Part B
print("\n=== Part B: constants of Theorem 8 at S = 0, delta = 0 ===")
C1, C2 = 72 * 4 * np.sqrt(2), 24 / (1 - 2 ** -0.5) * 2 ** 1.5
log2 = np.log2

RADIUS_FREE = 36 * np.sqrt(2)     # tube radius factor of the frame-free Theorem 8(b) (Lemma 15)
RADIUS_SUPP = np.sqrt(2)          # radius factor when prefix/suffix supports are bounded (Remark)

def numbers(eps, gamma=0.25, nf=0.0, det_rate_two=True, radius=RADIUS_FREE):
    L = log2(1 / eps)
    Q = float(np.floor(np.pi / (2 * np.arcsin(2 * eps))))
    K = Q - 1
    k = log2(K)
    eta = np.sqrt(eps / gamma)
    # rate one: window from Lemma 11
    w = float(np.floor((Q / np.pi) * np.arcsin(min(1.0, (1 + np.sqrt(3)) * eta))))
    win1 = log2(2 * w + 1)
    bits1 = k - win1
    log36 = log2(36)
    # rho_m(T)/m at the solution: solve T' = bits1 - log36 - rho(T')  (per coord, m->inf)
    def rho_per(T):  # per-coordinate, large-m limit: log2(e(1+T)) + o(1)
        return log2(np.e * (1 + T))
    Tp = fixpoint(lambda T: bits1 - log36 - rho_per(T))
    thm8a = (1 - gamma) * Tp          # T_t/m >= (1-gamma) T'
    # rate two: tube radius sqrt(2) eta, window from nearest-grid decoding
    eps_tube = radius * eta
    w2 = float(np.floor((2 * Q / np.pi) * np.arcsin(min(1.0, eps_tube))))
    win2 = log2(2 * w2 + 1)
    bits2 = k - win2 - 1   # one bit for the unbalanced cells of the coarse partition (cells of size h or h+1)
    tau_star = int(np.ceil(2 * bits2)) if bits2 > 0 else 0
    Z = sum((C1 * 2 ** t * eps_tube ** 2 + C2 * (t + 1) * 2 ** (t / 2)) * 2 ** (-t / 2)
            for t in range(tau_star + 1))
    A2p = 1 + log2(Z) if Z > 0 else 0.0
    # one extra bit per coordinate: the decoding window can straddle two cells
    # of the coarse partition used in Lemma 8 (Appendix F)
    thm8b = (1 - gamma) * 2 * (bits2 - 1 - A2p - nf)
    return dict(L=L, Q=Q, k=k, eta=eta, win1=win1, bits1=bits1, thm8a=thm8a,
                win2=win2, bits2=bits2, A2p=A2p, thm8b=thm8b)

hdr = f"{'eps':>10} {'L':>7} {'log2K':>7} {'coarse bits(a)':>14} {'Thm8a T/m':>10} | {'bits(b) FREE':>12} {'A2p':>6} {'Thm8b FREE':>10} | {'bits(b) SUPP':>12} {'Thm8b SUPP':>10}"
print(hdr)
for eps in [1e-10, 2.0 ** -45, 2.0 ** -60, 2.0 ** -80, 2.0 ** -120, 2.0 ** -200]:
    n = numbers(eps); m = numbers(eps, radius=RADIUS_SUPP)
    print(f"{eps:10.2e} {n['L']:7.1f} {n['k']:7.2f} {n['bits1']:14.2f} {n['thm8a']:10.2f} | {n['bits2']:12.2f} {n['A2p']:6.2f} {n['thm8b']:10.2f} | {m['bits2']:12.2f} {m['thm8b']:10.2f}")

# onset accuracies
for name, key, rad in [("rate one (Thm 8a)", "thm8a", RADIUS_FREE), ("rate two (Thm 8b) FREE", "thm8b", RADIUS_FREE), ("rate two (Thm 8b) SUPP", "thm8b", RADIUS_SUPP)]:
    Ls = np.arange(20, 600)
    vals = [numbers(2.0 ** -float(L), radius=rad)[key] for L in Ls]
    onset = next((L for L, v in zip(Ls, vals) if v > 0), None)
    print(f"{name}: first positive at L = {onset}")

# comparison with deterministic CW bounds at the same accuracies (paper's Thm 1 / Thm 3 numbers)
print("\nDeterministic CW reference at eps=1e-10 (paper): Thm1 21.75, Thm3 25.72, Thm5 52.9, achievable 105.34 (T per coordinate per round)")
n = numbers(1e-10)
print(f"Mixed model at eps=1e-10, gamma=1/4: coarse bits {n['bits1']:.2f}, Thm8a {n['thm8a']:.2f}, Thm8b {n['thm8b']:.2f} (vacuous)")

# asymptotic ratio check: Thm8b / L -> (1-gamma)
for L in [200, 400, 800]:
    n = numbers(2.0 ** -L)
    print(f"L={L}: Thm8b/L = {n['thm8b']/L:.3f}  (-> 1-gamma = {1-0.25})")

# gamma optimisation at moderate L
print("\nBest gamma for Thm 8b at L=120:")
for g in [0.05, 0.1, 0.25, 0.5]:
    n = numbers(2.0 ** -120, gamma=g)
    print(f"  gamma={g}: coarse bits {n['bits2']:.2f}, A2'={n['A2p']:.2f}, bound {n['thm8b']:.2f}")

# ---------------------------------------------------------------- Part C
print("\n=== Part C: Corollary 6 (rate three and quantization under mixing, given Conjecture 1) ===")
c0, c1, c = 8.0, 2.0, 48.0          # constants of Conjecture 1 as used in Theorem 5
c0p, c1p, cp = 9 * c0, 9 * c1, 36 * c   # after the nine-grid cover of Lemma 13

def cor6(eps, gamma=0.25, nf=0.0, S=0.0):
    n = numbers(eps, gamma=gamma, nf=nf)
    L = n['L']; Q = n['Q']; K = n['k']
    eta = np.sqrt(eps / gamma)
    eps_t = RADIUS_FREE * eta                 # tube radius of the frame-free Corollary 6
    assert eps_t <= 2 ** -6, "Corollary 6 needs eps' <= 2^-6"
    Lp = log2(1 / eps_t)
    tau0 = int(np.floor(2 * Lp - log2(cp)))
    N0 = c0p + c1p * tau0 + 1
    a0 = 1 + log2(N0)
    b = 2 * Lp - log2(cp) - 1
    w2 = float(np.floor((2 * Q / np.pi) * np.arcsin(eps_t)))
    Kp = np.ceil((2 ** K) / (2 * w2 + 1)); kbar = log2(Kp)
    kappa = 1 + b / kbar
    D = (n['bits2'] - 1) - nf - S            # per coordinate, delta = 0
    # (13'): T >= (1-gamma)[kappa (D - a0) - 3 log2(1 + T/((1-gamma)))]  (per coordinate)
    T = fixpoint(lambda T: (1 - gamma) * (kappa * (D - a0) - 3 * log2(1 + T / (1 - gamma))))
    frac = max(0.0, (D - a0) / kbar)          # (14'): fraction of coordinates above tau0
    return dict(L=L, Lp=Lp, tau0=tau0, N0=N0, a0=a0, b=b, kbar=kbar, kappa=kappa, D=D, T=T, frac=frac)

print("Theorem 8(b) numbers: FREE = frame-free version (radius 36*sqrt(2)*eta, no support assumption); SUPP = bounded-support version (radius sqrt(2)*eta, n_f=0)")
print(f"{'L':>6} {'L_tube':>7} {'tau0p':>6} {'a0p':>6} {'kappa':>6} {'D/m':>7} {'(13p) T/m':>10} {'(14p) frac>tau0':>16}")
for L in [33.219, 45, 60, 80, 120, 200]:
    r = cor6(2.0 ** -L)
    print(f"{r['L']:6.1f} {r['Lp']:7.2f} {r['tau0']:6d} {r['a0']:6.2f} {r['kappa']:6.2f} {r['D']:7.2f} {r['T']:10.2f} {r['frac']:16.3f}")
Ls = np.arange(30, 400)
onset = next((L for L in Ls if cor6(2.0 ** -float(L))['T'] > 0), None)
print(f"(13') first positive at L = {onset}; at eps=1e-10: tau0' = {cor6(1e-10)['tau0']}, "
      f"fraction of coordinates that must commit a representative of cost > tau0' >= {cor6(1e-10)['frac']:.3f}")
for L in [400, 800, 1000]:
    r = cor6(2.0 ** -L)
    print(f"L={L}: kappa'={r['kappa']:.3f}, (13')/L = {r['T']/L:.3f}  (-> 1.5(1-gamma) = {1.5*(1-0.25)})")

# ---------------------------------------------------------------- Part D
print("\n=== Part D: Proposition 4 (measurement-adaptive protocols on n qubits, rate 1/(2n+2)) ===")
Cn = {1: 24.0, 2: 11520.0}   # |C_n| modulo phase
def prop4(eps, n=2, gamma=0.25):
    nn = numbers(eps, gamma=gamma)
    bits1 = nn['bits1']
    def rho_per(N): return log2(np.e * (1 + N))
    N = fixpoint(lambda N: (bits1 - log2(Cn[n]) - rho_per(N)) / (2 * n + 2))
    return (1 - gamma) * N                      # (T + M)/m per round
for L in [33.219, 60, 120, 200, 400]:
    print(f"L={L:6.1f}: n=2 (one ancilla)  (T+M)/m >= {prop4(2.0**-L):7.2f};   n=1 crude count  T/m >= {prop4(2.0**-L, n=1):7.2f}")
onset = next((L for L in np.arange(20, 400) if prop4(2.0 ** -float(L)) > 0), None)
print(f"n=2: first positive at L = {onset}; asymptotic slope (1-gamma)/(2n+2)/2 = {(0.75/6/2):.4f} per L (floor L/16 per share at gamma=1/4)")

# ---------------------------------------------------------------- Part E
print("\n=== Part E: Lemma 13 (nine shifted grids of modulus Q_{2 eps} cover the tube of radius eps) ===")
worst = 0.0
for trial in range(20000):
    eps = rng.choice([1 / 16, 1 / 32, 1 / 100, 1e-3, 1e-5])
    Q1 = float(np.floor(np.pi / (2 * np.arcsin(4 * eps))))          # largest Q with sin(pi/2Q) >= 4 eps
    theta = rng.uniform(0, 2 * np.pi)
    grid = 2 * np.pi * (np.arange(9 * Q1)) / (9 * Q1)                # nine shifted grids of modulus Q1
    d = np.min(np.abs((theta - grid + np.pi) % (2 * np.pi) - np.pi))
    worst = max(worst, 2 * np.sin(d / 4) / eps)                       # dproj(Rz(theta), nearest) / eps
print(f"max dproj(nearest shifted-grid rotation)/eps = {worst:.4f}  (must be <= 1, so tube radius eps -> grid radius 2 eps)")

# ---------------------------------------------------------------- Part F
print("\n=== Part F: achievable numbers used in Table 4 (from [18] as quoted in [21, Eqs. (9),(13)]) ===")
print("Proposition 3 with Campbell Thm 2 (full diamond <= 10 eta^2): word accuracy eta = sqrt(eps/(5 r)):")
for eps_, r_ in [(1e-10, 2), (1e-10, 4)]:
    eta_ = np.sqrt(eps_ / (5 * r_)); print(f"  eps={eps_:.0e}, r={r_}: eta = {eta_:.3e}, RS estimate 3 log2(1/eta) + 2.7 = {3*log2(1/eta_)+2.7:.2f} T per share  (1.5 L + 1.5 log2(5r) = {1.5*log2(1/eps_)+1.5*log2(5*r_):.2f})")
for delta in [1e-10, 5e-11]:
    print(f"delta={delta:.0e}: mixed-diagonal 1.52 log2(1/delta) - 0.01 = {1.52*log2(1/delta)-0.01:6.2f};  "
          f"mixed-fallback 0.53 log2(1/delta) + 4.86 = {0.53*log2(1/delta)+4.86:6.2f};  "
          f"two-word RS mixing estimate 3 log2(1/sqrt(delta/2)) + 2.7 = {3*log2(1/np.sqrt(delta/2))+2.7:6.2f}")
# Corollary 5 boundary from [21, Eq. (10)] with their convention R_Z(theta) = e^{i theta Z} (theta = half our angle)
def T_small(theta, delta):
    Kc = (2 * np.sqrt(2 * np.e ** 3) / 3) ** (2 / 3)
    al = delta / (2 * theta) + theta
    ph = max(al - al / np.log(Kc / al), theta)
    return 3 * theta / (al + 2 * ph) * log2(12 / ((al - ph) ** 2 * (al + 2 * ph)))
eps = 1e-10; delta = eps
print("Corollary 5(i): expected T of [21, Eq.(10)] for a share whose low ell bits are shed (angle 8 eps 2^ell, delta = eps):")
for ell in [10, 12, 13, 14, 15, 16]:
    ang = 8 * eps * 2 ** ell            # our rotation angle
    th = ang / 2                        # [21] convention
    print(f"   ell={ell:2d}: angle={ang:.2e}, T_small-angle = {T_small(th, delta):6.2f}   (L/2 = {log2(1/eps)/2:.1f})")
