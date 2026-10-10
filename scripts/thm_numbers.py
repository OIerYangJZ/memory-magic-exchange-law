"""thm_numbers.py -- the bounds of Theorems 1, 3 and 8 per coordinate and per pre-final round,
as functions of the bits shed per coordinate b = log2 K - S/m, at eps = 1e-10, delta = 0, m = 1e4.
Reproduces Table I (b = log2 K) and the curves of Fig. 3(b).
Theorem 8 is evaluated assuming Conjecture 1 with (c0, c1, c) = (8, 2, 256): c = 256 is above the largest
requirement the enumeration has found (c >= 163.9 at the identity coset, Sec. IV) but is not proven to suffice.
c = 48, used in earlier versions, is excluded by that requirement and is no longer evaluated."""
import numpy as np

EPS = 1e-10
L = np.log2(1 / EPS)
Q = int(np.floor(np.pi / (2 * np.arcsin(2 * EPS))))
K = Q - 1
k = np.log2(K)
m = 1e4

# constants C1, C2 of Theorem 2 and A_2 of Theorem 3 (eqs. (5), (6))
C1 = 72 * 4 * np.sqrt(2)
C2 = 24 / (1 - 2 ** -0.5) * 2 ** 1.5
tau_star = int(np.ceil(2 * k))
Z = sum((C1 * 2 ** t * EPS ** 2 + C2 * (t + 1) * 2 ** (t / 2)) * 2 ** (-t / 2) for t in range(tau_star + 1))
A2 = 1 + np.log2(Z)

# Theorem 8 constants under Conjecture 1 with (c0, c1, c) = (8, 2, c); Table I and Fig. 3 use c = 256
c0, c1 = 8, 2
beta = 2
def set_c(c_):
    global c, tau0, N0, a0, b_, kappa
    c = c_
    tau0 = int(np.floor(beta * L - np.log2(c)))
    N0 = c0 + c1 * tau0 + 1
    a0 = 1 + np.log2(N0)
    b_ = beta * L - np.log2(c) - 1
    kappa = 1 + b_ / k
C_ASSUMED = 256
C_REQUIRED = 163.9   # largest requirement found: 232 words at the identity coset, eps = 0.024532, tau = 11
set_c(C_ASSUMED)

def fixed_point(f, x0=1.0):
    """largest x >= 0 with x = f(x), f decreasing; bisection on g(x) = x - f(x)"""
    if f(0.0) <= 0:
        return 0.0
    lo, hi = 0.0, 1.0
    while hi - f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - f(mid) < 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def thm1(b):   # eq. (3): S + T + m log2 36 + rho_m(T) >= m log2 K, per coordinate, m -> large
    return fixed_point(lambda x: b - np.log2(36) - np.log2(np.e * (1 + x)) - (2 * np.log2(1 + m * x) + 1) / m)

def thm3(b):
    return max(0.0, 2 * (b - A2))

def thm8(b):   # eq. (14) with D_t/m = b at delta = 0
    return fixed_point(lambda x: kappa * (b - a0) - 3 * np.log2(1 + x))

if __name__ == '__main__':
    print(f'eps={EPS:g}  L={L:.3f}  Q_eps={Q}  log2K={k:.3f}  (Theorem 8 at c={c})')
    print(f'tau*={tau_star}  A2={A2:.3f}')
    print(f'tau0={tau0}  N0={N0}  log2N0={np.log2(N0):.3f}  a0={a0:.3f}  kappa={kappa:.4f}  fee=2L-{np.log2(c)+1:.2f}')
    for name, f, slope in (('Thm 1', thm1, 1), ('Thm 3', thm3, 2), ('Thm 8', thm8, kappa)):
        x = f(k)
        print(f'{name}: T_t/m = {x:.3f}   per bit = {x / k:.3f}   slope = {slope:.3f}')
    for name, f in (('Thm 1', thm1), ('Thm 3', thm3), ('Thm 8', thm8)):
        h = 1e-4
        print(f'{name}: finite-precision slope -dT/dS at S=0 = {(f(k) - f(k - h)) / h:.4f}')
    print(f'achievable (gridsynth mean at eps, 102.32): per bit = {102.32 / k:.3f}')
    print(f'achievable (gridsynth mean at eps/2, 105.34): per bit = {105.34 / k:.3f}, slope = {105.34 / np.log2(Q):.4f}')
    print(f'ratio achievable(eps/2) / Thm 8 = {105.34 / thm8(k):.2f};  / Thm 1 = {105.34 / thm1(k):.2f}')
    print(f'eq. (M1) offset: 3L - log2(4c/pi) = 3L - {np.log2(4 * c / np.pi):.2f}')
    # cross-over of Thm 1 and Thm 3 in S/m
    bs = np.linspace(0, k, 2001)
    d = np.array([thm3(b) - thm1(b) for b in bs])
    i = np.argmax(d > 0)
    print(f'Thm 3 exceeds Thm 1 for S/m < {k - bs[i]:.2f} bits')
    # kappa at 2^-80
    L80 = 80.0; K80 = np.floor(np.pi / (2 * np.arcsin(2 * 2.0 ** -80))) - 1
    print(f'kappa at eps=2^-80: {1 + (2 * L80 - np.log2(c) - 1) / np.log2(K80):.3f}')
    # sensitivity: Theorem 8 at the largest requirement found, c = 163.9
    x256 = thm8(k)
    set_c(C_REQUIRED)
    x = thm8(k)
    print(f'c={c}: tau0={tau0}  a0={a0:.3f}  kappa={kappa:.4f}  Thm 8: T_t/m = {x:.3f}   per bit = {x / k:.3f}'
          f'   ({100 * (x / x256 - 1):.1f}% above c={C_ASSUMED})')
    set_c(C_ASSUMED)
