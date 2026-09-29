"""Committed T-count per share for block-batched rotations with clean ancillas.

A block of g data qubits receives the diagonal phase exp(i sum_{j in B} theta_j x_j).
  1. unary-iteration QROM loads phi(x_B) = sum_j theta_j x_j (b-bit word) : 4(2^g - 1) T
     (Babbush et al. 2018; data are written with CNOTs, so the cost does not depend on b)
  2. add the b-bit register into a catalytic phase-gradient register    : 4(b - 1) T
     (Gidney 2018 adder, temporary-AND, measurement-based uncomputation)
  3. uncompute the QROM, conservatively at full cost                    : 4(2^g - 1) T
b = ceil(log2(2 pi g / eps)) + 1 so that the g rounding errors sum to <= eps/2.
Compared with ancilla-free per-rotation synthesis at 3 log2(1/eps) (Ross-Selinger mean).
"""
import math

def per_share(eps, g):
    b = math.ceil(math.log2(2 * math.pi * g / eps)) + 1
    return (8 * (2 ** g - 1) + 4 * (b - 1)) / g, b

for eps in (1e-10, 1e-15, 1e-20, 2.0 ** -100, 2.0 ** -1000):
    L = math.log2(1 / eps)
    best = min((per_share(eps, g)[0], g) for g in range(1, 40))
    c, g = best
    b = per_share(eps, g)[1]
    print(f"L={L:7.1f}  best g={g:2d}  T/share={c:8.1f}  ancilla-free 3L={3*L:8.1f}  "
          f"ratio={c/(3*L):.3f}  T per bit={c/L:.3f}  (4/log2 L={4/math.log2(L):.3f})  qubits~{g+(g-1)+3*b}")
