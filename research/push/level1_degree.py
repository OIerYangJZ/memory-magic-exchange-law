"""Exact level-1 box condition for auxiliary degree d on S^3 (jets: transverse order j,
s-orders 0..2(d-j), multiplicity j+1; basis degrees sum_n n (n+1)^2), and the rate obtained
if level 2 on a degree-d auxiliary surface behaved like the worst round sphere
(MODEL ASSUMPTION, not proved for d >= 2)."""
import numpy as np
exec(open('rigid_plane_copy.py').read().split('grid = np.round')[0])

def l1_coeffs(d):
    A = sum(j * (j + 1) * (2 * (d - j) + 1) for j in range(d + 1))
    B = sum((j + 1) * (2 * (d - j)) * (2 * (d - j) + 1) // 2 for j in range(d + 1))
    S = sum(n * (n + 1) ** 2 for n in range(d + 1))
    return A, B, 2 * S          # condition: A a' + B b' > 2S

for d in (1, 2, 3, 5, 10, 40):
    A, B, S2 = l1_coeffs(d)
    print(f"d={d:3d}: {A} a' + {B} b' > {S2}   (normalised: a'*{A/S2*3:.3f} + b'*{B/S2*3:.3f} > 3)")

def count_exponent_d(a0, grid, d):
    A, B, S2 = l1_coeffs(d)
    best = (np.inf, None)
    for ap in grid[grid >= a0][:50]:
        bmin = max((S2 - A * ap) / B, 0)
        for bp in grid[grid >= bmin + 1e-9][:40]:
            E1 = 2 * (ap - a0) + bp
            if E1 >= best[0]:
                continue
            E2 = max(level2(ap, bp, g, grid) for g in np.linspace(0, 3.5, 71))
            if E1 + E2 < best[0]:
                best = (E1 + E2, (round(ap, 3), round(bp, 3), round(E2, 3)))
    return best

grid = np.round(np.arange(0, 6.0001, 0.02), 4)
for d in (2, 4, 40):
    out = []
    for tau in np.arange(2.2, 2.62, 0.04):
        E, arg = count_exponent_d(4 / tau, grid, d)
        out.append((round(tau, 2), round(E * tau / 4, 3)))
    ts = [t for t, b in out if b >= 1]
    print(f"d={d}: bits profile {out};  bits reach L at tau ~ {ts[0] if ts else '>2.6'}")
