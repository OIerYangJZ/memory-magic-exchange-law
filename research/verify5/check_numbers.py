"""Referee spot checks for App. F (thm:dioph): numerical constants and the F5(c) minimisation."""
from fractions import Fraction as F
import numpy as np, math
a0 = a = F(28, 17); b = F(157, 100); T = F(11183, 6800)
print("2a+3b =", float(2*a+3*b), " a0+a-2b =", float(a0+a-2*b), " 2b-a =", float(2*b-a), " a-b =", float(a-b))
kD = T/4
print("kappa_D =", kD, float(kD), " 7/17-kD =", float(F(7,17)-kD), " a0-T =", float(a0-T))
print("worst 67097/40800 =", float(F(67097,40800)), " T-worst =", float(T-F(67097,40800)))
print("(T-4)/2 =", float((T-4)/2), ">= -a0?", (T-4)/2 >= -a0, "  T-4+2a0 =", float(T-4+2*a0), " (T-4)/2+a0 =", float((T-4)/2+a0))
# tau* <= 17/7 (L+2):  17/7 (L+log2(pi/2)) + 1 <= 17/7 (L+2)
print("tau* slack:", 17/7*(2-math.log2(math.pi/2)) - 1)
# F5(c): min over region {X>=0, tau<=0, X+4+2max(tau,-a0) >= T} of X+max(0,tau+gamma) equals G
Tf, a0f = float(T), float(a0)
for gam in [0.1, 0.5, 1.0, 1.2, 1.4, 1.5, 1.6, a0f]:
    X = np.linspace(0, 3, 1201)[:, None]; tau = np.linspace(-4, 0, 1601)[None, :]
    ok = X + 4 + 2*np.maximum(tau, -a0f) >= Tf - 1e-12
    val = np.where(ok, X + np.maximum(0, tau + gam), np.inf).min()
    G = max(0, gam + Tf/2 - 2)
    print(f"gamma={gam:.3f}  brute min={val:.4f}  G={G:.4f}")
