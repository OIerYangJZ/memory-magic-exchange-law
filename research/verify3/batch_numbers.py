"""Referee check of Prop. batch numbers at eps=1e-10, m=1e4, r=2."""
import math
eps, m, r = 1e-10, 10**4, 2
L = math.log2(1 / eps)
def cost(b, g):
    return (8 * (2 ** g - 1) + 4 * (b - 1)) / g
def best(b):
    return min((cost(b, g), g) for g in range(1, 20))
crit = {
  "paper total, op-norm, eps/2: b=ceil(L+log2(2 pi m r))": math.ceil(L + math.log2(2 * math.pi * m * r)),
  "total, dproj via Lem loc, full eps: pi m r 2^-b/2<=eps": math.ceil(L + math.log2(math.pi * m * r / 2)),
  "blockwise (g=4) op-norm eps/2 over r rounds: ceil(L+log2(2 pi g r))": math.ceil(L + math.log2(2 * math.pi * 4 * r)),
  "coordinatewise op-norm eps/2: ceil(L+log2(2 pi r))": math.ceil(L + math.log2(2 * math.pi * r)),
  "coordinatewise dproj full eps: pi r 2^-b /2 <= eps": math.ceil(L + math.log2(math.pi * r / 2)),
}
for k, b in crit.items():
    c, g = best(b)
    print(f"{k:75s} b={b:3d}  best g={g}  T/share={c:.2f}   (g=4: {cost(b,4):.2f}, g=3: {cost(b,3):.2f})")
print("asymptotics: per-share cost / (4L/log2 L) with 2^g ~ L/log L, log(mr)=o(L)")
for LL in (50, 100, 1000, 10**4, 10**6):
    b = LL + 20
    c, g = best(b) if LL < 2000 else (min((cost(b, g), g) for g in range(1, 40)))
    print(f"L={LL:8d} best g={g:2d} cost={c:10.1f}  ratio to 4L/log2L = {c/(4*LL/math.log2(LL)):.3f}  cost/L={c/LL:.4f}")
