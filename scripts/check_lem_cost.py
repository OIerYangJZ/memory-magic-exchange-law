"""check_lem_cost.py -- random-instance check of the cost-weighted entropy lemma
(Lemma 13 of the manuscript, eq. (40)), and of the role of its hypothesis tau0 >= 1.

For an alphabet obeying  #{t(W) <= tau} <= c0 + c1 tau + c 2^tau eps^beta  the lemma bounds
    H(W) <= E[t] - q (beta L - log2 c - 1) + h2(q) + log2 N0 + 3 log2(1 + E t),
with tau0 = floor(beta L - log2 c), N0 = c0 + c1 tau0 + 1, q = Pr[t > tau0].

We build the densest alphabet allowed by the profile, put a distribution on it, and compare
the two sides.  With tau0 >= 1 the bound holds on every instance we try; without it the
bound is violated, which is why the hypothesis is now stated.
"""
import itertools, math, random

def shells(c0, c1, c, eps, beta, tmax):
    """n(t) = densest shell sizes consistent with the cumulative profile."""
    N = [c0 + c1 * t + c * (2.0 ** t) * eps ** beta for t in range(tmax + 1)]
    n, prev = [], 0.0
    for t in range(tmax + 1):
        cur = math.floor(max(N[t], 0.0))
        n.append(max(0, int(cur - prev))); prev = max(prev, cur)
    return n

def check(c0, c1, c, eps, beta, p_shape, tmax=40):
    L = math.log2(1 / eps)
    tau0 = math.floor(beta * L - math.log2(c))
    N0 = c0 + c1 * tau0 + 1
    if N0 <= 0:
        return tau0, None, None, None          # lemma's log2 N0 undefined
    n = shells(c0, c1, c, eps, beta, tmax)
    w = [p_shape(t) * n[t] for t in range(tmax + 1)]
    Z = sum(w)
    if Z <= 0:
        return tau0, None, None, None
    pt = [x / Z for x in w]                     # probability of shell t
    Et = sum(t * pt[t] for t in range(tmax + 1))
    # max entropy given the shell probabilities: uniform inside each shell
    H = sum(-pt[t] * math.log2(pt[t] / n[t]) for t in range(tmax + 1) if pt[t] > 0)
    q = sum(pt[t] for t in range(tmax + 1) if t > tau0)
    h2 = (-q * math.log2(q) - (1 - q) * math.log2(1 - q)) if 0 < q < 1 else 0.0
    rhs = Et - q * (beta * L - math.log2(c) - 1) + h2 + math.log2(N0) + 3 * math.log2(1 + Et)
    return tau0, H, rhs, H - rhs

random.seed(916)
bad_ok = bad_bad = 0
worst_ok = worst_bad = -1e9
for _ in range(20000):
    c0 = random.choice([0, 4, 8]); c1 = random.choice([0, 1, 2])
    c = 2.0 ** random.uniform(1, 14)
    eps = 2.0 ** -random.uniform(1.0, 12.0)
    beta = 2
    lam = random.uniform(0.05, 3.0)
    shape = lambda t, lam=lam: math.exp(-lam * t)
    tau0, H, rhs, gap = check(c0, c1, c, eps, beta, shape)
    if H is None:
        continue
    if tau0 >= 1:
        bad_ok += gap > 1e-9; worst_ok = max(worst_ok, gap)
    else:
        bad_bad += gap > 1e-9; worst_bad = max(worst_bad, gap)

print(f"tau0 >= 1  (hypothesis holds):  violations {bad_ok},  worst H-RHS = {worst_ok:+.4f} bits")
print(f"tau0 <= 0  (hypothesis fails):  violations {bad_bad},  worst H-RHS = {worst_bad:+.4f} bits")
print()
print("Violating instances found in the tau0 <= 0 region (hypothesis absent):")
random.seed(916)
shown = 0
for _ in range(20000):
    c0 = random.choice([0, 4, 8]); c1 = random.choice([0, 1, 2])
    c = 2.0 ** random.uniform(1, 14)
    eps = 2.0 ** -random.uniform(1.0, 12.0)
    lam = random.uniform(0.05, 3.0)
    tau0, H, rhs, gap = check(c0, c1, c, eps, 2, lambda t, lam=lam: math.exp(-lam * t))
    if H is not None and tau0 <= 0 and gap > 1e-9 and shown < 5:
        print(f"  c0={c0} c1={c1} c=2^{math.log2(c):5.2f} eps=2^-{-math.log2(eps):5.2f} "
              f"decay={lam:.2f}: tau0={tau0:3d}  H={H:.4f} > RHS={rhs:.4f}  (by {gap:.4f} bits)")
        shown += 1
if shown == 0:
    print("  none in this seed's draw")
