"""Identity coset: words within eps of some R_z(theta), fibred over the off-diagonal
entry t.  Checks the two inputs of Theorem A: (i) #t in the tube box vs (1+eps 2^k)^2,
(ii) max number of words over one t (the divisor factor), and the total vs eps^2 2^tau."""
import collections, sys, math
from zomega import all_words, val, canon, mul, omega_pow
tmax = int(sys.argv[1])
ws = all_words(tmax)
print("words", len(ws))
for L in (3, 4, 5, 6):
    eps = 2.0 ** -L
    fib = collections.defaultdict(set)      # (k, t-class mod units) -> set of words
    cum = collections.Counter()
    for tc, (k, m) in ws:
        t = m[1][0]
        b2 = abs(val(t)) ** 2 / 2 ** k
        if b2 <= eps * eps * (1 - eps * eps / 4) + 1e-15:
            # canonical t up to the 8 phases (u,t) -> omega^j (u,t)
            tk = min(tuple(mul(omega_pow(j), t)) for j in range(8))
            fib[(k, tk)].add(canon((k, m)))
            cum[tc] += 1
    per_k_t = collections.Counter(k for (k, _) in fib)
    maxfib = collections.defaultdict(int)
    for (k, tk), s in fib.items():
        maxfib[k] = max(maxfib[k], len(s))
    tot = 0
    rows = []
    for tc in range(tmax + 1):
        tot += cum[tc]
        rows.append((tc, tot, round(72 * eps * eps * 2 ** tc, 1)))
    print(f"L={L}: cumulative count vs 72 eps^2 2^tau:", [r for r in rows if r[1] or r[2] >= 1])
    print("   per k: #t-classes, bound (1+eps 2^k)^2, max words per t:",
          [(k, per_k_t[k], round((1 + eps * 2 ** k) ** 2, 1), maxfib[k]) for k in sorted(per_k_t)])
