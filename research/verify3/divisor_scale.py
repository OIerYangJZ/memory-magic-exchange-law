"""Scale of the worst-case divisor factor in Thm clifford at eps=1e-10:
max d(n) = #ideal divisors of an ideal n of Z[sqrt2] with N(n) <= X (exponents
non-increasing along prime-ideal norms, the standard highly-composite reduction)."""
import math
def primerange(a, b):
    s = bytearray([1]) * b; s[0:2] = b"\x00\x00"
    for i in range(2, int(b ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(a, b) if s[i]]
def prime_ideal_norms(limit):
    out = [2]
    for p in primerange(3, limit):
        out += [p, p] if p % 8 in (1, 7) else [p * p]
    return sorted(out)
def maxdiv(X, norms):
    best = [1, None]
    def dfs(i, rem, val, emax, exps):
        if val > best[0]:
            best[0], best[1] = val, list(exps)
        if i >= len(norms) or emax == 0:
            return
        q, rr = norms[i], rem
        for e in range(1, emax + 1):
            rr //= q
            if rr < 1:
                break
            exps.append((q, e)); dfs(i + 1, rr, val * (e + 1), e, exps); exps.pop()
    dfs(0, X, 1, 400, [])
    return best
for tau in (66, 99, 3 * 33):
    X = 16 * 2 ** tau
    d, ex = maxdiv(X, prime_ideal_norms(2000))
    print(f"tau={tau}: max d(n) over N(n)<=16*2^tau is {d} (log2 {math.log2(d):.1f}); exps {ex[:6]}...")
    print(f"   log2(576*8*d) = {math.log2(576*8*d):.1f}  vs log2 K_eps = 32.87 at eps=1e-10")
