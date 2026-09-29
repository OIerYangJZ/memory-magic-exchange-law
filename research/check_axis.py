"""Axis vectors nu(V) = (u*t + u t*, (u*t - u t*)/i, u u* - t t*) of Ad(V) z-hat, exact.
Checks: nu in Z[sqrt2]^3, <nu,nu> = 4^k in both embeddings; the map V -> nu has fibres of
size exactly 8 (right multiplication by powers of T) on each T-count ball."""
import collections, sys
from zomega import all_words, mul, conj, add, neg, val, _sig2
def real_part2(x):          # x + conj(x), lies in Z[sqrt2] (as element of Z[omega])
    return add(x, conj(x))
def imag_part2(x):          # (x - conj(x))/i : multiply by -i = omega^6 = -omega^2
    d = add(x, neg(conj(x)))
    return mul(d, (0, 0, -1, 0))
tmax = int(sys.argv[1])
ws = all_words(tmax)
fib = collections.Counter()
bad = 0
for tc, (k, m) in ws:
    u, t = m[0][0], m[1][0]
    ut = mul(conj(u), t)
    nu = (real_part2(ut), imag_part2(ut), add(mul(u, conj(u)), neg(mul(t, conj(t)))))
    for c in nu:                       # element of Z[sqrt2] <=> a1 = -a3 and a2 = 0
        if not (c[2] == 0 and c[1] == -c[3]):
            bad += 1
    n1 = sum(abs(val(c)) ** 2 for c in nu)
    n2 = sum(abs(val(_sig2(c))) ** 2 for c in nu)
    if abs(n1 - 4 ** k) > 1e-6 * 4 ** k or abs(n2 - 4 ** k) > 1e-6 * 4 ** k:
        bad += 1
    fib[(k, nu)] += 1
print("words", len(ws), "violations", bad)
print("fibre sizes over the ball:", collections.Counter(fib.values()))
