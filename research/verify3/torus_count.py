"""Referee check of Thm clifford / Cor typical: words near the torus {R_z}, their
off-diagonal numerators B=2c+2di, and fibre sizes; compare with tau+1+eps^2 2^tau."""
import sys, os, math, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
from zomega import all_words, mul, add, conj, omega_pow, SQ2, ONE, OMEGA, val

tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 13
ws = all_words(tmax)
def sqrt2pow(e):
    r = ONE
    for _ in range(e):
        r = mul(r, SQ2)
    return r
def quat_form(X):
    A, B = X[0][0], X[0][1]; C, D = X[1][0], X[1][1]
    return C == tuple(-x for x in conj(B)) and D == conj(A)
onepw = add(ONE, OMEGA)
recs = []   # (t, dist_to_torus, A, B)
for t, (k, M) in ws:
    e = t // 2 + 2 - k if t % 2 == 0 else (t + 3) // 2 - k
    pre = sqrt2pow(e) if t % 2 == 0 else mul(sqrt2pow(e), onepw)
    for j in range(8):
        p = mul(pre, omega_pow(j))
        X = [[mul(p, M[r][c]) for c in range(2)] for r in range(2)]
        if quat_form(X):
            break
    A, B = X[0][0], X[0][1]
    a, b = abs(val(A)), abs(val(B))
    n = math.hypot(a, b)
    # min_theta dproj(W,R_z) = sqrt(2-2|a_hat|)
    d = math.sqrt(max(0.0, 2 - 2 * a / n))
    # canonical sign class of (A,B)
    key = min((A, B), (tuple(-x for x in A), tuple(-x for x in B)))
    recs.append((t, d, key[0], key[1]))

for L in (3, 4, 5, 6):
    eps = 2.0 ** -L
    near = [r for r in recs if r[1] <= eps]
    print(f"L={L} eps={eps}")
    for tau in range(0, tmax + 1):
        sel = [r for r in near if r[0] <= tau]
        nB = len(set(r[3] for r in sel))
        fib = collections.Counter(r[3] for r in sel)
        nonTpow = sum(1 for r in sel if r[3] != (0, 0, 0, 0))
        mint = min((r[0] for r in near if r[3] != (0, 0, 0, 0)), default=None)
        print(f"  tau={tau:2d} N={len(sel):6d} nonTpow={nonTpow:6d} distinctB={nB:5d} maxfibre={max(fib.values()) if fib else 0:3d}"
              f"  ref(tau+1+eps^2 2^tau)={tau+1+eps*eps*2**tau:9.1f}  ratio={len(sel)/(tau+1+eps*eps*2**tau):7.2f}  Haar72eps^2 2^tau={72*eps*eps*2**tau:8.1f}")
    print("  least T-count of a non-T-power word within eps of torus:", mint, " 2L-8 =", 2 * L - 8)
