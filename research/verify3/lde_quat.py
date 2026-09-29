"""Referee check of App. clifford: denominator exponent vs T-count, and the
quaternion numerator 2w in Z[sqrt2,i] (exact, over Z[omega])."""
import sys, collections, math
sys.path.insert(0, '/Users/yangjinsey/Desktop/memory-magic-exchange-law/research')
from zomega import all_words, mul, add, conj, omega_pow, SQ2, ONE, OMEGA, val

tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
ws = all_words(tmax)
print("words", len(ws), "expected", 72 * 2 ** tmax - 48)

def sqrt2pow(e):
    r = ONE
    for _ in range(e):
        r = mul(r, SQ2)
    return r

def in_Zsqrt2_i(a):   # a0 + a1 w + a2 w^2 + a3 w^3 with a1 = a3 mod 2
    return (a[1] - a[3]) % 2 == 0

def quat_form(X):
    # X = [[A,B],[C,D]] ; need C = -conj(B), D = conj(A)
    A, B = X[0][0], X[0][1]
    C, D = X[1][0], X[1][1]
    return C == tuple(-x for x in conj(B)) and D == conj(A)

onepw = add(ONE, OMEGA)
tab = collections.Counter()
viol = collections.Counter()
bad_quat = collections.Counter()
maxk = collections.defaultdict(int)
for t, (k, M) in ws:
    tab[(t, k)] += 1
    maxk[t] = max(maxk[t], k)
    if k > t / 2 + 3:
        viol[t] += 1
    # quaternion numerator check
    if t % 2 == 0:
        e = t // 2 + 2 - k
        pre = sqrt2pow(e) if e >= 0 else None
    else:
        e = (t + 3) // 2 - k
        pre = mul(sqrt2pow(e), onepw) if e >= 0 else None
    if pre is None:
        bad_quat[(t, 'neg_exp')] += 1
        continue
    ok = False
    for j in range(8):
        p = mul(pre, omega_pow(j))
        X = [[mul(p, M[r][c]) for c in range(2)] for r in range(2)]
        if quat_form(X) and all(in_Zsqrt2_i(X[r][c]) for r in range(2) for c in range(2)):
            # reduced norm check in sigma1: |A|^2+|B|^2 = 4 n_t
            A, B = X[0][0], X[0][1]
            nr = abs(val(A)) ** 2 + abs(val(B)) ** 2
            nt = 2 ** (t / 2) if t % 2 == 0 else (2 + math.sqrt(2)) * 2 ** ((t - 1) / 2)
            if abs(nr - 4 * nt) < 1e-6 * nt:
                ok = True
                break
    if not ok:
        bad_quat[t] += 1

for t in range(tmax + 1):
    print(t, "k-distribution", sorted((kk, c) for (tt, kk), c in tab.items() if tt == t),
          "max k", maxk[t], "ceil(t/2)+1 =", -(-t // 2) + 1, "t/2+3 =", t / 2 + 3)
print("violations of k<=t/2+3:", dict(viol))
print("quaternion-numerator failures:", dict(bad_quat))
