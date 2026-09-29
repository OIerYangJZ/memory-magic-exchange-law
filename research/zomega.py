"""Exact single-qubit Clifford+T arithmetic over Z[omega], omega = e^{i pi/4}.

A unitary is stored as (k, M) with M a 2x2 matrix over Z[omega] and the unitary equal
to M / sqrt2^k.  Elements of Z[omega] are 4-tuples (a0,a1,a2,a3) = sum a_i omega^i,
omega^4 = -1.  Everything here is exact; floats are used only for embeddings.
"""
import cmath
import itertools

W8 = cmath.exp(1j * cmath.pi / 4)


def mul(a, b):
    c = [0, 0, 0, 0]
    for i in range(4):
        if a[i] == 0:
            continue
        for j in range(4):
            s = i + j
            if s < 4:
                c[s] += a[i] * b[j]
            else:
                c[s - 4] -= a[i] * b[j]
    return tuple(c)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def conj(a):  # complex conjugation: omega -> omega^7 = -omega^3
    return (a[0], -a[3], -a[2], -a[1])


def sig2(a):  # Galois sqrt2 -> -sqrt2 : omega -> omega^3 (then |.| is the sigma2 absolute value)
    return (a[0], -a[3], a[2], -a[1]) if False else _sig2(a)


def _sig2(a):
    # omega -> omega^3: omega^0->1, omega^1->omega^3, omega^2->omega^6=-omega^2, omega^3->omega^9=omega
    return (a[0], a[3], -a[2], a[1])


def val(a, which=1):
    z = W8 if which == 1 else W8 ** 3
    return sum(a[i] * z ** i for i in range(4))


SQ2 = (0, 1, 0, -1)          # omega - omega^3 = sqrt2
ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)
OMEGA = (0, 1, 0, 0)


def div_sqrt2(a):
    """a / sqrt2 if in Z[omega], else None.  a/sqrt2 = a*sqrt2/2."""
    b = mul(a, SQ2)
    if all(x % 2 == 0 for x in b):
        return tuple(x // 2 for x in b)
    return None


def mat_mul(A, B):
    (ka, a), (kb, b) = A, B
    m = [[ZERO, ZERO], [ZERO, ZERO]]
    for i in range(2):
        for j in range(2):
            m[i][j] = add(mul(a[i][0], b[0][j]), mul(a[i][1], b[1][j]))
    return normalize((ka + kb, m))


def normalize(A):
    k, m = A
    while k > 0:
        d = [[div_sqrt2(m[i][j]) for j in range(2)] for i in range(2)]
        if any(d[i][j] is None for i in range(2) for j in range(2)):
            break
        m, k = d, k - 1
    return (k, [list(r) for r in m])


def omega_pow(j):
    j %= 8
    s = 1 if j < 4 else -1
    e = [0, 0, 0, 0]
    e[j % 4] = s
    return tuple(e)


def canon(A):
    """Canonical key modulo global phase omega^j (the only phases in the group)."""
    k, m = A
    best = None
    for j in range(8):
        p = omega_pow(j)
        key = (k,) + tuple(x for i in range(2) for jj in range(2) for x in mul(p, m[i][jj]))
        if best is None or key < best:
            best = key
    return best


H = normalize((1, [[ONE, ONE], [ONE, neg(ONE)]]))
T = (0, [[ONE, ZERO], [ZERO, OMEGA]])
S = (0, [[ONE, ZERO], [ZERO, omega_pow(2)]])
I2 = (0, [[ONE, ZERO], [ZERO, ONE]])


def cliffords():
    seen = {canon(I2): I2}
    frontier = [I2]
    while frontier:
        nf = []
        for A in frontier:
            for G in (H, S):
                B = mat_mul(G, A)
                c = canon(B)
                if c not in seen:
                    seen[c] = B
                    nf.append(B)
        frontier = nf
    return list(seen.values())


def ma_words(tmax):
    """Yield (tcount, unitary) for every Clifford+T unitary mod phase with minimal
    T-count <= tmax, via the Matsumoto-Amano normal form (T|e)(HT|SHT)^* C."""
    CL = cliffords()
    assert len(CL) == 24
    HT = mat_mul(H, T)
    SHT = mat_mul(S, HT)
    layer = [(0, C) for C in CL]
    for C in CL:
        yield 0, C
    for n in range(1, tmax + 1):
        new = []
        for _, W in layer:
            for G in (HT, SHT):
                new.append((n, mat_mul(G, W)))
        layer = new
        for n_, W in layer:
            yield n_, W            # no leading T: T-count n
        if n - 1 >= 0:
            pass
    # leading-T words handled separately below


def all_words(tmax):
    """All words with T-count <= tmax, including the leading-T branch."""
    CL = cliffords()
    HT = mat_mul(H, T)
    SHT = mat_mul(S, HT)
    out = [(0, C) for C in CL]
    layer = list(CL)
    for n in range(0, tmax + 1):
        if n > 0:
            layer = [mat_mul(G, W) for W in layer for G in (HT, SHT)]
            out.extend((n, W) for W in layer)
        if n + 1 <= tmax:                       # leading T on top of n syllables
            out.extend((n + 1, mat_mul(T, W)) for W in layer)
    return out
