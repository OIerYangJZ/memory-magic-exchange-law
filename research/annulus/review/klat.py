"""Referee toolkit: exact arithmetic in K = Q(sqrt2), K^4, the maximal order O and O*.

A vector x in K^4 is stored in Q-coordinates as a tuple of 8 Fractions (p1..p4, q1..q4),
x_i = p_i + q_i sqrt2.  The R^8 metric (sigma1 (+) sigma2, standard in the basis 1,i,j,k) is the
trace form Tr<x,y> = 2 sum(p p' + 2 q q'), so all Gram determinants (covolume^2) are exact rationals.
Written independently of research/verify4.
"""
from fractions import Fraction as F
from math import gcd, sqrt
import itertools, random
import mpmath as mp

mp.mp.dps = 40
S2 = mp.sqrt(2)

# ------------------------------------------------------------------ K elements
class Kel:
    __slots__ = ('p', 'q')
    def __init__(self, p=0, q=0): self.p = F(p); self.q = F(q)
    def __add__(s, o): o = K(o); return Kel(s.p + o.p, s.q + o.q)
    __radd__ = __add__
    def __sub__(s, o): o = K(o); return Kel(s.p - o.p, s.q - o.q)
    def __rsub__(s, o): return K(o) - s
    def __neg__(s): return Kel(-s.p, -s.q)
    def __mul__(s, o): o = K(o); return Kel(s.p * o.p + 2 * s.q * o.q, s.p * o.q + s.q * o.p)
    __rmul__ = __mul__
    def norm(s): return s.p * s.p - 2 * s.q * s.q
    def inv(s):
        n = s.norm(); assert n != 0
        return Kel(s.p / n, -s.q / n)
    def __truediv__(s, o): return s * K(o).inv()
    def iszero(s): return s.p == 0 and s.q == 0
    def s1(s): return mp.mpf(s.p.numerator) / s.p.denominator + S2 * (mp.mpf(s.q.numerator) / s.q.denominator)
    def s2(s): return mp.mpf(s.p.numerator) / s.p.denominator - S2 * (mp.mpf(s.q.numerator) / s.q.denominator)
    def isint(s): return s.p.denominator == 1 and s.q.denominator == 1
    def __repr__(s): return f"({s.p}+{s.q}r2)"
def K(x): return x if isinstance(x, Kel) else Kel(x)
R2 = Kel(0, 1)
UNIT = Kel(1, 1)

# K^4 vectors <-> Q^8
def kv(x8): return [Kel(x8[i], x8[i + 4]) for i in range(4)]
def q8(v): return tuple([F(c.p) for c in v] + [F(c.q) for c in v])
def kdot(u, v): return sum((a * b for a, b in zip(u, v)), Kel(0))
def kscale(c, v): return [c * x for x in v]
def kadd(u, v): return [a + b for a, b in zip(u, v)]
def ksub(u, v): return [a - b for a, b in zip(u, v)]
def sig1(v): return mp.matrix([c.s1() for c in v])
def sig2(v): return mp.matrix([c.s2() for c in v])
def nrm(m): return mp.sqrt(sum(m[i] ** 2 for i in range(len(m))))
def h(v): return nrm(sig1(v)) * nrm(sig2(v))
def trace_form(x8, y8):
    return 2 * sum(x8[i] * y8[i] for i in range(4)) + 4 * sum(x8[i] * y8[i] for i in range(4, 8))

# ------------------------------------------------------------------ integer linear algebra
def hnf_rows(rows):
    """Z-basis (list of integer rows) of the Z-span of integer rows (row echelon, gcd steps)."""
    rows = [list(r) for r in rows if any(r)]
    if not rows: return []
    ncol = len(rows[0]); out = []; piv = 0
    for c in range(ncol):
        while True:
            nz = [r for r in rows if r[c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(r[c]))
            p = nz[0]
            new = [p]
            for r in rows:
                if r is p: continue
                if r[c] != 0:
                    k = r[c] // p[c]
                    r = [a - k * b for a, b in zip(r, p)]
                if any(r): new.append(r)
            rows = new
        nz = [r for r in rows if r[c] != 0]
        if nz:
            p = nz[0]; out.append(p); rows = [r for r in rows if r is not p]
    return out

def zspan_basis(vecs):
    """Z-basis of the Z-span of rational vectors."""
    den = 1
    for v in vecs:
        for c in v: den = den * c.denominator // gcd(den, c.denominator)
    rows = [[int(c * den) for c in v] for v in vecs]
    B = hnf_rows(rows)
    return [tuple(F(c, den) for c in r) for r in B]

def int_kernel(M):
    """integer kernel {z in Z^n : sum_i z_i M[i] = 0} of rational rows M[i] (n rows)."""
    n = len(M); m = len(M[0])
    den = 1
    for r in M:
        for c in r: den = den * F(c).denominator // gcd(den, F(c).denominator)
    rows = [[int(F(c) * den) for c in M[i]] + [1 if j == i else 0 for j in range(n)] for i in range(n)]
    B = hnf_rows(rows)
    return [r[m:] for r in B if not any(r[:m])]

def combo(z, basis):
    return tuple(sum((F(zi) * b[k] for zi, b in zip(z, basis)), F(0)) for k in range(len(basis[0])))

def gram(basis): return [[trace_form(x, y) for y in basis] for x in basis]
def det(Mx):
    M = [[F(c) for c in r] for r in Mx]; n = len(M); d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f: M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d
def covol2(basis): return det(gram(basis))           # exact covol^2
def covol(basis): return mp.sqrt(mp.mpf(covol2(basis).numerator) / covol2(basis).denominator)

# ------------------------------------------------------------------ the order
h_ = F(1, 2)
OBAS = [  # O_K-basis of O, K-coordinates in (1,i,j,k)
    [Kel(1), Kel(0), Kel(0), Kel(0)],
    [Kel(0, h_), Kel(0, h_), Kel(0), Kel(0)],
    [Kel(0, h_), Kel(0), Kel(0, h_), Kel(0)],
    [Kel(h_), Kel(h_), Kel(h_), Kel(h_)],
]
O_Z = [q8(b) for b in OBAS] + [q8(kscale(R2, b)) for b in OBAS]     # Z-basis of O
OK4_Z = [tuple(F(int(i == j)) for j in range(8)) for i in range(8)]   # Z-basis of O_K^4

def dual_Z(Zbasis_kvecs_OK):
    """O* = {x : <x, b> in O_K for the O_K-basis b}; returns Z-basis in Q^8."""
    # map x -> (<x,b_i>)_i  in Q^8 (p parts then q parts); O* = preimage of Z^8
    cols = []
    for e in OK4_Z:
        v = kv(e); vals = [kdot(v, b) for b in Zbasis_kvecs_OK]
        cols.append([c.p for c in vals] + [c.q for c in vals])
    # cols[k] = image of k-th standard vector; matrix A with x -> x A
    A = cols
    # invert A (8x8 rational)
    n = 8; M = [list(A[i]) + [F(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    Ainv = [row[n:] for row in M]
    # x A = e_k  <=>  x = e_k Ainv ; rows of Ainv are the basis of O*
    return [tuple(r) for r in Ainv]
OSTAR_Z = dual_Z(OBAS)

def in_lattice(x8, basis):
    """is x8 in the Z-span of basis (full-rank case solved by rational elimination)"""
    B = [list(b) for b in basis]
    n = len(B)
    # solve z B = x (least squares exact via normal eq on Gram with standard dot)
    G = [[sum(a * b for a, b in zip(B[i], B[j])) for j in range(n)] for i in range(n)]
    r = [sum(a * b for a, b in zip(x8, B[i])) for i in range(n)]
    M = [G[i] + [r[i]] for i in range(n)]
    for c in range(n):
        p = next(rr for rr in range(c, n) if M[rr][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for rr in range(n):
            if rr != c and M[rr][c] != 0:
                f = M[rr][c]; M[rr] = [a - f * b for a, b in zip(M[rr], M[c])]
    z = [M[i][n] for i in range(n)]
    back = combo(z, basis)
    return all(z_.denominator == 1 for z_ in z) and back == tuple(x8)

# ------------------------------------------------------------------ K-subspaces
def kgram_solve(G, rhs):
    """solve G c = rhs over K (small dense)"""
    n = len(G); M = [list(G[i]) + [rhs[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if not M[r][c].iszero()); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and not M[r][c].iszero():
                f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]

def proj(x, basisK):
    """K-orthogonal projection of x onto span_K(basisK)"""
    G = [[kdot(u, v) for v in basisK] for u in basisK]
    rhs = [kdot(u, x) for u in basisK]
    c = kgram_solve(G, rhs)
    out = [Kel(0)] * 4
    for ci, u in zip(c, basisK): out = kadd(out, kscale(ci, u))
    return out

def kperp(basisK):
    """K-basis of the orthogonal complement of span_K(basisK) in K^4"""
    out = []
    for e in range(4):
        v = [Kel(int(i == e)) for i in range(4)]
        w = ksub(v, proj(v, basisK)) if basisK else v
        if all(c.iszero() for c in w): continue
        cand = basisK + out
        w2 = ksub(v, proj(v, cand)) if cand else v
        if not all(c.iszero() for c in w2): out.append(w2)
        if len(out) + len(basisK) == 4: break
    return out

def qspan_of(basisK):
    """Q-basis (Q^8 rows) of span_K(basisK)"""
    return [q8(b) for b in basisK] + [q8(kscale(R2, b)) for b in basisK]

def lattice_cap_subspace(Lbasis, basisK):
    """Z-basis of L cap span_K(basisK): integer z with z.L orthogonal (trace form) to the complement"""
    comp = qspan_of(kperp(basisK))
    M = [[trace_form(x, c) for c in comp] for x in Lbasis]
    ker = int_kernel(M)
    return [combo(z, Lbasis) for z in ker]

def lattice_proj(Lbasis, basisK):
    """Z-basis of P_{span basisK}(L)"""
    imgs = [q8(proj(kv(x), basisK)) for x in Lbasis]
    return zspan_basis(imgs)

def wedge(vs):
    """Pluecker coordinates (K) of v1 ^ ... ^ vk"""
    k = len(vs); out = []
    for I in itertools.combinations(range(4), k):
        M = [[vs[r][c] for c in I] for r in range(k)]
        out.append(kdet(M))
    return out
def kdet(M):
    n = len(M)
    if n == 1: return M[0][0]
    s = Kel(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        t = M[0][j] * kdet(minor)
        s = s + t if j % 2 == 0 else s - t
    return s

def rand_OK(rng, B=3): return Kel(rng.randint(-B, B), rng.randint(-B, B))
def rand_vec(rng, B=3): return [rand_OK(rng, B) for _ in range(4)]
def rand_lattice_vec(rng, basis, B=3):
    z = [rng.randint(-B, B) for _ in basis]
    return combo(z, basis)

def fn_index(n8, Lbasis):
    """index of the Z-module {<n,x> : x in L} in O_K (0 if not full rank)."""
    n = kv(n8)
    vals = [kdot(n, kv(x)) for x in Lbasis]
    assert all(v.isint() for v in vals), "values not in O_K"
    rows = [[int(v.p), int(v.q)] for v in vals]
    B = hnf_rows(rows)
    if len(B) < 2: return 0
    return abs(B[0][0] * B[1][1] - B[0][1] * B[1][0])

# ------------------------------------------------------------------ enumeration (floats)
def emb8(x8):
    v = kv(x8); a = sig1(v); b = sig2(v)
    return [float(a[i]) for i in range(4)] + [float(b[i]) for i in range(4)]

import numpy as np
def lll(B, delta=0.99):
    """float LLL on rows of B (numpy)"""
    B = np.array(B, dtype=float).copy(); n = B.shape[0]
    def gso(B):
        Bs = np.zeros_like(B); mu = np.zeros((n, n))
        for i in range(n):
            v = B[i].copy()
            for j in range(i):
                mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j]); v -= mu[i, j] * Bs[j]
            Bs[i] = v
        return Bs, mu
    k = 1; Bs, mu = gso(B)
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k, j])
            if q: B[k] -= q * B[j]; Bs, mu = gso(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]): k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]; Bs, mu = gso(B); k = max(k - 1, 1)
    return B

def enum_quadform(M, bound):
    """all integer z (as tuples) with z^T M z <= bound, M positive definite (numpy n x n)"""
    n = M.shape[0]
    L = np.linalg.cholesky(M)            # M = L L^T ; z^T M z = |L^T z|^2
    U = L.T                              # upper triangular
    out = []
    z = [0] * n
    def rec(i, partial):
        # coordinates i..n-1 of U z ; choose z_i given z_{i+1..}
        c = sum(U[i, j] * z[j] for j in range(i + 1, n))
        r2 = bound - partial
        if r2 < -1e-12: return
        r = np.sqrt(max(r2, 0.0)) / abs(U[i, i])
        center = -c / U[i, i]
        lo = int(np.ceil(center - r - 1e-12)); hi = int(np.floor(center + r + 1e-12))
        for zi in range(lo, hi + 1):
            z[i] = zi
            val = (U[i, i] * zi + c) ** 2
            if i == 0:
                if partial + val <= bound + 1e-9: out.append(tuple(z))
            else:
                rec(i - 1, partial + val)
        z[i] = 0
    rec(n - 1, 0.0)
    return out

def enum_box(M, bound, maxpts=4_000_000):
    """all integer z with z^T M z <= bound (vectorised: LLL-reduce w.r.t. M, then scan a box).
    Returns an (k x n) int array, or None if the box is too large."""
    n = M.shape[0]
    U = np.linalg.cholesky(M).T                       # Q(z) = |U z|^2
    Bred = lll(U.T)                                   # rows = reduced basis vectors (U T)^T
    T = np.rint(np.linalg.solve(U, Bred.T)).astype(np.int64)   # U T = Bred^T
    G = Bred @ Bred.T
    Gi = np.linalg.inv(G)
    r = np.floor(np.sqrt(bound * np.diag(Gi)) + 1e-9).astype(int)
    size = np.prod(2 * r + 1)
    if size > maxpts: return None
    axes = [np.arange(-ri, ri + 1) for ri in r]
    W = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, n)
    q = np.einsum('ij,jk,ik->i', W, G, W)
    W = W[q <= bound * (1 + 1e-9)]
    return W @ T.T

def reduce_exact(Lb, delta=F(99, 100)):
    """exact LLL (Fractions, on the trace-form Gram matrix) of a Z-basis given as Q^8 rows"""
    B = [list(x) for x in Lb]; n = len(B)
    def ip(x, y): return trace_form(x, y)
    def gso():
        Bs = []; mu = [[F(0)] * n for _ in range(n)]; nb = []
        for i in range(n):
            v = list(B[i])
            for j in range(i):
                mu[i][j] = ip(B[i], Bs[j]) / nb[j]
                v = [a - mu[i][j] * c for a, c in zip(v, Bs[j])]
            Bs.append(v); nb.append(ip(v, v))
        return Bs, mu, nb
    k = 1; Bs, mu, nb = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * c for a, c in zip(B[k], B[j])]; Bs, mu, nb = gso()
        if nb[k] >= (delta - mu[k][k - 1] ** 2) * nb[k - 1]: k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu, nb = gso(); k = max(k - 1, 1)
    return [tuple(x) for x in B]
