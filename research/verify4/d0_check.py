"""Referee check of D0 and of the D6(b) duality covol(V^perp cap O*) ~ H_V.
Exact arithmetic in K = Q(sqrt2): elements are pairs (p,q) = p + q*sqrt2 with rationals."""
from fractions import Fraction as F
import itertools, random, math
import numpy as np

s2 = math.sqrt(2)
# O_K-basis of the maximal order, coordinates in (1,i,j,k); entries p+q sqrt2 stored as (p,q)
h = F(1, 2)
B = [[(1,0),(0,0),(0,0),(0,0)],
     [(0,h),(0,h),(0,0),(0,0)],          # (1+i)/sqrt2
     [(0,h),(0,0),(0,h),(0,0)],          # (1+j)/sqrt2
     [(h,0),(h,0),(h,0),(h,0)]]          # (1+i+j+k)/2
def emb(v, sgn):
    return np.array([float(p) + sgn * float(q) * s2 for (p, q) in v])
def mulsqrt2(v):  # multiply K-vector by sqrt2
    return [(2 * q, p) for (p, q) in v]
Zbasis = []
for b in B:
    Zbasis.append(b); Zbasis.append(mulsqrt2(b))
def R8(v): return np.concatenate([emb(v, 1), emb(v, -1)])
MO = np.array([R8(v) for v in Zbasis])           # rows: Z-basis of O in R^8
covO = abs(np.linalg.det(MO))
# trace dual
Odual = np.linalg.inv(MO).T                      # rows: Z-basis of O^vee
# O* = delta * O^vee with delta = 2 sqrt2 : sigma1 -> *2sqrt2, sigma2 -> *(-2sqrt2)
D = np.diag([2*s2]*4 + [-2*s2]*4)
Ostar = Odual @ D
covOs = abs(np.linalg.det(Ostar))
print("covol O =", covO, " covol O* =", covOs, " ratio", covO/covOs)
# check O* subset O_K^4 : coordinates of each basis vector are in O_K  <=> (x1+x2)/2 in Z, (x1-x2)/(2 sqrt2) in Z
def in_OK4(x):
    a = (x[:4] + x[4:]) / 2; b = (x[:4] - x[4:]) / (2 * s2)
    return np.allclose(a, np.round(a), atol=1e-9) and np.allclose(b, np.round(b), atol=1e-9)
print("O* subset O_K^4:", all(in_OK4(x) for x in Ostar))
# 4O subset O* : coordinates of 4O in basis of O* are integers
C = np.linalg.solve(Ostar.T, (4 * MO).T)
print("4O subset O*:", np.allclose(C, np.round(C), atol=1e-8))
C2 = np.linalg.solve(Ostar.T, (2 * MO).T)
print("2O subset O*:", np.allclose(C2, np.round(C2), atol=1e-8))
# O_K^4 subset O
E = []
for m in range(4):
    e = [(0,0)]*4; e[m] = (1,0); E.append(R8(e)); E.append(R8(mulsqrt2(e)))
C3 = np.linalg.solve(MO.T, np.array(E).T)
print("O_K^4 subset O:", np.allclose(C3, np.round(C3), atol=1e-8))

# ---- D0: M = O cap n^perp, covol vs h(n), n primitive in O*
def lattice_covol(rows):
    G = rows @ rows.T
    return math.sqrt(abs(np.linalg.det(G)))
def ok_coords(x):  # R^8 vector with O_K coordinates -> integer coords (a,b) for each of 4 coords
    a = (x[:4] + x[4:]) / 2; b = (x[:4] - x[4:]) / (2 * s2)
    return np.round(a).astype(int), np.round(b).astype(int)
random.seed(1)
# express O* basis vectors integrally
res = []
for trial in range(12):
    c = [random.randint(-6, 6) for _ in range(8)]
    n = sum(ci * Ostar[i] for i, ci in enumerate(c))
    if np.allclose(n, 0): continue
    # f_n : O -> O_K ;  matrix of Z-linear map Z^8 -> Z^2 (O_K = Z + Z sqrt2)
    vals = []
    for x in MO:
        v1 = float(np.dot(n[:4], x[:4])); v2 = float(np.dot(n[4:], x[4:]))
        p = (v1 + v2) / 2; q = (v1 - v2) / (2 * s2)
        vals.append((round(p), round(q)))
        assert abs(p - round(p)) < 1e-7 and abs(q - round(q)) < 1e-7
    A = np.array(vals).T  # 2 x 8
    # image index: gcd of 2x2 minors = index of f_n(O) in O_K
    minors = [abs(A[0,i]*A[1,j]-A[0,j]*A[1,i]) for i in range(8) for j in range(i+1,8)]
    idx = 0
    for m in minors: idx = math.gcd(idx, int(m))
    # primitive part: divide n by generator of the ideal -> we only test primitive ones (idx = 1)
    hn = np.linalg.norm(n[:4]) * np.linalg.norm(n[4:])
    # kernel lattice via sympy smith normal form
    rows = [list(map(int, A[:, i])) + [1 if j == i else 0 for j in range(8)] for i in range(8)]
    # integer row reduction on first 2 columns
    def reduce(rows, col, start):
        piv = start
        while True:
            nz = [r for r in range(piv, len(rows)) if rows[r][col] != 0]
            if not nz: return piv
            r0 = min(nz, key=lambda r: abs(rows[r][col]))
            rows[piv], rows[r0] = rows[r0], rows[piv]
            done = True
            for r in range(piv + 1, len(rows)):
                if rows[r][col]:
                    qq = rows[r][col] // rows[piv][col]
                    rows[r] = [a - qq * b for a, b in zip(rows[r], rows[piv])]
                    if rows[r][col]: done = False
            if done: return piv + 1
    p1 = reduce(rows, 0, 0); p2 = reduce(rows, 1, p1)
    ker = [r[2:] for r in rows[p2:]]
    Kmat = np.array(ker) @ MO   # rows in R^8
    cM = lattice_covol(Kmat)
    res.append((idx, hn, cM, cM / hn * idx))
for r in res: print("index f_n(O)=%d  h(n)=%.3f  covol(M)=%.4f  covol(M)*idx/h(n)=%.6f" % r)
print("predicted covol(O)/covol(O_K) =", covO / (2 * s2))
