"""Representations of N in Z[sqrt2] as x1^2+x2^2+x3^2 with x_i in Z[sqrt2], exact.
Elements a+b*sqrt2 stored as integer pairs (a,b).  sigma1 = a+b*s, sigma2 = a-b*s.
We pick N = sum of squares of a random triple with small sigma1 and large sigma2 parts,
then enumerate ALL representations, and study clustering of the sigma1 points at scales
below the mean spacing: are clusters of >= 4 points always coplanar (i.e. concyclic)?"""
import math, itertools, random, sys
import numpy as np
S = math.sqrt(2)

def elems(X, Y):
    """all a+b sqrt2 with |a+b s|<=X, |a-b s|<=Y"""
    out = []
    bmax = int((X + Y) / (2 * S)) + 1
    for b in range(-bmax, bmax + 1):
        lo = max(-X - b * S, -Y + b * S); hi = min(X - b * S, Y + b * S)
        for a in range(math.ceil(lo), math.floor(hi) + 1):
            out.append((a, b))
    return out

def sq(x):  # (a+b s)^2 = a^2+2b^2 + 2ab s
    a, b = x; return (a * a + 2 * b * b, 2 * a * b)

def is_square(n):
    """return all y in Z[sqrt2] with y^2 = n (n=(n0,n1)), exact"""
    n0, n1 = n
    s1 = n0 + n1 * S; s2 = n0 - n1 * S
    if s1 < -1e-9 or s2 < -1e-9: return []
    y1 = math.sqrt(max(s1, 0)); y2 = math.sqrt(max(s2, 0))
    res = []
    for e1 in (1, -1):
        for e2 in (1, -1):
            a = (e1 * y1 + e2 * y2) / 2; b = (e1 * y1 - e2 * y2) / (2 * S)
            A, B = round(a), round(b)
            if sq((A, B)) == (n0, n1) and (A, B) not in res: res.append((A, B))
    return res

def reps(N, X, Y):
    E = elems(X, Y)
    sqs = {x: sq(x) for x in E}
    out = set()
    for x1 in E:
        s1 = sqs[x1]
        for x2 in E:
            s2 = sqs[x2]
            rem = (N[0] - s1[0] - s2[0], N[1] - s1[1] - s2[1])
            for x3 in is_square(rem):
                out.add((x1, x2, x3))
    return sorted(out)

if __name__ == "__main__":
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    r, r2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6, float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
    # random triple with |x|_1 <= r/sqrt3, |x|_2 <= r2/sqrt3
    E = elems(r / math.sqrt(3), r2 / math.sqrt(3))
    tri = [random.choice(E) for _ in range(3)]
    N = tuple(map(sum, zip(*[sq(x) for x in tri])))
    R = reps(N, math.sqrt(N[0] + N[1] * S) + 1e-9, math.sqrt(N[0] - N[1] * S) + 1e-9)
    P = np.array([[a + b * S for (a, b) in x] for x in R])
    print("N =", N, " sigma1(N) =", N[0] + N[1] * S, " sigma2(N) =", N[0] - N[1] * S, " #points =", len(R))
    np.save("pts.npy", P); import pickle; pickle.dump((N, R), open("reps.pkl", "wb"))
