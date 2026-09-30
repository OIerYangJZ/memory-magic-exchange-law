#!/usr/bin/env python3
"""Referee check of NOTES.md sec. 3-4: exact arithmetic of sphere sections Y_{g,s} (needs numpy; run on the server).

O_K = Z[sqrt2] elements are integer pairs (a, b) = a + b sqrt2.  Quaternions are stored in doubled coordinates
X = 2x in O_K^4 (basis 1, i, j, k).  Membership of x = X/2 in the maximal order O of main.tex eq. (order):
    X1 - X3 in sqrt2 O_K,  X2 - X3 in sqrt2 O_K,  X0 - X1 - X2 + X3 in 2 O_K.
For l in O_K^+, s = S/2 in (1/2)O_K and n = n_tau:  Y_{g,s} = {w in O : nrd w = n, <w,g> = s},  nrd g = l.
z = w gbar has Re z = s, nrd z = n l;  conversely w = z g / l.  We enumerate Z = {z in O : Re z = s, nrd z = nl}
exactly and count |Y_{g,s}| = #{z in Z : z g / l in O} for every g of norm l.

Usage:
  check_sections.py known          # the two clusters of NOTES sec. 4 (tau = 21 and 26)
  check_sections.py census Hmax mu tau0 tau1
"""
import sys, math
import numpy as np

S2 = math.sqrt(2.0)


def haar(rho):
    th = 2 * math.asin(rho / 2)
    return 2 * (th - math.sin(th) * math.cos(th)) / math.pi


def rho_for(mu, tau):
    N = 72.0 * 2.0 ** tau - 48
    lo, hi = 0.0, 1.4
    for _ in range(200):
        md = (lo + hi) / 2
        if N * haar(md) < mu:
            lo = md
        else:
            hi = md
    return hi


def n_tau(t):
    return (2 ** (t // 2), 0) if t % 2 == 0 else (2 ** ((t + 1) // 2), 2 ** ((t - 1) // 2))


def s1(x): return x[0] + x[1] * S2
def s2(x): return x[0] - x[1] * S2
def mul(x, y): return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def nrm(x): return x[0] * x[0] - 2 * x[1] * x[1]


def cands(r1, r2):
    """all (a, b) with |a + b sqrt2| <= r1 and |a - b sqrt2| <= r2."""
    out = []
    bmax = int((r1 + r2) / (2 * S2)) + 1
    for b in range(-bmax, bmax + 1):
        lo = math.ceil(max(-r1 - b * S2, -r2 + b * S2) - 1e-9)
        hi = math.floor(min(r1 - b * S2, r2 + b * S2) + 1e-9)
        for a in range(lo, hi + 1):
            out.append((a, b))
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def sqrt_ok(p, q):
    """vectorised: returns list of (mask, x, y) with (x + y sqrt2)^2 = p + q sqrt2, both square roots up to sign."""
    a1 = p + q * S2; a2 = p - q * S2
    ok = (a1 >= -0.5) & (a2 >= -0.5)
    r1 = np.sqrt(np.maximum(a1, 0)); r2 = np.sqrt(np.maximum(a2, 0))
    res = []
    for sg in (1, -1):
        x = np.rint((r1 + sg * r2) / 2).astype(np.int64)
        y = np.rint((r1 - sg * r2) / (2 * S2)).astype(np.int64)
        m = ok & (x * x + 2 * y * y == p) & (2 * x * y == q)
        res.append((m, x, y))
    return res


def inO(X):
    """X: array (..., 4, 2) doubled coordinates."""
    c1 = (X[..., 1, 0] - X[..., 3, 0]) % 2 == 0
    c2 = (X[..., 2, 0] - X[..., 3, 0]) % 2 == 0
    t = X[..., 0, :] - X[..., 1, :] - X[..., 2, :] + X[..., 3, :]
    return c1 & c2 & (t[..., 0] % 2 == 0) & (t[..., 1] % 2 == 0)


def primitive(X):
    """x = X/2 not in sqrt2*O."""
    ev = np.all(X[..., :, 0] % 2 == 0, axis=-1)
    Y = np.stack([X[..., :, 1], X[..., :, 0] // 2], axis=-1)
    return ~(ev & inO(Y))


def reps(target, first=None):
    """all X in O_K^4 (doubled) with sum X_i^2 = target (in O_K) and X/2 in O; optionally X0 fixed = first."""
    T1, T2 = s1(target), s2(target)
    if first is not None:   # the remaining three coordinates are bounded by the remainder, not by the target
        T1 -= s1(first) ** 2; T2 -= s2(first) ** 2
    C = cands(math.sqrt(max(T1, 0)) + 1e-9, math.sqrt(max(T2, 0)) + 1e-9)
    c1 = C[:, 0] + C[:, 1] * S2; c2 = C[:, 0] - C[:, 1] * S2
    sq = np.stack([C[:, 0] ** 2 + 2 * C[:, 1] ** 2, 2 * C[:, 0] * C[:, 1]], axis=1)
    X0s = [first] if first is not None else [tuple(c) for c in C]
    out = []
    for x0 in X0s:
        r0 = (target[0] - (x0[0] ** 2 + 2 * x0[1] ** 2), target[1] - 2 * x0[0] * x0[1])
        R1, R2 = s1(r0), s2(r0)
        if R1 < -1e-9 or R2 < -1e-9:
            continue
        for i in range(len(C)):
            if c1[i] ** 2 > R1 + 1e-9 or c2[i] ** 2 > R2 + 1e-9:
                continue
            p = r0[0] - sq[i, 0] - sq[:, 0]; q = r0[1] - sq[i, 1] - sq[:, 1]
            seen = set()
            for m, x, y in sqrt_ok(p, q):
                for sgn in (1, -1):
                    for j in np.nonzero(m)[0]:
                        key = (j, sgn * x[j], sgn * y[j])
                        if key in seen:
                            continue
                        seen.add(key)
                        out.append([x0, tuple(C[i]), tuple(C[j]), (sgn * x[j], sgn * y[j])])
    if not out:
        return np.zeros((0, 4, 2), dtype=np.int64)
    X = np.array(out, dtype=np.int64)
    return X[inO(X)]


def okmul(x, y):  # arrays (..., 2)
    return np.stack([x[..., 0] * y[..., 0] + 2 * x[..., 1] * y[..., 1], x[..., 0] * y[..., 1] + x[..., 1] * y[..., 0]], -1)


def qmul(A, B):  # Hamilton product, arrays (..., 4, 2)
    a = [A[..., i, :] for i in range(4)]; b = [B[..., i, :] for i in range(4)]
    m = okmul
    r0 = m(a[0], b[0]) - m(a[1], b[1]) - m(a[2], b[2]) - m(a[3], b[3])
    r1 = m(a[0], b[1]) + m(a[1], b[0]) + m(a[2], b[3]) - m(a[3], b[2])
    r2 = m(a[0], b[2]) - m(a[1], b[3]) + m(a[2], b[0]) + m(a[3], b[1])
    r3 = m(a[0], b[3]) + m(a[1], b[2]) - m(a[2], b[1]) + m(a[3], b[0])
    return np.stack([r0, r1, r2, r3], -2)


def section_counts(n, l, S, Gs):
    """Z-set and per-g counts |Y_{g,s}| (all w, and P-primitive w) for s = S/2."""
    nl = mul(n, l)
    target = (4 * nl[0], 4 * nl[1])            # sum of squares of doubled coords = 4 nrd
    Z = reps(target, first=S)
    if len(Z) == 0:
        return Z, np.zeros(len(Gs), int), np.zeros(len(Gs), int)
    # w = z g / l  ->  W = 2w = Z*G / (2 l)
    d = (2 * l[0], 2 * l[1]); Nd = nrm(d); dc = np.array([d[0], -d[1]], dtype=np.int64)
    allc, primc = [], []
    for G in Gs:
        P = qmul(Z, G[None, :, :])               # (nz, 4, 2)
        num = okmul(P, dc[None, None, :])
        ok = np.all(num % Nd == 0, axis=(1, 2))
        W = num // Nd
        ok &= inO(W)
        allc.append(int(ok.sum())); primc.append(int((ok & primitive(W)).sum()))
    return Z, np.array(allc), np.array(primc)


def find_S(v1, bound2):
    """S in O_K with sigma1(S) ~ v1 and |sigma2(S)| <= bound2."""
    out = []
    dmax = int((abs(v1) + bound2) / (2 * S2)) + 2
    for dd in range(-dmax, dmax + 1):
        c = round(v1 - dd * S2)
        if abs(c + dd * S2 - v1) < 1e-6 and abs(c - dd * S2) <= bound2:
            out.append((c, dd))
    return out


def report(tau, l, S, Gs, rho):
    n = n_tau(tau); nl = mul(n, l)
    c = s1(S) / 2 / math.sqrt(s1(nl)); dist = math.sqrt(max(0.0, 2 - 2 * c))
    m4 = (4 * nl[0] - (S[0] ** 2 + 2 * S[1] ** 2), 4 * nl[1] - 2 * S[0] * S[1])   # 4m
    Nm = nrm(m4) / 16.0
    Z, ca, cp = section_counts(n, l, S, Gs)
    zp = int(primitive(Z).sum()) if len(Z) else 0
    print(f"  tau={tau} l={l} N(l)={nrm(l)} S=2s={S} dist(xi_g, Y)={dist:.6g} (rho={rho:.6g}, ratio {dist/rho:.3f}) "
          f"N(m)={Nm:.4g} sqrtN(m)={math.sqrt(Nm):.3g} #Z={len(Z)} #Z_prim={zp} #g={len(Gs)} "
          f"|Y_g| all: max {ca.max()} mean {ca.mean():.2f}; T-count tau: max {cp.max()} mean {cp.mean():.2f} "
          f"| NOTES formula rho*2^(tau/2)/sqrt(H) = {rho * 2 ** (tau / 2) / math.sqrt(nrm(l)):.3g}")
    return ca, cp


def lclasses(Hmax):
    """totally positive l, one per ideal of norm <= Hmax (ratio sigma1/sigma2 in [1/5.83, 5.83))."""
    lam = (1 + S2) ** 2
    out = []
    for a in range(1, int(3 * math.sqrt(Hmax)) + 3):
        for b in range(-a, a + 1):
            x = (a, b); N = nrm(x)
            if N <= 0 or N > Hmax or s2(x) <= 0 or s1(x) <= 0:
                continue
            r = s1(x) / s2(x)
            if 1 / lam <= r < lam:
                out.append(x)
    return out


def known():
    print("tau=21 cluster at I (clus21.txt): l = 1")
    G1 = reps((4, 0))
    print("  #units of norm 1 in O:", len(G1))
    t = 21; n = n_tau(t); v1 = 2 * 0.99998523409434648 * math.sqrt(s1(n))
    Ss = find_S(v1, 2 * math.sqrt(s2(n)))
    print("  S candidates:", Ss)
    report(t, (1, 0), Ss[0], G1, rho_for(16, t))
    print("tau=26 cluster (clus26.txt): g = ((3+2sqrt2)+i+j-k)/2, l = 5+3sqrt2")
    l = (5, 3)
    Gl = reps((20, 12))
    print("  #g of norm l:", len(Gl))
    gi = [k for k in range(len(Gl)) if tuple(map(tuple, Gl[k])) == ((3, 2), (1, 0), (1, 0), (-1, 0))]
    for t, cval in ((25, 0.9999974727403747), (26, 0.9999967820005790)):
        n = n_tau(t); nl = mul(n, l)
        Ss = find_S(2 * cval * math.sqrt(s1(nl)), 2 * math.sqrt(s2(nl)))
        print("  tau", t, "S candidates:", Ss)
        ca, cp = report(t, l, Ss[0], Gl, rho_for(16, 26))
        if gi:
            print(f"    at the g of the note: |Y_g| all {ca[gi[0]]}, T-count {t}: {cp[gi[0]]}")


def census(Hmax, mu, t0, t1):
    Ls = lclasses(Hmax)
    Gcache = {}
    for t in range(t0, t1 + 1):
        rho = rho_for(mu, t); n = n_tau(t)
        cmin = 1 - rho * rho / 2
        pred = 0.0; ev = []
        for l in Ls:
            nl = mul(n, l)
            pred += S2 * rho * rho * math.sqrt(nrm(nl))
            A1, A2 = 2 * math.sqrt(s1(nl)), 2 * math.sqrt(s2(nl))
            for S in map(tuple, cands(A1, A2)):
                if s1(S) >= A1 * cmin and s1(S) < A1 and abs(s2(S)) < A2:
                    ev.append((l, S))
        print(f"tau={t} mu={mu} rho={rho:.5g} Hmax={Hmax}: {len(ev)} events (heuristic expectation {pred:.2f})")
        for l, S in ev:
            if l not in Gcache:
                Gcache[l] = reps((4 * l[0], 4 * l[1]))
            report(t, l, S, Gcache[l], rho)


if __name__ == "__main__":
    if sys.argv[1] == "known":
        known()
    else:
        census(int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
