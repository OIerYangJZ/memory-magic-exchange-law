#!/usr/bin/env python3
"""Sphere-section clusters at large T-count, found arithmetically (no word enumeration).

Setting (NOTES.md sec. 3): n = 2^{tau/2} (tau even), so the words of T-count tau, tau-2, ... are the w in O with
nrd(w) = n.  For l in O_K^+ and s in (1/2)O_K put m = n l - s^2.  For every g in O with nrd(g) = l the section
    Y_g = { w in O : nrd(w) = n, <w, g> = s }
lies in the dproj-ball of radius d(l,s) = sqrt(2 - 2 sqrt(1 - sigma1(m)/sigma1(nl))) around g/|g| (at sigma1).
We look for (l, s) with d(l,s) <= rho, where rho is the radius at which the Haar expectation of #Lambda_tau in a
ball is mu (mu = 1: the mean spacing).

Counting.  z = w gbar maps pairs (w, g) with <w,g> = s bijectively, up to the 48 norm-one units (w,g) -> (w e, g e),
onto the P-primitive z in O with nrd(z) = n l and Re(z) = s (non-primitive z only add pairs).  Hence
    sum_g |Y_g| >= 48 * #Z_prim,    max_g |Y_g| >= 48 * #Z_prim / r_O(l),
with r_O(l) = #{g in O : nrd g = l}.  The right-hand side is a lower bound for the largest ball count.

Coordinates: X = 2x in O_K^4 ("doubled"), O_K elements as integer pairs (a, b) = a + b sqrt2.
x in O  <=>  X1 - X3 in sqrt2 O_K,  X2 - X3 in sqrt2 O_K,  X0 - X1 - X2 + X3 in 2 O_K     (see NOTES / eq. (order)).

Usage: clusters_arith.py tau Hmax [mu ...]
"""
import sys, math
import numpy as np

S2 = math.sqrt(2.0)
U = 1 + S2


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


def in_O(X0, X1, X2, X3):
    """X_i = (a_i, b_i) integer arrays (doubled coordinates). Vectorized membership test for X/2 in O."""
    c1 = ((X1[0] - X3[0]) % 2 == 0)
    c2 = ((X2[0] - X3[0]) % 2 == 0)
    t0 = X0[0] - X1[0] - X2[0] + X3[0]
    t1 = X0[1] - X1[1] - X2[1] + X3[1]
    return c1 & c2 & (t0 % 2 == 0) & (t1 % 2 == 0)


def div_sqrt2(X):
    a, b = X
    return (b, a // 2)  # valid when a is even: (a + b sqrt2)/sqrt2 = b + (a/2) sqrt2


def primitive(X0, X1, X2, X3):
    """z = X/2 is NOT in sqrt2 * O."""
    ev = (X0[0] % 2 == 0) & (X1[0] % 2 == 0) & (X2[0] % 2 == 0) & (X3[0] % 2 == 0)
    Y = [div_sqrt2(X) for X in (X0, X1, X2, X3)]
    return ~(ev & in_O(*Y))


def cands(r1, r2):
    """All X = (a,b) in O_K with |a + b sqrt2| <= r1, |a - b sqrt2| <= r2 (r1 small)."""
    bmax = int(math.floor((r1 + r2) / (2 * S2))) + 1
    b = np.arange(-bmax, bmax + 1, dtype=np.int64)
    lo = np.ceil(-r1 - b * S2).astype(np.int64)
    hi = np.floor(r1 - b * S2).astype(np.int64)
    A, B = [], []
    for k in range(int(max(0, (hi - lo).max())) + 1):
        a = lo + k
        ok = a <= hi
        A.append(a[ok]); B.append(b[ok])
    a = np.concatenate(A); b = np.concatenate(B)
    s1 = a + b * S2; s2 = a - b * S2
    ok = (np.abs(s1) <= r1) & (np.abs(s2) <= r2)
    return a[ok], b[ok]


def sqrt_OK(p, q):
    """Vectorized square roots in O_K of T = p + q sqrt2.  The roots are {R, -R} (or {0}); returns [(good, x, y)]
    with one root R = x + y sqrt2 per entry where good; callers add -R themselves (skipping R = 0)."""
    s1 = p + q * S2; s2 = p - q * S2
    ok = (s1 >= -0.5) & (s2 >= -0.5)
    r1 = np.sqrt(np.maximum(s1, 0)); r2 = np.sqrt(np.maximum(s2, 0))
    good = np.zeros(np.shape(p), dtype=bool)
    X = np.zeros(np.shape(p), dtype=np.int64); Y = np.zeros(np.shape(p), dtype=np.int64)
    for sg in (1, -1):
        for dx in (0, -1, 1):          # guard against rounding at the boundary
            x = np.rint((r1 + sg * r2) / 2).astype(np.int64) + dx
            y = np.rint((r1 - sg * r2) / (2 * S2)).astype(np.int64)
            g = ok & ~good & (x * x + 2 * y * y == p) & (2 * x * y == q)
            X = np.where(g, x, X); Y = np.where(g, y, Y); good |= g
    return [(good, X, Y)]


def r_O(la, lb):
    """#{g in O : nrd g = l}, l = la + lb sqrt2, by enumeration in doubled coordinates (sum X_i^2 = 4l)."""
    L1 = 4 * (la + lb * S2); L2 = 4 * (la - lb * S2)
    a, b = cands(math.sqrt(L1) + 1e-9, math.sqrt(L2) + 1e-9)
    s1 = a + b * S2; s2 = a - b * S2
    cnt = 0
    n = len(a)
    for i in range(n):
        for j in range(n):
            if s1[i] ** 2 + s1[j] ** 2 > L1 + 1e-9 or s2[i] ** 2 + s2[j] ** 2 > L2 + 1e-9:
                continue
            # X0 = (a_i,b_i), X1 = (a_j,b_j); X2 over all candidates, X3 solved
            p = 4 * la - (a[i] ** 2 + 2 * b[i] ** 2) - (a[j] ** 2 + 2 * b[j] ** 2) - (a * a + 2 * b * b)
            q = 4 * lb - 2 * a[i] * b[i] - 2 * a[j] * b[j] - 2 * a * b
            for good, x, y in sqrt_OK(p, q):
                # both sign variants in the loop give X3 and its conjugate-sign partner; +-X3 counted via x,y and -x,-y
                for sgn in (1, -1):
                    X0 = (np.full(n, a[i]), np.full(n, b[i])); X1 = (np.full(n, a[j]), np.full(n, b[j]))
                    X2 = (a, b); X3 = (sgn * x, sgn * y)
                    m = good & in_O(X0, X1, X2, X3)
                    # avoid double counting when x = y = 0 (X3 = -X3) or when the two sigma2-sign branches coincide
                    if sgn == -1:
                        m &= ~((x == 0) & (y == 0))
                    cnt += int(m.sum())
    # the two sigma2-sign branches coincide exactly when sigma2-root = 0, i.e. p - q sqrt2 = 0: then only one root pair
    return cnt


def mul_OK(x, y):
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def divides(d, l):
    """exact division l / d in Z[sqrt2]; returns quotient or None."""
    Nd = d[0] * d[0] - 2 * d[1] * d[1]
    num = mul_OK(l, (d[0], -d[1]))
    if num[0] % Nd or num[1] % Nd:
        return None
    return (num[0] // Nd, num[1] // Nd)


_PRIMES = None


def primes_OK(limit):
    """One generator per prime ideal of Z[sqrt2] with norm <= limit."""
    out = [(0, 1)]  # sqrt2, norm 2
    for p in range(3, limit + 1, 2):
        if any(p % q == 0 for q in range(3, int(p ** 0.5) + 1, 2)):
            continue
        if p % 8 in (1, 7):
            # find a^2 - 2 b^2 = +-p
            found = None
            for b in range(0, p):
                for sgn in (1, -1):
                    a2 = sgn * p + 2 * b * b
                    if a2 >= 0:
                        a = int(round(a2 ** 0.5))
                        if a * a == a2:
                            found = (a, b); break
                if found:
                    break
            out.append(found); out.append((found[0], -found[1]))
        elif p * p <= limit:
            out.append((p, 0))
    return out


def sigma_ideal(la, lb, limit):
    """sum over ideals d | (l) of N(d), by trial division."""
    global _PRIMES
    if _PRIMES is None:
        _PRIMES = primes_OK(limit)
    l = (la, lb); tot = 1
    for pi in _PRIMES:
        Np = abs(pi[0] * pi[0] - 2 * pi[1] * pi[1])
        e = 0
        while True:
            q = divides(pi, l)
            if q is None:
                break
            l = q; e += 1
        if e:
            tot *= sum(Np ** k for k in range(e + 1))
    assert abs(l[0] * l[0] - 2 * l[1] * l[1]) == 1, (la, lb, l)
    return tot


def count_Z(tau, la, lb, S, r1sq, r2sq):
    """#Z_prim: z = X/2 with X0 = S, X1^2+X2^2+X3^2 = 4m (m = nl - S^2/4), X/2 in O, primitive."""
    n = 2 ** (tau // 2)
    # 4m = 4 n l - S^2
    Sa, Sb = S
    pm = 4 * n * la - (Sa * Sa + 2 * Sb * Sb)
    qm = 4 * n * lb - 2 * Sa * Sb
    r1 = math.sqrt(max(r1sq, 0)) + 1e-12; r2 = math.sqrt(r2sq) + 1e-9
    a, b = cands(r1, r2)
    s1 = a + b * S2; s2 = a - b * S2
    total = 0
    X0 = None
    for i in range(len(a)):
        rest1 = r1sq - s1[i] ** 2; rest2 = r2sq - s2[i] ** 2
        sel = (s1 ** 2 <= rest1 + 1e-12) & (s2 ** 2 <= rest2 + 1e-6)
        if not sel.any():
            continue
        a2 = a[sel]; b2 = b[sel]; k = len(a2)
        p = pm - (a[i] ** 2 + 2 * b[i] ** 2) - (a2 * a2 + 2 * b2 * b2)
        q = qm - 2 * a[i] * b[i] - 2 * a2 * b2
        seen = np.zeros(k, dtype=bool)
        for good, x, y in sqrt_OK(p, q):
            for sgn in (1, -1):
                X0 = (np.full(k, Sa, dtype=np.int64), np.full(k, Sb, dtype=np.int64))
                X1 = (np.full(k, a[i], dtype=np.int64), np.full(k, b[i], dtype=np.int64))
                X2 = (a2, b2); X3 = (sgn * x, sgn * y)
                m = good & in_O(X0, X1, X2, X3) & primitive(X0, X1, X2, X3)
                if sgn == -1:
                    m &= ~((x == 0) & (y == 0))
                total += int(m.sum())
        del seen
    return total


def check_rO():
    """r_O(l) = 48 * sum_{d | l} N(d): brute force for small l."""
    for la, lb in [(1, 0), (2, 1), (3, 1), (3, 0), (4, 1), (5, 2), (5, 3), (7, 3), (2, 0), (4, 2), (6, 1), (9, 4)]:
        if la <= abs(lb) * S2:
            continue
        bf = r_O(la, lb); fo = 48 * sigma_ideal(la, lb, 200)
        print(f"check r_O l=({la},{lb}) N={la*la-2*lb*lb}: brute {bf} formula {fo}", flush=True)


def main():
    if sys.argv[1] == "check":
        check_rO(); return
    tau = int(sys.argv[1]); Hmax = int(sys.argv[2])
    mus = [float(x) for x in sys.argv[3:]] or [1.0]
    assert tau % 2 == 0
    n = 2 ** (tau // 2)
    rhos = {mu: rho_for(mu, tau) for mu in mus}
    rho = max(rhos.values())
    print(f"# tau={tau} n=2^{tau//2} Hmax={Hmax} rhos={rhos}", flush=True)
    # l = la + lb sqrt2 totally positive, N(l) <= Hmax, one per unit-square class: ratio sigma1/sigma2 in [U^-2, U^2)
    ls = []
    B = int(math.sqrt(Hmax) * U ** 2) + 2
    for lb in range(-B, B + 1):
        for la in range(1, 4 * B + 2):
            N = la * la - 2 * lb * lb
            if N < 1 or N > Hmax or la <= abs(lb) * S2:
                continue
            r = (la + lb * S2) / (la - lb * S2)
            if U ** -2 <= r < U ** 2:
                ls.append((la, lb, N))
    ls.sort(key=lambda t: t[2])
    print(f"# {len(ls)} classes of l", flush=True)
    events = []
    for la, lb, N in ls:
        A1 = 4 * n * (la + lb * S2); A2 = 4 * n * (la - lb * S2)   # sigma_i(4 n l)
        # S = c + d sqrt2 with A1 (1 - rho^2) <= sigma1(S)^2 < A1 and sigma2(S)^2 < A2
        hi1 = math.sqrt(A1); lo1 = math.sqrt(A1 * (1 - rho * rho * (1 - rho * rho / 4)))
        # d range from sigma1 - sigma2 = 2 d sqrt2 with |sigma2| < sqrt(A2)
        dlo = int(math.floor((lo1 - math.sqrt(A2)) / (2 * S2))) - 1
        dhi = int(math.ceil((hi1 + math.sqrt(A2)) / (2 * S2))) + 1
        d = np.arange(dlo, dhi + 1, dtype=np.int64)
        c = np.floor(hi1 - d * S2).astype(np.int64)
        s1 = c + d * S2
        ok = (s1 >= lo1 - 1e-9) & (np.abs(c - d * S2) < math.sqrt(A2))
        for cc, dd in zip(c[ok], d[ok]):
            cc = int(cc); dd = int(dd)
            # exact: 4m = 4nl - S^2 = pm + qm sqrt2, need sigma1 > 0, sigma2 > 0
            pm = 4 * n * la - (cc * cc + 2 * dd * dd); qm = 4 * n * lb - 2 * cc * dd
            m1 = pm + qm * S2; m2 = pm - qm * S2
            if not (m1 > 0 and m2 > 0):
                continue
            ratio = m1 / (4 * n * (la + lb * S2))       # sigma1(m)/sigma1(nl)
            dist = math.sqrt(2 - 2 * math.sqrt(1 - ratio))
            if dist > rho:
                continue
            events.append((N, la, lb, cc, dd, m1, m2, dist))
    print(f"# {len(events)} events with d <= {rho:.4g}", flush=True)
    best = {mu: (0, None) for mu in mus}
    rO_cache = {}
    for N, la, lb, cc, dd, m1, m2, dist in events:
        nz = count_Z(tau, la, lb, (cc, dd), m1, m2)
        key = (la, lb)
        if key not in rO_cache:
            rO_cache[key] = 48 * sigma_ideal(la, lb, Hmax)
        ro = rO_cache[key]
        avg = 48 * nz / ro if ro else float('nan')
        print(f"N(l)={N} l=({la},{lb}) S=({cc},{dd}) dist={dist:.4g} N(m)~{m1*m2/16:.3g} #Z_prim={nz} r_O(l)={ro} "
              f"avg|Y_g|>={avg:.2f} pred~{math.sqrt(m1*m2)/4/N * 1:.1f}", flush=True)
        for mu in mus:
            if dist <= rhos[mu] and avg > best[mu][0]:
                best[mu] = (avg, (N, la, lb, cc, dd))
    for mu in mus:
        print(f"RESULT tau={tau} mu={mu} rho={rhos[mu]:.4g} best_avg_cluster={best[mu][0]:.2f} at {best[mu][1]}")


if __name__ == "__main__":
    main()
