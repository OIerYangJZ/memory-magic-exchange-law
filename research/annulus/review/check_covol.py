"""Referee checks (exact over Q(sqrt2)) of the covolume identities and height bounds used in
A1 (step 5), P2 (step 3), P3 (step 4) and the rank-3 route N (Hodge star), on the explicit order O.

All covolume^2 are exact rationals (trace-form Gram determinants); heights h(v)^2 = N_{K/Q}<v,v>
are exact too.  Floats (mpmath, 40 digits) only for printing ratios and for short-vector search."""
import sys, random
from fractions import Fraction as F
import numpy as np
import mpmath as mp
from klat import *

def hsq(v):             # h(v)^2 = sigma1<v,v> sigma2<v,v> = N<v,v>, exact
    return kdot(v, v).norm()
def hsq_wedge(vs):      # h(v1^...^vk)^2 = N(sum of squares of Pluecker coords), exact
    w = wedge(vs)
    return kdot(w, w).norm() if len(w) == 4 else sum((c * c for c in w), Kel(0)).norm()
def fsqrt(x): return mp.sqrt(mp.mpf(x.numerator) / x.denominator)

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 25

# ---------------------------------------------------------------- basic facts
assert covol2(O_Z) == 16, covol2(O_Z)
c2s = covol2(OSTAR_Z)
print("covol(O)^2 =", covol2(O_Z), " covol(O*)^2 =", c2s, "-> covol(O*) =", fsqrt(c2s))
print("O_K^4 in O:", all(in_lattice(e, O_Z) for e in OK4_Z),
      " O* in O_K^4:", all(all(c.denominator == 1 for c in x) for x in OSTAR_Z),
      " 2O in O*:", all(in_lattice(tuple(2 * c for c in x), OSTAR_Z) for x in O_Z),
      " O* in 2O:", all(in_lattice(tuple(c / 2 for c in x), O_Z) for x in OSTAR_Z))

def primitive_in(L_basis, n8):
    """is n primitive in the lattice L (i.e. Kn cap L = O_K n)?  compute Kn cap L."""
    line = lattice_cap_subspace(L_basis, [kv(n8)])
    return covol2(line) == covol2([n8, q8(kscale(R2, kv(n8)))])

def rand_primitive_Ostar_in(subK, tries=200):
    """random primitive n in O* lying in span_K(subK)"""
    Lb = lattice_cap_subspace(OSTAR_Z, subK)
    for _ in range(tries):
        n8 = rand_lattice_vec(rng, Lb, 2)
        if all(c == 0 for c in n8): continue
        if fn_index(n8, O_Z) == 1:          # F1: primitive <=> <n,O> = O_K
            return n8
    return None

def minh(Lb, H):
    """min of h over nonzero vectors of lattice Lb (Q^8 rows), searched among |x|^2 <= 2.83*H*1.2"""
    E = np.array([emb8(x) for x in Lb]); B = lll(E)
    M = B @ B.T
    pts = enum_quadform(M, 2.9 * H)
    best = None
    for z in pts:
        if not any(z): continue
        x = np.array(z) @ B
        hv = np.linalg.norm(x[:4]) * np.linalg.norm(x[4:])
        if best is None or hv < best: best = hv
    return best

# ---------------------------------------------------------------- planes V (A1, P2)
print("\n== K-planes V: A1 step 5, P2 step 3 ==")
worst = dict(lamV=[1e9, -1e9], lamVp=[1e9, -1e9], hmin=[1e9, -1e9], KB=[1e9, -1e9])
for t in range(NT):
    B = rng.choice([2, 3, 5])
    v1, v2 = rand_vec(rng, B), rand_vec(rng, B)
    if hsq_wedge([v1, v2]) == 0: continue
    V = [v1, v2]; Vp = kperp(V)
    S = lattice_cap_subspace(OK4_Z, V)                    # saturated O_K^4 cap V
    assert len(S) == 4
    HV2 = covol2(S) / 64                                   # paper: covol(S) = 8 H_V
    # independent computation of H_V = h(s1^s2): h(v1^v2)/[S : O_K v1 + O_K v2]
    sub = [q8(v1), q8(kscale(R2, v1)), q8(v2), q8(kscale(R2, v2))]
    idx2 = covol2(sub) / covol2(S)                         # index^2
    HV2_alt = hsq_wedge(V) / idx2
    assert HV2 == HV2_alt, (HV2, HV2_alt)
    HV = fsqrt(HV2)
    # O cap V, P_{Vperp}(O), product = covol(O)
    OV = lattice_cap_subspace(O_Z, V); Lam = lattice_proj(O_Z, Vp)
    assert len(OV) == 4 and len(Lam) == 4
    assert covol2(OV) * covol2(Lam) == 16, "covol(O cap V) covol(P_Vperp O) != covol O"
    r = covol(Lam) * HV
    worst['lamVp'] = [min(worst['lamVp'][0], r), max(worst['lamVp'][1], r)]
    # P_V(O) (for P2)
    OVp = lattice_cap_subspace(O_Z, Vp); LamV = lattice_proj(O_Z, V)
    assert covol2(OVp) * covol2(LamV) == 16
    Sp = lattice_cap_subspace(OK4_Z, Vp)
    assert covol2(Sp) == covol2(S), "height of V-perp != height of V"
    r2 = covol(LamV) * HV
    worst['lamV'] = [min(worst['lamV'][0], r2), max(worst['lamV'][1], r2)]
    # A1 step 5b: min h over P_{Vperp}(O) \ 0  >=  1/(C H_V)
    hm = minh(Lam, float(8 / HV) + 4)
    if hm is not None:
        worst['hmin'] = [min(worst['hmin'][0], hm * HV), max(worst['hmin'][1], hm * HV)]
    # P2 step 3: K_B = P_V(O) cap n_B-perp, n_B primitive in O*, n_B in V
    nB = rand_primitive_Ostar_in(V)
    if nB is not None:
        assert fn_index(nB, LamV) == 1, "f(P_V O) != O_K"
        KB = lattice_cap_subspace(LamV, [kv(nB)] + Vp) if False else None
        # K_B: vectors of P_V(O) orthogonal to n_B
        comp = qspan_of([kv(nB)] + Vp)
        M = [[trace_form(x, c) for c in comp] for x in LamV]
        KB = [combo(z, LamV) for z in int_kernel(M)]
        assert len(KB) == 2
        # exact sequence: covol(K_B)^2 = covol(P_V O)^2 h(n_B)^2 / 8
        assert covol2(KB) == covol2(LamV) * hsq(kv(nB)) / 8, "exact sequence P2 fails"
        rKB = covol(KB) * HV / fsqrt(hsq(kv(nB)))
        worst['KB'] = [min(worst['KB'][0], rKB), max(worst['KB'][1], rKB)]
print("  covol(S)=8H_V (vs h(v1^v2)/index): exact on all trials")
print("  covol(O cap V) covol(P_Vperp O) = covol(O) = 4: exact;  height(V-perp) = height(V): exact")
print("  covol(P_Vperp O) * H_V in [%.4f, %.4f]   (claimed ~1, bounds [1/2, 8])" % tuple(worst['lamVp']))
print("  covol(P_V O)     * H_V in [%.4f, %.4f]" % tuple(worst['lamV']))
print("  min_{x in P_Vperp O} h(x) * H_V in [%.4f, %.4f]   (A1 claims >= 1/C)" % tuple(worst['hmin']))
print("  covol(K_B) H_V / h(n_B) in [%.4f, %.4f];  covol(K_B)^2 = covol(P_V O)^2 h(n_B)^2/8 exact;  f(P_V O)=O_K" % tuple(worst['KB']))

# ---------------------------------------------------------------- lines m (P3, rank-3 route N)
print("\n== K-lines m, W = m-perp: P3 step 4, Hodge star ==")
worst = dict(LW=[1e9, -1e9], LB=[1e9, -1e9], hodge=[1e9, -1e9])
done = 0
while done < NT:
    B = rng.choice([2, 3, 5])
    m = rand_vec(rng, B)
    if all(c.iszero() for c in m): continue
    m8 = q8(m)
    if not primitive_in(OK4_Z, m8): continue
    done += 1
    W = kperp([m])
    # W cap O_K^4 has height h(m): covol = (2 sqrt2)^3 h(m)
    SW = lattice_cap_subspace(OK4_Z, W)
    assert covol2(SW) == 512 * hsq(m), "Hodge star / covol(W cap O_K^4) != 16 sqrt2 h(m)"
    # three K-independent v in W cap O*: h(v1^v2^v3) is |N xi| h(m), xi in O_K \ 0
    WO = lattice_cap_subspace(OSTAR_Z, W)
    vs = [kv(rand_lattice_vec(rng, WO, 2)) for _ in range(3)]
    hw = hsq_wedge(vs)
    if hw != 0:
        ratio = hw / hsq(m)                                 # = N(xi)^2, should be a nonzero square integer
        assert ratio.denominator == 1 and ratio >= 1, ratio
        worst['hodge'] = [min(worst['hodge'][0], ratio), max(worst['hodge'][1], ratio)]
    Om = lattice_cap_subspace(O_Z, [m]); LW = lattice_proj(O_Z, W)
    assert covol2(Om) * covol2(LW) == 16
    r = covol(LW) * fsqrt(hsq(m))
    worst['LW'] = [min(worst['LW'][0], r), max(worst['LW'][1], r)]
    nB = rand_primitive_Ostar_in(W)
    if nB is None: continue
    assert fn_index(nB, LW) == 1, "f(Lambda_W) != O_K"
    comp = qspan_of([kv(nB), m])
    M = [[trace_form(x, c) for c in comp] for x in LW]
    LB = [combo(z, LW) for z in int_kernel(M)]
    assert len(LB) == 4
    assert covol2(LB) == covol2(LW) * hsq(kv(nB)) / 8, "exact sequence P3 fails"
    rr = covol(LB) * fsqrt(hsq(m)) / fsqrt(hsq(kv(nB)))
    worst['LB'] = [min(worst['LB'][0], rr), max(worst['LB'][1], rr)]
print("  covol(W cap O_K^4) = 16 sqrt2 h(m) exact (height of W = h(m));  h(v1^v2^v3)/h(m) = |N xi| >= 1: range of ratio^2", worst['hodge'])
print("  covol(O cap Km) covol(Lambda_W) = 4 exact")
print("  covol(Lambda_W) * h(m) in [%.4f, %.4f]   (claimed ~1)" % tuple(worst['LW']))
print("  covol(L_B) h(m)/h(n_B) in [%.4f, %.4f];  covol(L_B)^2 = covol(Lambda_W)^2 h(n_B)^2/8 exact;  f(Lambda_W)=O_K" % tuple(worst['LB']))
