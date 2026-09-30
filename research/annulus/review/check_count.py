"""Referee check of the lattice-counting steps on small random instances.

A1 step 4:  Lambda = P_{Vperp}(O) (rank-4 Z-lattice), D = {sigma1 x in [-l,l]x[-4d,4d] (random frame of
            V1perp), |sigma2 x| <= 4R}.  Checks: lambda2 <= 2.4143 lambda1, lambda4 <= 2.4143 lambda3
            (successive minima w.r.t. D), and #(Lambda cap D) <= C (1 + H_V l R + H_V l d R^2).
P3 step 5:  L_B = Lambda_W cap n_B-perp, D = {sigma1 x in random rectangle eb x eu, |sigma2 x| <= 4R}.
            Checks: if mu2 <= 1 then # <= C vol/covol;  if mu2 > 1 then all points of L_B cap D lie in one
            K-line (exact test);  and # <= C(1 + eb eu R^2 h(m)/h(n_B)).
Lattices are exact (klat); the enumeration is float with a safety factor."""
import sys, random
import numpy as np
from fractions import Fraction as F
from klat import *
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
nrng = np.random.default_rng(3)
def hsq(v): return kdot(v, v).norm()
def hsq_wedge(vs):
    w = wedge(vs); return sum((c * c for c in w), Kel(0)).norm()
def fsqrt(x): return mp.sqrt(mp.mpf(x.numerator) / x.denominator)
def primitive_in(L_basis, n8):
    line = lattice_cap_subspace(L_basis, [kv(n8)])
    return covol2(line) == covol2([n8, q8(kscale(R2, kv(n8)))])
def rand_primitive_Ostar_in(subK, tries=200):
    Lb = lattice_cap_subspace(OSTAR_Z, subK)
    for _ in range(tries):
        n8 = rand_lattice_vec(rng, Lb, 2)
        if all(c == 0 for c in n8): continue
        if fn_index(n8, O_Z) == 1: return n8
    return None

def frame_sigma(Lb, which):
    """orthonormal basis (4 x 2) of the real span of sigma_which(L)"""
    E = np.array([emb8(x) for x in Lb]); X = E[:, :4] if which == 1 else E[:, 4:]
    U, s, _ = np.linalg.svd(X.T); return U[:, :2]

def body_count(Lb, e1, e2, A, Bw, R2r):
    """points of L in {|<s1x,e1>|<=A, |<s1x,e2>|<=Bw, |s2x|<=R2r}; returns list of (z, x8float), plus norm fn"""
    E = np.array([emb8(x) for x in Lb])
    f1 = E[:, :4] @ e1 / A; f2 = E[:, :4] @ e2 / Bw; G = E[:, 4:] / R2r
    M = np.outer(f1, f1) + np.outer(f2, f2) + G @ G.T
    def ND(z):
        z = np.array(z); return max(abs(z @ f1), abs(z @ f2), np.linalg.norm(z @ G))
    return M, ND

def succ_minima(Lb, M, ND, tmax=64.0):
    t = 1.0
    while True:
        Z = enum_box(M, 3 * t * t)
        if Z is None: return None, None
        pts = [tuple(z) for z in Z if any(z)]
        pts.sort(key=ND)
        chosen = []
        for z in pts:
            if ND(z) > t: break
            if np.linalg.matrix_rank(np.array(chosen + [z], dtype=float)) == len(chosen) + 1:
                chosen.append(z)
                if len(chosen) == 4: return [ND(c) for c in chosen], chosen
        t *= 2
        if t > tmax: return None, None

def count_in(M, ND):
    Z = enum_box(M, 3.0)
    if Z is None: return None
    return [tuple(z) for z in Z if ND(z) <= 1 + 1e-12]

def same_Kline(Lb, zs):
    vs = [kv(combo(z, Lb)) for z in zs if any(z)]
    if len(vs) <= 1: return True
    return all(all(c.iszero() for c in wedge([vs[0], v])) for v in vs[1:])

# ------------------------------------------------------------ A1 step 4
print("== A1 step 4 ==")
rat = []; mins_ok = True; big3 = []
for trial in range(30):
    v1, v2 = rand_vec(rng, rng.choice([1, 2, 3])), rand_vec(rng, 2)
    if hsq_wedge([v1, v2]) == 0: continue
    V = [v1, v2]; Vp = kperp(V)
    S = lattice_cap_subspace(OK4_Z, V); HV = float(fsqrt(covol2(S) / 64))
    Lam = reduce_exact(lattice_proj(O_Z, Vp))
    U1 = frame_sigma(Lam, 1)
    for rep in range(3):
        ang = nrng.uniform(0, np.pi); e1 = U1 @ [np.cos(ang), np.sin(ang)]; e2 = U1 @ [-np.sin(ang), np.cos(ang)]
        Rr = float(10 ** nrng.uniform(-0.3, 1.0)); d1 = float(10 ** nrng.uniform(-2.5, -0.3)) ; l = d1 * float(10 ** nrng.uniform(0.6, 1.3))
        # D = B - B for the rectangle l x 4 d1 (half-sides l, 4 d1) and |s2| <= 2R  ->  |s2| <= 4R
        M, ND = body_count(Lam, e1, e2, l, 4 * d1, 4 * Rr)
        lam, _ = succ_minima(Lam, M, ND)
        if lam is not None:
            if not (lam[1] <= 2.4143 * lam[0] + 1e-9 and lam[3] <= 2.4143 * lam[2] + 1e-9): mins_ok = False
        pts = count_in(M, ND)
        if pts is None: continue
        bound = 1 + HV * l * Rr + HV * l * d1 * Rr ** 2
        volcov = (2 * l) * (8 * d1) * np.pi * (4 * Rr) ** 2 / float(covol(Lam))
        rat.append((len(pts) / bound, len(pts) / (1 + HV * l * Rr + volcov)))
        if lam is not None and lam[2] > 1:
            big3.append((len(pts) / max(1, lam[0] ** -2), lam[0] ** -2 / max(1e-9, HV * l * Rr)))
print(f"  successive minima: lambda2<=2.4143 lambda1 and lambda4<=2.4143 lambda3 on all instances: {mins_ok}")
rr = np.array(rat)
print(f"  #(Lambda cap D) / (1 + H_V l R + H_V l d R^2): max {rr[:,0].max():.3f};  #/(1 + H_V l R + vol(D)/covol): max {rr[:,1].max():.3f}  over {len(rr)} instances")
if big3:
    b = np.array(big3)
    print(f"  lambda3>1 cases ({len(b)}): #/max(1,lambda1^-2) <= {b[:,0].max():.2f};  lambda1^-2/(H_V l R) <= {b[:,1].max():.3f}")

# ------------------------------------------------------------ P3 step 5
print("== P3 step 5 ==")
rat = []; line_ok = True; nline = 0; nlat = 0; ratvol = []
done = 0
while done < 25:
    m = rand_vec(rng, rng.choice([1, 2, 3]))
    if all(c.iszero() for c in m) or not primitive_in(OK4_Z, q8(m)): continue
    W = kperp([m]); LW = lattice_proj(O_Z, W)
    nB = rand_primitive_Ostar_in(W)
    if nB is None: continue
    done += 1
    comp = qspan_of([kv(nB), m])
    MM = [[trace_form(x, c) for c in comp] for x in LW]
    LB = reduce_exact([combo(z, LW) for z in int_kernel(MM)])
    cv = float(covol(LB)); hm = float(fsqrt(hsq(m))); hn = float(fsqrt(hsq(kv(nB))))
    U1 = frame_sigma(LB, 1)
    for rep in range(3):
        ang = nrng.uniform(0, np.pi); e1 = U1 @ [np.cos(ang), np.sin(ang)]; e2 = U1 @ [-np.sin(ang), np.cos(ang)]
        Rr = float(10 ** nrng.uniform(-0.3, 1.0)); eu = float(10 ** nrng.uniform(-1.5, 0.2)); eb = eu * float(10 ** nrng.uniform(0, 1.0))
        M, ND = body_count(LB, e1, e2, eb, eu, 4 * Rr)
        lam, ch = succ_minima(LB, M, ND)
        pts = count_in(M, ND)
        if pts is None: continue
        vol = (2 * eb) * (2 * eu) * np.pi * (4 * Rr) ** 2
        rat.append(len(pts) / (1 + vol / cv))
        # K-minima: mu1 = lambda1; mu2 <= 1 iff two K-independent vectors in D
        nz = [z for z in pts if any(z)]
        indep = not same_Kline(LB, nz)
        if indep:
            nlat += 1; ratvol.append(len(pts) / (vol / cv))
        else:
            nline += 1
print(f"  #(L_B cap D) / (1 + vol(D)/covol(L_B)): max {max(rat):.3f} over {len(rat)} instances  (covol(L_B) ~ h(n_B)/h(m) checked exactly in check_covol.py)")
print(f"  instances with two K-independent vectors in D (mu2<=1): {nlat}, #/(vol/covol) <= {max(ratvol) if ratvol else float('nan'):.3f}")
print(f"  instances with all of L_B cap D in one K-line (mu2>1): {nline}  (then the u's lie on an affine K-line -> circle, E5)")
