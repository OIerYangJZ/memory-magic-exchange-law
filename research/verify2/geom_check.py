"""Numerical stress tests of Lemma E2 (|J| bound, <=2 intervals) and Lemma E3 (3-volume bound) in R^4.
Units: R' = 1. Core circle C(s) = (cos s, sin s, 0, 0); transverse offsets are coords 3,4."""
import numpy as np
rng = np.random.default_rng(1)

def sphere_section(p0, r, rng):
    """random hyperplane H with Y=S^3 cap H of radius r passing through p0 (|p0|=1)."""
    c = np.sqrt(1 - r*r)
    # nu = c p0 + r omega, omega unit perpendicular to p0  => <nu,p0> = c
    w = rng.normal(size=4); w -= w.dot(p0)*p0; w /= np.linalg.norm(w)
    nu = c*p0 + r*w
    return nu, c

def sample_Y(nu, c, r, center_pt, spread, n, rng):
    ctr = c*nu
    # orthonormal basis of nu^perp
    M = np.linalg.svd(np.eye(4) - np.outer(nu, nu))[0][:, :3]
    # point center_pt on Y: direction
    d0 = (center_pt - ctr) / r
    v = d0 @ M
    pts = []
    for _ in range(n):
        u = v + rng.normal(size=3)*spread/r
        u /= np.linalg.norm(u)
        pts.append(ctr + r*(M @ u))
    return np.array(pts)

def vol3(d):
    G = d @ d.T
    return np.sqrt(max(np.linalg.det(G), 0))

# ---- Lemma E3 ----
worst = 0
for trial in range(4000):
    XL = rng.uniform(-3, -0.3); XS = rng.uniform(max(2*XL - 0.0, XL - 2.5), XL)  # 2XL-1<=XS with R'=... use scale
    Rp = 1.0
    xL = 10**XL; xS = max(10**XS, xL*xL/Rp)
    if xS > xL: continue
    r = 10**rng.uniform(np.log10(2*xL), 0) if 2*xL < 0.5 else 0.5
    if 2*xL > r: continue
    s0 = rng.uniform(0, 1); y0 = rng.uniform(-0.01, 0.01, 2)
    def incell(x):
        s = np.arctan2(x[:, 1], x[:, 0])
        return (s >= s0) & (s <= s0 + xL) & (x[:, 2] >= y0[0]) & (x[:, 2] <= y0[0]+xS) & (x[:, 3] >= y0[1]) & (x[:, 3] <= y0[1]+xS)
    # a point in the cell on unit sphere
    sm = s0 + rng.uniform(0, xL); ym = y0 + rng.uniform(0, xS, 2)
    cp = np.sqrt(1 - ym.dot(ym))
    p0 = np.array([cp*np.cos(sm), cp*np.sin(sm), ym[0], ym[1]])
    nu, c = sphere_section(p0, r, rng)
    P = sample_Y(nu, c, r, p0, 2*xL, 6000, rng)
    P = P[incell(P)]
    if len(P) < 4: continue
    best = 0
    for _ in range(300):
        idx = rng.choice(len(P), 4, replace=False)
        best = max(best, vol3(P[idx[1:]] - P[idx[0]]))
    bound = xL*xS*min(xS, xL*xL/r)
    worst = max(worst, best/bound)
print('E3: max observed V / (xL xS min(xS, xL^2/r)) =', worst)

# ---- Lemma E2 ----
worstJ = 0; maxint = 0
for trial in range(3000):
    rho = 10**rng.uniform(-4, -1.5); eps = rho*rng.uniform(1, 3); lam = 10**rng.uniform(-2.5, -0.3)
    r = 10**rng.uniform(np.log10(rho), np.log10(0.5))
    y0 = rng.uniform(-eps, eps, 2)
    sm = rng.uniform(0, lam); ym = y0 + rng.uniform(0, rho, 2)
    cp = np.sqrt(1 - ym.dot(ym)); p0 = np.array([cp*np.cos(sm), cp*np.sin(sm), ym[0], ym[1]])
    nu, c = sphere_section(p0, r, rng)
    # scan core parameters on a fine grid; for each s, check if Y meets the transverse square fibre
    ss = np.linspace(0, lam, 4001); hit = np.zeros_like(ss, bool)
    for k, s in enumerate(ss):
        # fibre: x = (sqrt(1-|y|^2) C(s), y), y in square [y0,y0+rho]^2; <nu,x> = c ?
        g = np.linspace(0, rho, 41)
        Y1, Y2 = np.meshgrid(y0[0]+g, y0[1]+g)
        cp = np.sqrt(1 - Y1**2 - Y2**2)
        f = nu[0]*cp*np.cos(s) + nu[1]*cp*np.sin(s) + nu[2]*Y1 + nu[3]*Y2 - c
        hit[k] = (f.min() <= 0 <= f.max())
    if not hit.any(): continue
    J = hit.mean()*lam
    nint = np.sum(np.diff(hit.astype(int)) == 1) + hit[0]
    maxint = max(maxint, nint)
    bound = min(lam, np.sqrt(rho*(r+rho)), r)
    worstJ = max(worstJ, J/bound)
print('E2: max observed |J| / min(lam, sqrt(rho(r+rho)), r) =', worstJ, '; max #intervals =', maxint)
