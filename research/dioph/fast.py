"""Vectorised float optimiser (same model as an earlier optimiser, opt5.py, that was not kept), used to locate good parameters.
Exponents in units of log R'.  See DIOPH_NOTES.md for the lemmas behind every term."""
import numpy as np

T_GRID = np.arange(0.0, 2.2, 0.004)          # XS = top - t
MARGIN = 0.002

def cell_cost_vec(S0, U0, R0, eta, DEL=1e-9):
    """S0,U0,R0: arrays over g.  min_{XL,XS} N + max(0,P)."""
    S0 = S0[:, None]; U0 = U0[:, None]; R0 = R0[:, None]
    top = np.minimum(U0, R0 - DEL)
    XS = top - T_GRID[None, :]
    best = np.full(S0.shape[0], np.inf)
    cands = (S0 - U0 + XS, (XS + R0) / 2, S0 + 0 * XS, XS, (1 + XS) / 2, R0 - DEL + 0 * XS,
             eta - 3 - 2 * XS, (eta - 3 - XS + R0) / 3)
    for XL in cands:
        XL = np.minimum(np.minimum(XL, S0), np.minimum(R0 - DEL, (1 + XS) / 2))
        XL = np.maximum(XL, XS)
        ok = (XL < R0) & (2 * XL - 1 <= XS + 1e-12)
        N = np.maximum(np.maximum(np.maximum(0, S0 - XL), np.maximum(U0 - XS, S0 + U0 - XL - XS)), 2 * (U0 - XS))
        P = (XL + XS + np.minimum(XS, 2 * XL - R0) - eta) / 3 + 1
        c = np.where(ok, N + np.maximum(0, P), np.inf)
        best = np.minimum(best, c.min(axis=1))
    return best

def patch_vec(a, b, g, zero=False):
    S0 = np.full_like(g, 1 - b) if zero else np.minimum(np.minimum(1 - b, (2 - a - np.minimum(g, a)) / 2), 1 - g)
    return S0, np.minimum(1 - a, 1 - g)

def cost_vec(a, b, S0, U0, R0, eta):
    c = cell_cost_vec(S0, U0, R0, eta)
    if 2 * b >= a:
        c = np.minimum(c, np.maximum(0.0, (S0 + 2 * U0 - eta) / 3 + 1))
    return c

def normals(a0, gl, eta, Fstar, dichotomy=True):
    gam = np.minimum(gl, a0)
    G = np.maximum(0.0, gam + Fstar / 2 - 2) if dichotomy else 0.0
    return np.maximum(np.maximum(eta, 2 * eta - G), np.maximum(3 * eta - gam, 4 * eta - 2 * gam))

def low_sphere_vec(a0, a, b, gl, gr, eta_n, eta_c, nD=8):
    lde = np.maximum(1 - gl - a0, 1 - 2 * a0)
    ldb = np.maximum(1 - gl - a, 1 - 2 * a)
    Dmax = np.maximum(1 - 2 * gl, lde)
    UT = np.minimum(1 - a0, 1 - gl); U0 = np.minimum(1 - a, 1 - gl)
    worst = np.full_like(gl, -np.inf)
    # tangent class, then crossing classes D in [lde, Dmax] (intervals: offs at right end, per at left end)
    fr = np.linspace(0, 1, nD + 1)
    pieces = [("tan", None, None)] + [("x", fr[i], fr[i + 1]) for i in range(nD)]
    for kind, f0, f1 in pieces:
        if kind == "tan":
            Doff = lde
            ST = np.minimum((1 + lde) / 2, 1 - gl)
            S0 = np.minimum(np.minimum(1 - b, (1 + ldb) / 2), 1 - gl)
        else:
            Dl = lde + f0 * (Dmax - lde); Dr = lde + f1 * (Dmax - lde)
            Doff = Dr
            ST = np.minimum(np.minimum((1 + lde) / 2, lde + (1 - Dl) / 2), 1 - gl)
            S0 = np.minimum(np.minimum(1 - b, (1 + ldb) / 2), np.minimum(ldb + (1 - Dl) / 2, 1 - gl))
        offs = np.maximum(0, eta_n + Doff + 1)
        boxes = np.maximum(0, ST - (1 - b)) + 2 * np.maximum(0, UT - (1 - a))
        per = np.minimum(boxes + cost_vec(a, b, S0, U0, 1 - gr, eta_c),
                         np.maximum(0.0, (ST + 2 * UT - eta_c) / 3 + 1))
        worst = np.maximum(worst, offs + per)
    return worst

def E(a0, a, b, den=100, gmax=2.02, deta=0.01, etamax=1.6, dichotomy=True, nD=8, split=True, profile=False):
    if not (a >= a0 and a >= b >= 0 and a0 + a >= 2 * b): return np.inf, None
    Fstar = a0 - MARGIN
    k = max(0.0, 2 - a / 2 - 3 * b / 4)
    E1 = 2 * (a - a0) + b + k
    g0 = np.array([0.0])
    S0z, U0z = patch_vec(a, b, g0, zero=True)
    zero = E1 + cost_vec(a, b, S0z, U0z, 1 - np.array([1 / den]), 0.0)[0]
    gs = np.arange(0, gmax + 1e-9, 1 / den)
    gl, gr = gs[:-1], gs[1:]
    S0, U0 = patch_vec(a, b, gl)
    best = E1 + cost_vec(a, b, S0, U0, 1 - gr, 0.0)
    arg = np.zeros_like(gl)
    if split:
        run = np.full_like(gl, -np.inf)
        etas = np.arange(0, etamax + 1e-9, deta)
        alive = np.ones_like(gl, dtype=bool)
        for j in range(1, len(etas)):
            run = np.maximum(run, normals(a0, gl, etas[j], Fstar, dichotomy)
                             + low_sphere_vec(a0, a, b, gl, gr, etas[j], etas[j - 1], nD))
            alive &= run < best
            if not alive.any(): break
            hi = E1 + cost_vec(a, b, S0, U0, 1 - gr, etas[j])
            v = np.maximum(hi, run)
            upd = alive & (v < best) & (gl > 0)
            best = np.where(upd, v, best); arg = np.where(upd, etas[j], arg)
    i = int(np.argmax(best))
    if profile:
        return best, arg, gl
    w = max(zero, best[i])
    return w, (round(gl[i], 3), round(arg[i], 3), round(zero, 4), round(E1, 4))

def best_ab(alpha, a_range=(0, 0.06), b_range=(1.45, 1.62), step=0.01, **kw):
    a0 = 4 / alpha; out = None
    for da in np.arange(a_range[0], a_range[1] + 1e-9, step):
        a = a0 + da
        for b in np.arange(b_range[0], min(a, (a0 + a) / 2, b_range[1]) + 1e-9, step):
            e, w = E(a0, a, b, **kw)
            if out is None or e - a0 < out[0]:
                out = (round(e - a0, 5), round(a, 4), round(b, 4), round(e, 5), w)
    return out

if __name__ == "__main__":
    import sys, time
    t = time.time()
    print(E(4 / 2.4, 4 / 2.4, 1.55), time.time() - t)
