"""Float (non-rigorous) re-optimisation of (a,b) using EXACTLY the lemma set written in main.tex
(E1: 2a+3b>8, a>=max(a0,b), a0+a>=2b; E2 patch; E3 hyps; E4 N with 2(U0-XS)), to test the '~2.23' claim."""
import numpy as np
XL = np.arange(-3.0, 1.0, 0.01)[:, None]; XS = np.arange(-4.0, 1.0, 0.01)[None, :]
gs = np.concatenate([[0.0], np.arange(0.01, 2.0001, 0.01)])

def level2(a, b):
    worst = 0.0
    for g in gs:
        if g == 0: S0 = 1 - b
        else: S0 = min(1 - b, (2 - a - min(g, a)) / 2, 1 - g)
        U0 = min(1 - a, 1 - g)
        ok = (XS <= XL) & (XL < 1 - g) & (2 * XL - 1 <= XS) & (XL + XS + np.minimum(XS, 2 * XL - (1 - g)) + 3 < 0)
        N = np.maximum.reduce([np.zeros_like(XL + XS), S0 - XL + 0 * XS, U0 - XS + 0 * XL, S0 + U0 - XL - XS, 2 * (U0 - XS) + 0 * XL])
        N = np.where(ok, N, np.inf)
        worst = max(worst, N.min())
    return worst

def E(a0):
    best = (9, None)
    for a in np.arange(a0, a0 + 0.3, 0.005):
        for b in np.arange(max(0, (8 - 2 * a) / 3 + 1e-4), min(a, (a0 + a) / 2) + 1e-9, 0.005):
            if 2 * a + 3 * b <= 8: continue
            e1 = 2 * (a - a0) + b
            if e1 >= best[0]: continue
            val = e1 + level2(a, b)
            if val < best[0]: best = (val, (round(a, 3), round(b, 3)))
    return best

for al in [2.20, 2.22, 2.23, 2.24]:
    a0 = 4 / al
    e, ab = E(a0)
    print(f"alpha={al}: a0={a0:.4f} E={e:.4f} E<a0? {e < a0}  (a,b)={ab}")
print("fixed (1.82,1.46) level2 =", level2(1.82, 1.46))
