"""Statistics for the outputs of sl.c.  Usage: python3 analyze.py prefix [prefix ...]

Test A: Z = S_l(W_j) / (8 sqrt(S_0(W_j))).  Under random phases E Z^2 ~ 1 and max|Z| ~ sqrt(2 ln #tests).
        |S_l| <= L^{1-theta} uniformly gives rate 1 + 1/(1-theta) (NOTES §10(A)); theta = 1/3 is the 5/2 line.
Test B: counts of orbit pairs in (z, phi) cells of side 2/nb, |z| <= 0.8, compared with Poisson.
"""
import math
import sys

import numpy as np


def pois_max(lam, ncells):
    """Smallest x with ncells * P(Pois(lam) >= x) <= 1."""
    x = int(lam)
    while True:
        # log P(X >= x) via the sum of the tail terms
        logs = [-lam + j * math.log(lam) - math.lgamma(j + 1) for j in range(x, x + 400)]
        m = max(logs)
        tail = math.exp(m) * sum(math.exp(v - m) for v in logs)
        if ncells * tail <= 1:
            return x
        x += 1


def windows(c, w):
    """Max over sliding w x w windows (no wrap)."""
    cs = np.cumsum(np.cumsum(np.pad(c, ((1, 0), (1, 0))), axis=0), axis=1)
    s = cs[w:, w:] - cs[:-w, w:] - cs[w:, :-w] + cs[:-w, :-w]
    return s.max(), s.size


for pre in sys.argv[1:]:
    info = {a: b for a, b in (line.split() for line in open(pre + ".txt"))}
    k, nb, N, nphi = int(info["k"]), int(info["nb"]), int(info["N"]), int(info["nphi"])
    npairs = int(info["npairs"])
    S = np.fromfile(pre + "_A.bin", dtype=np.float64).reshape(nb, N + 1)
    L = np.fromfile(pre + "_L.bin", dtype=np.int64)
    C = np.fromfile(pre + "_B.bin", dtype=np.uint32).reshape(nb, nphi).astype(np.int64)
    z = 2 * (np.arange(nb) + 0.5) / nb - 1
    bulk = np.abs(z) <= 0.8
    # Test A
    Z = S[bulk, 1:] / (8 * np.sqrt(S[bulk, :1]))
    ntest = Z.size
    j, n = np.unravel_index(np.argmax(np.abs(Z)), Z.shape)
    zj = z[bulk][j]
    Lm = L[bulk].mean()
    Smax = np.abs(S[bulk, 1:]).max()
    q = np.array_split(np.arange(N), 4)
    ez2q = [float((Z[:, qq] ** 2).mean()) for qq in q]
    print(f"[{pre}] k={k} nb={nb} ell<= {8 * N}  bands={bulk.sum()}  mean L={Lm:.0f}  mean S0/L={S[bulk, 0].mean() / Lm:.1f}")
    print(f"  A: E Z^2={float((Z**2).mean()):.3f} (by ell-quartile {', '.join(f'{v:.2f}' for v in ez2q)})"
          f"  E Z^4/(E Z^2)^2={float((Z**4).mean() / (Z**2).mean() ** 2):.2f}"
          f"  max|Z|={np.abs(Z).max():.2f} at z={zj:+.3f}, ell={8 * (n + 1)}  [gauss max ~ {math.sqrt(2 * math.log(ntest)):.2f}]")
    print(f"     max|S_l|={Smax:.3g}  log max|S| / log L = {math.log(Smax) / math.log(Lm):.3f}"
          f"  max|S|/sqrt(64 S0)~{np.abs(Z).max():.2f}")
    # Test B (drop the partial last phi column)
    Cb = C[bulk, : nphi - 1]
    lam = npairs * (2 / nb) ** 2 / (2 * math.pi / 4)
    mean, var = Cb.mean(), Cb.var()
    out = [f"  B: lambda={lam:.2f} (empirical {mean:.2f})  var/mean={var / mean:.3f}"]
    for w in (1, 2, 4):
        mx, nw = windows(Cb, w)
        out.append(f"w={w}: max {mx} vs Poisson-max {pois_max(lam * w * w, nw)} (mean {lam * w * w:.1f})")
    print("  ".join(out))
    sys.stdout.flush()
