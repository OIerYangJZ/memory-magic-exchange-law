"""Cell-count moments from sl.c Test-B outputs.  Usage: python3 moments.py prefix [prefix ...]
Factorial moments F_m = E[N(N-1)...(N-m+1)] equal the m-point correlation (ordered m-tuples of distinct
points in a common cell).  Poisson: F_m = lam^m.  Centred Poisson moments: mu3 = lam, mu4 = 3lam^2+lam,
mu6 = 15lam^3+25lam^2+lam.  Bulk |z| <= 0.8, last (partial) phi column dropped, as in analyze.py.
"""
import math, sys
import numpy as np
for pre in sys.argv[1:]:
    info = {a: b for a, b in (l.split() for l in open(pre + ".txt"))}
    k, nb, N, nphi = int(info["k"]), int(info["nb"]), int(info["N"]), int(info["nphi"])
    npairs = int(info["npairs"])
    C = np.fromfile(pre + "_B.bin", dtype=np.uint32).reshape(nb, nphi).astype(np.float64)
    z = 2 * (np.arange(nb) + 0.5) / nb - 1
    bulk = np.abs(z) <= 0.8
    Cb = C[bulk, : nphi - 1]
    lam = npairs * (2 / nb) ** 2 / (2 * math.pi / 4)
    m = Cb.mean()
    fm = []
    acc = np.ones_like(Cb)
    for j in range(8):
        acc = acc * (Cb - j)
        if j >= 1:
            fm.append(acc.mean() / m ** (j + 1))
    d = Cb - m
    c3 = (d ** 3).mean() / m
    c4 = (d ** 4).mean() / (3 * m ** 2 + m)
    c6 = (d ** 6).mean() / (15 * m ** 3 + 25 * m ** 2 + m)
    # same with the top 0.01% of cells removed (robustness against rich circles)
    q = np.quantile(Cb, 0.9999)
    Ct = Cb[Cb <= q]
    mt = Ct.mean()
    acc = np.ones_like(Ct); fmt = []
    for j in range(8):
        acc = acc * (Ct - j)
        if j >= 1:
            fmt.append(acc.mean() / mt ** (j + 1))
    print(f"[{pre}] k={k} nb={nb} cells={Cb.size} lam={lam:.2f} mean={m:.2f} var/mean={Cb.var() / m:.3f} max={int(Cb.max())}")
    print("   F_m/mean^m (m=2..8):        " + " ".join(f"{v:.3f}" for v in fm))
    print("   same, top 0.01% removed:    " + " ".join(f"{v:.3f}" for v in fmt))
    print(f"   centred/Poisson: mu3 {c3:.3f}  mu4 {c4:.3f}  mu6 {c6:.3f}")
    sys.stdout.flush()
