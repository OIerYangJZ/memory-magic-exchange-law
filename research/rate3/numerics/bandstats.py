"""Fluctuations of the band count B_W (= 64*S_0) across bulk bands at beta = 4/5. Pure Python (no numpy)."""
import sys, math
from array import array
S = sys.argv[1]
def stats(xs):
    n = len(xs); m = sum(xs)/n
    v = sum((x-m)**2 for x in xs)/(n-1)
    return m, math.sqrt(v)
for k in map(int, sys.argv[2:]):
    info = dict(l.split() for l in open(f"{S}/k{k}.txt"))
    nb, N = int(info["nb"]), int(info["N"])
    A = array('d'); A.frombytes(open(f"{S}/k{k}_A.bin","rb").read())
    L = array('q'); L.frombytes(open(f"{S}/k{k}_L.bin","rb").read())
    S0 = [A[j*(N+1)] for j in range(nb)]
    z = [2.0*(j+0.5)/nb - 1.0 for j in range(nb)]
    bulk = [j for j in range(nb) if abs(z[j]) <= 0.8]
    s0 = [S0[j] for j in bulk]; lb = [float(L[j]) for j in bulk]
    m, s = stats(s0); mL, _ = stats(lb)
    # local fluctuation: neighbour differences remove smooth z-trends
    d = [s0[i+1]-s0[i] for i in range(len(s0)-1)]
    _, sd = stats(d); loc = sd/math.sqrt(2)
    # per-band Poisson scale: sqrt(mean count); divisor-weight variance proxy from the S0/L ratio spread
    ratio = [s0[i]/lb[i] for i in range(len(s0))]
    mr, sr = stats(ratio)
    mx = max(abs(s0[i]-m) for i in range(len(s0)))/loc
    print(f"k={k} nb={nb} bulk={len(bulk)} meanS0={m:.1f} meanL={mL:.1f} S0/L={mr:.4f}±{sr:.4f}"
          f" relstd_raw={s/m:.4f} relstd_local={loc/m:.4f} relstd_local*sqrt(S0)={loc/m*math.sqrt(m):.3f}"
          f" relstd_local*sqrt(L)={loc/m*math.sqrt(mL):.3f} max|dev|/local={mx:.2f}")
