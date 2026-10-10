"""Re-derives the numbers of Sec. IV, Fig. 2 and Table S1 (Supplemental Sec. S1) from data/enum_res2.pkl:
the (t, L) exponent fit, the Haar-normalized ratio ranges, and the table rows.

    python3 enum_analysis.py [--tmax 22]
"""
import argparse
import os
import pickle

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(HERE, os.pardir))

ap = argparse.ArgumentParser()
ap.add_argument("--tmax", type=int, default=None)
ap.add_argument("--thresh", type=int, default=30)
ap.add_argument("--rows", type=int, nargs="*", default=[12, 15, 19, 22])
args = ap.parse_args()

res = pickle.load(open("data/enum_res2.pkl", "rb"))


def tuned_Q(eps):
    """Q_eps = floor(pi / (2 arcsin 2 eps)); 12, 25, 50, 100, 201 for L = 4..8."""
    return int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))


def c_grid(eps):
    """Exact Haar constant for the grid count.

    The grid acceptance region occupies 1/3 of the tube, and the tube carries
    36*2^t*eps^2 words; the grid spacing is 2*pi/Q_eps, not exactly 8 eps, so
    c = 36 * (1/3) * 8 eps / (2 pi / Q_eps) = 48 eps Q_eps / pi, which is 12
    only up to the rounding of Q_eps.
    """
    return 48 * eps * tuned_Q(eps) / np.pi


C_TUBE = 36.0

ts = sorted({k[0] for k in res})
if args.tmax is not None:
    ts = [t for t in ts if t <= args.tmax]
epss = sorted({k[2] for k in res}, reverse=True)
names = ["I,I", "haar0", "haar1"]
haar = ["haar0", "haar1"]
print(f"t = {min(ts)}..{max(ts)}   eps = 2^-{[int(-np.log2(e)) for e in epss]}")
print(f"words at t={max(ts)}: {res[(max(ts), 'haar0', epss[0])][2]:,}")


def fit(kind, use, thresh):
    """log2 N = a t + b L + c over cells with N >= thresh.  kind 0=grid, 1=tube."""
    A, y = [], []
    for t in ts:
        for nm in use:
            for e in epss:
                N = res[(t, nm, e)][kind]
                if N >= thresh:
                    A.append([t, -np.log2(e), 1.0]); y.append(np.log2(N))
    A = np.array(A); y = np.array(y)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef, len(y)


for kind, lab in [(0, "grid"), (1, "tube")]:
    for use, ulab in [(haar, "haar only"), (names, "all three")]:
        (a, b, c), n = fit(kind, use, args.thresh)
        print(f"{lab:5s} {ulab:10s}: log2 N = {a:.3f} t {b:+.3f} L + {c:.3f}"
              f"   ({n} cells >= {args.thresh})")

# Ratio ranges must be cut on the EXPECTED count, not the observed one: selecting
# cells with N >= c conditions on the upper Poisson tail and biases the ratio high,
# which is why an observed-count cut gives an asymmetric interval around the Haar
# value while an expected-count cut does not.
print("\nexact Haar constants:")
print("  c_grid = 48 eps Q_eps / pi = "
      + ", ".join(f"{c_grid(e):.4f}" for e in epss)
      + f"   (Q = {[tuned_Q(e) for e in epss]});  c_tube = 36 for every eps")

print("\nN / (Haar prediction), Haar pairs, by cut on the expected count:")
for kind, lab, cst in [(0, "grid", c_grid), (1, "tube", lambda e: C_TUBE)]:
    for E in (0, 100, 300, 1000, 3000):
        r = [res[(t, nm, e)][kind] / (cst(e) * 2.0 ** t * e ** 2)
             for t in ts for nm in haar for e in epss
             if cst(e) * 2.0 ** t * e ** 2 >= E]
        if len(r) < 3:
            continue
        print(f"  {lab} expected >= {E:5d}: [{min(r):.4f}, {max(r):.4f}]  "
              f"mean {np.mean(r):.4f}  n={len(r):3d}")
    r = [(res[(t, nm, e)][kind] / (cst(e) * 2.0 ** t * e ** 2), res[(t, nm, e)][kind])
         for t in ts if t >= 12 for nm in haar for e in epss]
    r2 = [x for x, N in r if N >= args.thresh]
    print(f"  {lab} [biased] observed >= {args.thresh}, t >= 12: "
          f"[{min(r2):.4f}, {max(r2):.4f}]  mean {np.mean(r2):.4f}  n={len(r2)}")

tmax = max(ts)
print(f"\nat t={tmax}, pooling the two Haar pairs:  N/(2^t eps^2) vs exact prediction")
for kind, lab, cst in [(0, "grid", c_grid), (1, "tube", lambda e: C_TUBE)]:
    obs = [sum(res[(tmax, nm, e)][kind] for nm in haar) / (2 * 2.0 ** tmax * e ** 2)
           for e in epss]
    pre = [cst(e) for e in epss]
    print(f"  {lab} obs  " + ", ".join(f"{v:.2f}" for v in obs))
    print(f"  {lab} pred " + ", ".join(f"{v:.2f}" for v in pre)
          + "   dev " + ", ".join(f"{100 * (o / p - 1):+.2f}%"
                                  for o, p in zip(obs, pre)))

# the claim "Haar pairs, t >= 12 => grid ratio in [11.4,12.6]" as literally stated
bad = [(t, nm, int(-np.log2(e)), res[(t, nm, e)][0], round(12 * 2.0 ** t * e ** 2, 1),
        round(res[(t, nm, e)][0] / (2.0 ** t * e ** 2), 2))
       for t in ts if 12 <= t <= 19 for nm in haar for e in epss]
out = [b for b in bad if not 11.4 <= b[5] <= 12.6]
print(f"\ncells with t in [12,19] on Haar pairs whose grid ratio is outside "
      f"[11.4,12.6]: {len(out)} of {len(bad)}")
for b in sorted(out, key=lambda z: abs(z[5] - 12))[-5:]:
    print(f"    t={b[0]} {b[1]} L={b[2]}: N={b[3]} expected={b[4]} ratio={b[5]}")

print("\nTable (grid counts on haar0; exact prediction c_grid*2^(t-2L) in parentheses):")
hdr = " & ".join([r"$\varepsilon=2^{-%d}$" % int(-np.log2(epss[0]))]
                 + [r"$2^{-%d}$" % int(-np.log2(e)) for e in epss[1:]])
print(f"$t$ & {hdr}\\\\")
for t in args.rows:
    if t not in ts:
        continue
    cells = []
    for e in epss:
        N = res[(t, "haar0", e)][0]
        p = c_grid(e) * 2.0 ** t * e ** 2
        cells.append(f"{N} ({p:.0f})" if p >= 10 else f"{N} ({p:.1f})")
    print(f"{t} & " + " & ".join(cells) + r"\\")
print("  deviations from the exact prediction:")
for t in args.rows:
    if t not in ts:
        continue
    d = []
    for e in epss:
        N = res[(t, "haar0", e)][0]
        p = c_grid(e) * 2.0 ** t * e ** 2
        d.append(f"{100 * (N / p - 1):+.2f}% ({(N - p) / np.sqrt(p):+.2f} sd)")
    print(f"   t={t}: " + "  ".join(d))

print("\nIdentity coset (I,I), grid counts, onset of the fee:")
for e in epss:
    L = int(-np.log2(e))
    nz = [t for t in ts if res[(t, "I,I", e)][0] > 0 and t >= 2]
    print(f"  L={L}: first t>=2 with a word = {min(nz) if nz else None}  (2L+2 = {2 * L + 2})")
