"""Sec. 9 numbers for the arithmetic cosets: reads data/enum_arith.pkl.

    python3 arith_analysis.py

Prints the exponent fit, the ratio ranges against the Haar prediction, the t=22 row,
and the onsets, split by how far the coset pair sits from the identity in the tree.
The three reference pairs of enum2.py are recomputed by enum_arith.py --check and
must agree cell by cell with data/enum_res2.pkl (they do: 345/345).
"""
import os
import pickle

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(HERE, os.pardir))

res = pickle.load(open("data/enum_arith.pkl", "rb"))
ts = sorted({k[0] for k in res})
epss = sorted({k[2] for k in res}, reverse=True)
ARITH = ["arith5", "arith10", "arith15"]
FAR = ["arith10", "arith15"]
NEAR = ["cliff", "I,I"]
HAAR = ["haar0", "haar1"]


def tuned_Q(eps):
    return int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))


def c_grid(eps):
    return 48 * eps * tuned_Q(eps) / np.pi


PRED = {0: ("grid", c_grid), 1: ("tube", lambda e: 36.0)}


def fit(use, kind, thresh=30):
    A, y = [], []
    for t in ts:
        for nm in use:
            for e in epss:
                N = res[(t, nm, e)][kind]
                if N >= thresh:
                    A.append([t, -np.log2(e), 1.0])
                    y.append(np.log2(N))
    coef, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
    return coef, len(y)


print(f"t = {min(ts)}..{max(ts)};  coset pairs: {sorted({k[1] for k in res})}")
print("\nexponent fit, cells with count >= 30:")
for kind, (lab, _) in PRED.items():
    for use, ul in [(ARITH + ["cliff"], "arithmetic"), (ARITH, "arith, no Clifford"),
                    (HAAR, "Haar (reference)")]:
        (a, b, c), n = fit(use, kind)
        print(f"  {lab:5s} {ul:20s}: log2 N = {a:.3f} t {b:+.3f} L + {c:.3f}  ({n} cells)")

print("\nN / Haar prediction, all cells with expected count >= 100:")
for kind, (lab, cst) in PRED.items():
    for use, ul in [(FAR, "T-count 10, 15"), (["arith5"], "T-count 5"),
                    (NEAR, "Clifford = identity"), (HAAR, "Haar")]:
        r = [res[(t, nm, e)][kind] / (cst(e) * 2.0 ** t * e ** 2)
             for t in ts for nm in use for e in epss if cst(e) * 2.0 ** t * e ** 2 >= 100]
        print(f"  {lab:5s} {ul:20s}: [{min(r):.3f}, {max(r):.3f}]  mean {np.mean(r):.3f}"
              f"  n={len(r)}")

tmax = max(ts)
print(f"\nat t={tmax}, ratio to the Haar prediction, per coset pair:")
for nm in sorted({k[1] for k in res}):
    g = [res[(tmax, nm, e)][0] / (c_grid(e) * 2.0 ** tmax * e ** 2) for e in epss]
    tb = [res[(tmax, nm, e)][1] / (36 * 2.0 ** tmax * e ** 2) for e in epss]
    print(f"  {nm:8s} grid " + " ".join(f"{v:5.3f}" for v in g)
          + "   tube " + " ".join(f"{v:5.3f}" for v in tb))
allr = [res[(tmax, nm, e)][k] / (PRED[k][1](e) * 2.0 ** tmax * e ** 2)
        for k in (0, 1) for nm in {kk[1] for kk in res} for e in epss]
print(f"  all {len(allr)} cells at t={tmax}: [{min(allr):.3f}, {max(allr):.3f}]")

print("\nonset (least t >= 2 carrying a grid word), against 2L:")
for nm in sorted({k[1] for k in res}):
    row = []
    for e in epss:
        L = int(-np.log2(e))
        nz = [t for t in ts if t >= 2 and res[(t, nm, e)][0] > 0]
        row.append(f"L={L}: {min(nz) if nz else '-':>3}({2 * L})")
    print(f"  {nm:8s} " + " ".join(row))

same = all(res[(t, "cliff", e)] == res[(t, "I,I", e)] for t in ts for e in epss)
print(f"\nClifford pair identical to the identity coset in every cell: {same}")
