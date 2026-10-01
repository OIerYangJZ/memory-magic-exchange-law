"""Brute-force check of sl.c for small k: enumerates every u (not only orbit representatives).

Usage: python3 brute.py k beta prefix   (prefix = output prefix of ./sl k beta prefix)
"""
import math
import sys

import numpy as np

k, beta, pre = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
P = 2**k
nb = int(round(2 ** (beta * k)))
N = nb // 8
nphi = math.ceil(math.pi * nb / 8)
SQ2 = math.sqrt(2)
R = math.isqrt(P)

us = {}
nu = 0
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        for c in range(-R, R + 1):
            for d in range(-R, R + 1):
                A = a * a + b * b + c * c + d * d
                if A > P:
                    continue
                B = a * b + b * c + c * d - d * a
                if 2 * B * B > (P - A) ** 2:
                    continue
                nu += 1
                if A == 0:
                    continue
                us.setdefault((A, B), []).append(math.atan2(c + (b + d) / SQ2, a + (b - d) / SQ2))

S = np.zeros((nb, N + 1))
L = np.zeros(nb, dtype=np.int64)
cells = np.zeros((nb, nphi))
ells = 8 * np.arange(N + 1)
for (A, B), tu in us.items():
    tt = us.get((P - A, -B))
    if tt is None:
        continue
    j = min(int((A + B * SQ2) / P * nb), nb - 1)
    L[j] += 1
    tu = np.array(tu)
    tt = np.array(tt)
    diff = (tt[None, :] - tu[:, None]).ravel()
    S[j] += np.cos(np.outer(ells, diff)).sum(axis=1)
    ph = np.mod(diff, math.pi / 4)
    ip = np.minimum((ph * nb / 2).astype(int), nphi - 1)
    np.add.at(cells[j], ip, 1.0 / 64)

SA = np.fromfile(pre + "_A.bin", dtype=np.float64).reshape(nb, N + 1)
LA = np.fromfile(pre + "_L.bin", dtype=np.int64)
CB = np.fromfile(pre + "_B.bin", dtype=np.uint32).reshape(nb, nphi)
info = dict(line.split() for line in open(pre + ".txt"))
print("nu brute", nu, "C", info["nu"])
print("max |S_brute - S_C| =", np.abs(S - SA).max(), " (scale", np.abs(S).max(), ")")
print("L equal:", np.array_equal(L, LA))
print("max |cells_brute - cells_C| =", np.abs(cells - CB).max(), " total", cells.sum(), CB.sum())
bad = np.argwhere(np.abs(cells - CB) > 1e-9)
print("mismatched cells (z-bin, phi-bin):", sorted(set(map(tuple, bad.tolist())))[:12], "phi-bins used:", sorted(set(bad[:, 1].tolist())))
