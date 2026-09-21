"""Conjecture H on *arithmetic* cosets: G1 R_z(theta) G2 with G1, G2 in Gamma.

enum2.py tests the identity coset and two Haar-random coset pairs.  The cosets that
Corollary 2 actually feeds to Theorem 2 are G1 = (U^{<t})^dagger, G2 = (U^{>t})^dagger,
both Clifford+T words, so the uniformity of Conjecture H over G1, G2 has to be tested
on Gamma itself.  This script counts the same two quantities (words within eps of a
grid rotation on the coset, and words within eps of the whole coset) for coset pairs
built from random Matsumoto-Amano words of prescribed T-count, plus one Clifford pair.

    python3 enum_arith.py 0 22            # writes ../data/enum_arith.pkl, resumable
    python3 enum_arith.py 0 16 --check    # also recompute enum2's three pairs and diff

Counting is identical to enum2.py's count_near; the eps loop is hoisted out of the two
matrix products, which is what makes four extra coset pairs affordable.  --check
verifies the hoisted version against data/enum_res2.pkl cell by cell.
"""
import argparse
import os
import pickle
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data")
# reuse the Cliffords, the syllable matrices A, B and the core enumeration
exec(open(os.path.join(HERE, "enum2.py")).read().split("tmin,tmax=")[0])

ap = argparse.ArgumentParser()
ap.add_argument("tmin", type=int)
ap.add_argument("tmax", type=int)
ap.add_argument("--check", action="store_true",
                help="also compute enum2's (I,I) and two Haar pairs and diff them")
args = ap.parse_args()

EPSS = [2.0 ** -k for k in (4, 5, 6, 7, 8)]
QS = {e: int(np.floor(np.pi / (2 * np.arcsin(2 * e)))) for e in EPSS}


def ma_word(rng, t):
    """A random Matsumoto-Amano word of minimal T-count exactly t: (T|eps)(HT|SHT)^* C."""
    if t == 0:
        return CL[rng.integers(24)], 0
    lead = bool(rng.integers(2))          # leading T present?
    nsyl = t - 1 if lead else t
    W = np.eye(2, dtype=complex)
    for _ in range(nsyl):
        W = (A if rng.integers(2) == 0 else B) @ W
    if lead:
        W = T @ W
    return W @ CL[rng.integers(24)], t


def counts(M):
    """(grid count, tube count) per eps, from M = G1^dag W G2^dag for all words at once."""
    m00, m11 = M[:, 0, 0], M[:, 1, 1]
    a00, a11 = np.abs(m00), np.abs(m11)
    dcurve = np.sqrt(np.maximum(2 - a00 - a11, 0))          # distance to the whole coset
    theta = np.angle(m11) - np.angle(m00)                   # optimal rotation angle
    out = {}
    for e in EPSS:
        Q = QS[e]
        th = (2 * np.pi / Q) * np.round(theta / (2 * np.pi / Q))
        tr = np.abs(np.exp(0.5j * th) * m00 + np.exp(-0.5j * th) * m11)
        dg = np.sqrt(np.maximum(2 - tr, 0))
        out[e] = (int((dg <= e).sum()), int((dcurve <= e).sum()))
    return out


rng = np.random.default_rng(20260909)
pairs = [("arith5", *[ma_word(rng, 5)[0] for _ in range(2)]),
         ("arith10", *[ma_word(rng, 10)[0] for _ in range(2)]),
         ("arith15", *[ma_word(rng, 15)[0] for _ in range(2)]),
         ("cliff", CL[rng.integers(24)], CL[rng.integers(24)])]

if args.check:                      # enum2.py's own pairs, same seed, same order
    r7 = np.random.default_rng(7)

    def haar():
        z = r7.normal(size=(2, 2)) + 1j * r7.normal(size=(2, 2))
        q, r = np.linalg.qr(z)
        return q * (np.diag(r) / np.abs(np.diag(r)))

    pairs += [("I,I", I2, I2)] + [(f"haar{i}", haar(), haar()) for i in range(2)]

path = os.path.join(DATA, "enum_arith.pkl")
res = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
t0 = time.time()
for t in range(args.tmin, args.tmax + 1):
    Wc = I2[None] if t == 0 else np.concatenate(
        [np.einsum('ij,njk->nik', T, cores(t - 1)), cores(t)], 0)
    tot = {(nm, e): [0, 0] for nm, _, _ in pairs for e in EPSS}
    for name, G1, G2 in pairs:
        G1d, G2d = G1.conj().T, G2.conj().T
        for C in CL:
            M = np.einsum('ij,njk,kl->nil', G1d, np.einsum('nij,jk->nik', Wc, C), G2d)
            for e, (g, c) in counts(M).items():
                tot[(name, e)][0] += g
                tot[(name, e)][1] += c
    for (name, e), (g, c) in tot.items():
        res[(t, name, e)] = (g, c, len(Wc) * 24)
    pickle.dump(res, open(path, "wb"))
    print(f"t={t} words={len(Wc) * 24} {time.time() - t0:.0f}s", flush=True)

if args.check:
    ref = pickle.load(open(os.path.join(DATA, "enum_res2.pkl"), "rb"))
    bad = [(k, res[k], ref[k]) for k in res
           if k in ref and res[k] != ref[k]]
    n = sum(1 for k in res if k in ref)
    print(f"regression against enum_res2.pkl: {n - len(bad)}/{n} cells identical")
    for b in bad[:10]:
        print("   MISMATCH", b)
