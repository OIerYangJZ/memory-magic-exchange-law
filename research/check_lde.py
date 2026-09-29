"""Relation between minimal T-count and least denominator exponent k (unitary = M/sqrt2^k)."""
import collections, sys
from zomega import all_words, canon
tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
ws = all_words(tmax)
keys = set(canon(W) for _, W in ws)
print("words", len(ws), "distinct", len(keys), "expected", 72 * 2 ** tmax - 48)
tab = collections.Counter((t, W[0]) for t, W in ws)
for t in range(tmax + 1):
    print(t, sorted((k, c) for (tt, k), c in tab.items() if tt == t))
