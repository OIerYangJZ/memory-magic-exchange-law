"""Fig. 1: Conjecture H test from the exhaustive enumeration (data/enum_res2.pkl).

(a) raw counts against the Haar prediction 12*2^(t-2L);
(b) the same counts normalized by that prediction, which is what the conjecture
    asserts converges to 1.
"""
import argparse
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(HERE, os.pardir))

ap = argparse.ArgumentParser()
ap.add_argument("--name", default="haar0")
ap.add_argument("--out", default="fig1_tube_counts.png")
args = ap.parse_args()

res = pickle.load(open("data/enum_res2.pkl", "rb"))
tmax = max(k[0] for k in res)
ts = np.arange(tmax + 1)
ks = (4, 5, 6, 7, 8)


def c_grid(eps):
    """48 eps Q_eps / pi -- the Haar constant at the true grid spacing 2 pi / Q_eps,
    equal to 12 only up to the rounding of Q_eps (see enum_analysis.py)."""
    Q = int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))
    return 48 * eps * Q / np.pi

fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9))

ax = axes[0]
for i, k in enumerate(ks):
    e = 2.0 ** -k
    N = np.array([res[(t, args.name, e)][0] for t in ts], float)
    ax.semilogy(ts, np.where(N > 0, N, np.nan), "o-", ms=3, color=f"C{i}",
                label=rf"$\varepsilon=2^{{-{k}}}$")
    ax.semilogy(ts, c_grid(e) * 2.0 ** (ts - 2 * k), "--", color=f"C{i}",
                lw=0.7, alpha=0.6)
ax.set_ylim(0.5, None)
ax.set_xlabel(r"$T$-count $t$")
ax.set_ylabel(r"words within $\varepsilon$ of a grid rotation")
ax.set_xticks(np.arange(0, tmax + 1, 2))
ax.legend(fontsize=7)
ax.set_title(r"(a) counts vs. Haar prediction $c_{\rm grid}\,2^{t-2L}$", fontsize=9)

ax = axes[1]
for i, k in enumerate(ks):
    e = 2.0 ** -k
    N = np.array([res[(t, args.name, e)][0] for t in ts], float)
    pred = c_grid(e) * 2.0 ** (ts - 2 * k)
    r = np.where(pred >= 30, N / pred, np.nan)   # only where the count is meaningful
    ax.plot(ts, r, "o-", ms=3, color=f"C{i}", label=rf"$\varepsilon=2^{{-{k}}}$")
ax.axhline(1.0, color="k", lw=0.8)
ax.set_ylim(0.6, 1.4)
ax.set_xlim(axes[0].get_xlim())
ax.set_xticks(np.arange(0, tmax + 1, 2))
ax.set_xlabel(r"$T$-count $t$")
ax.set_ylabel(r"$N/(c_{\rm grid}\,2^{t-2L})$")
ax.set_title(r"(b) normalized, cells with prediction $\geq30$", fontsize=9)
ax.legend(fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig(args.out, dpi=150)
print(f"wrote {args.out}  (t <= {tmax})")
