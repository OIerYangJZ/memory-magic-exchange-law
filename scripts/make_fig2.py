"""make_fig2.py -- Fig. 2: the fractional-passthrough tradeoff family against the lower bounds.
(a) exact optimal synthesis from data/taumin_eps_half.json (L = 5, 6, m = 64, r = 2),
    at both error budgets: eps/2 (the certified r = 2 point of Thm. 9) and eps (a calibration);
(b) Ross-Selinger (gridsynth) calibration at eps = 1e-10, m = 1e4, r = 2, against Theorems 1, 3, 10
    (Theorem 10 assuming Conjecture 1 with (c0, c1, c) = (8, 2, 256), see thm_numbers.py).

Called from the repository root, like the other scripts:  .venv/bin/python scripts/make_fig2.py
Supersedes frontier_fig2.py (kept as .bak-rev1), which drew the pre-v6 panel (b)."""
import json, os, sys, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thm_numbers as tn

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))

# ------------------------------------------------ (a)
# Two error budgets.  The certified r = 2 point of Thm. 9 synthesizes each committed share
# to eps/2; evaluating tau_min at eps instead is a single-rotation calibration, not a point
# of the family.  data/taumin_eps_half.json carries both (scripts/frontier_eps_half.py);
# its taumin_eps column reproduces data/taumin_exact.json grid point by grid point.
R = json.load(open('data/taumin_eps_half.json'))
m = 64
colors = {'5': 'tab:blue', '6': 'tab:green'}
Bmax = 0
for Ls in ('5', '6'):
    d = R[Ls]; Q = d['Q']; K = d['K']; lK = np.log2(K)
    q = np.arange(0, m + 1)
    S = q * np.log2(Q); B = m * lK - S; ok = B >= 0
    mk = {'5': 'o', '6': 's'}[Ls]
    for col, lab, sty, al, lw in (('taumin_eps_half', r'\varepsilon/2', '-', 1.0, 1.5),
                                  ('taumin_eps', r'\varepsilon', '--', 0.55, 1.1)):
        tm = np.array(d[col]); Etau = tm.mean()
        ax1.plot(B[ok], ((m - q) * Etau)[ok], sty, lw=lw, color=colors[Ls], alpha=al,
                 marker=mk, ms=3.2, markevery=8, mfc='none',
                 label=r'$L=%s$ at $%s$: $\mathbb{E}_u\tau=%.2f$, slope %.3f'
                       % (Ls, lab, Etau, Etau / lK))
    # partial commitment of the low l bits, committed piece synthesized at eps
    tm = np.array(d['taumin_eps']); u = np.arange(Q)
    lmax = int(np.floor(np.log2(Q)))
    xs = [l * m for l in range(1, lmax + 1)]
    ys = [m * tm[u % 2 ** l].mean() for l in range(1, lmax + 1)]
    ax1.plot(xs, ys, 'x', color=colors[Ls], ms=7, mew=1.5,
             label=(r'partial commitment, low $\ell$ bits at $\varepsilon$' if Ls == '5' else None))
    Bmax = max(Bmax, m * lK)
ax1.plot([2 * m], [0], marker='*', mfc='none', mec='k', ms=11, ls='none', label=r'$S^k$ commitment (2 bits, free)')
ax1.plot([3 * m], [0.5 * m], 'k*', ms=11, label=r'$T^j$ commitment (3 bits, half a gate)')
Bl = np.linspace(0, Bmax, 50)
for a, ls in ((1, ':'), (2, '--'), (3, '-')):
    ax1.plot(Bl, a * Bl, ls, color='tab:red', lw=1.2, label=r'$\alpha=%d$ (leading order)' % a)
ax1.set_xlabel(r'bits shed $B=m\log_2K-S$'); ax1.set_ylabel(r'committed magic $T_{\mathrm{pre}}$')
ax1.set_title(r'(a) exact optimal synthesis, $m=64$, $r=2$, two error budgets', fontsize=10)
ax1.set_ylim(0, 1.22 * Bmax * 3)
ax1.legend(fontsize=6.0, loc='upper left', ncol=2, framealpha=0.92,
           columnspacing=1.0, handlelength=2.4)

# ------------------------------------------------ (b)
m = 1e4; k = tn.k; Q = tn.Q
bs = np.linspace(0, k, 400)
B = m * bs
ax2.plot(B, m * np.array([tn.thm1(b) for b in bs]), ':', color='tab:red', label='Thm. 1 (uncond., slope 1)')
ax2.plot(B, m * np.array([tn.thm3(b) for b in bs]), '--', color='tab:red', label='Thm. 3 (given (R), slope 2)')
ax2.plot(B, m * np.array([tn.thm9(b) for b in bs]), '-', color='tab:red',
         label=r'Thm. 10 (Conj. H with $c=%d$, slope $\kappa=%.2f$)' % (tn.c, tn.kappa))
ax2.plot(B, 3 * B, '-', color='0.6', lw=1, label='slope 3 (leading order)')
Etau = 105.34   # gridsynth mean at eps/2, the certified share accuracy for r = 2 (102.32 at eps)
q = np.linspace(0, m, 200); S = q * np.log2(Q); Bp = m * k - S; ok = Bp >= 0
ax2.plot(Bp[ok], (m - q[ok]) * Etau, '-', color='tab:blue', lw=2,
         label=r'$\mathbb{E}\tau=105.34$ at $\varepsilon/2$, slope %.3f' % (Etau / np.log2(Q)))
per_bit = [61, 46, 35, 28, 22, 19, 16, 14]
ax2.plot([l * m for l in range(1, 9)], [l * c * m for l, c in enumerate(per_bit, 1)], 'x', color='tab:blue', ms=7, mew=1.5,
         ls='none', label=r'partial commitment, low $\ell\leq8$ bits')
ax2.plot([2 * m], [0], marker='*', mfc='none', mec='k', ms=11, ls='none', label=r'$S^k$ commitment (2 bits, free)')
ax2.plot([3 * m], [0.5 * m], 'k*', ms=11, label=r'$T^j$ commitment (3 bits, half a gate)')
ax2.set_xlabel(r'bits shed $B=m\log_2K-S$'); ax2.set_ylabel(r'committed magic $T_{\mathrm{pre}}$')
ax2.set_title(r'(b) Ross--Selinger, $\varepsilon=10^{-10}$, $m=10^4$, $r=2$', fontsize=10)
ax2.set_xlim(0, 1.02 * m * k); ax2.set_ylim(0, 1.3e6); ax2.legend(fontsize=7.5, loc='upper left', bbox_to_anchor=(0.27, 0.995))
ax2.ticklabel_format(style='sci', scilimits=(0, 0), axis='both')

plt.tight_layout()
plt.savefig('fig2_frontier_twopanel.png', dpi=200)
print('wrote fig2_frontier_twopanel.png')
