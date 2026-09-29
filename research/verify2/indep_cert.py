"""Independent re-check of research/cert/certificate.py against the hypotheses AS STATED in main.tex.
Re-derives cells with the certificate's search, then re-verifies each with separately written checks:
  Lemma E3 hyps: XS <= XL, XL < 1-g, 2XL-1 <= XS, XL+XS+min(XS,2XL-(1-g))+3 < 0  (checked at g_right)
  N (Lemma E4) evaluated at g_left with S0,U0 of the Assembly paragraph.
Also: continuity/margin beyond g=2, the extra constraint 1-2a<=XS, and the heuristic ceilings 5/2, 8/3."""
from fractions import Fraction as F
import sys, importlib.util, io, contextlib

A0 = F(9, 5); a, b = F(91, 50), F(73, 50)
assert 2*a+3*b > 8 and a >= A0 and a >= b and A0 + a >= 2*b

def S0U0(g, zero=False):
    S0 = (1-b) if zero else min(1-b, (2-a-min(g, a))/2, 1-g)
    return S0, min(1-a, 1-g)

def E3ok(g, XL, XS):
    return XS <= XL and XL < 1-g and 2*XL-1 <= XS and XL+XS+min(XS, 2*XL-(1-g))+3 < 0

def Ncost(g, XL, XS, zero=False):
    S0, U0 = S0U0(g, zero)
    return max(F(0), S0-XL, U0-XS, S0+U0-XL-XS, 2*(U0-XS))

# load certificate's best_cell without running its main loop
src = open('/Users/yangjinsey/Desktop/memory-magic-exchange-law/research/cert/certificate.py').read()
head = src.split('TARGET =')[0]
ns = {}
sys.argv = ['x', '9/5']
exec(compile(head, 'cert_head', 'exec'), ns)
best_cell = ns['best_cell']

G = [F(k, 400) for k in range(0, 805)]   # extend to 2.01
worst = F(0); worst_at = None; fails = []; extra_binding = 0
c0 = best_cell(F(0), F(0), zero=True)
assert E3ok(F(0), c0[1], c0[2]); w0 = Ncost(F(0), c0[1], c0[2], zero=True)
print('g=0 cell', c0[1:], 'N', w0)
worst = w0
for gl, gr in zip(G, G[1:]):
    c = best_cell(gl, gr)
    if c is None: fails.append((gl, gr)); continue
    _, XL, XS = c
    if not E3ok(gr, XL, XS): fails.append(('E3', gl, gr))
    if not (1-2*a <= XS): extra_binding += 1
    n = Ncost(gl, XL, XS)
    if n != c[0]: fails.append(('Nmismatch', gl))
    if n > worst: worst, worst_at = n, (gl, gr, XL, XS)
print('fails:', fails[:5], len(fails))
print('worst N =', worst, float(worst), 'at', [float(x) for x in worst_at])
E = 2*(a-A0)+b+worst
print('E(9/5) =', E, float(E), ' <= 1.80501 ?', E <= F(180501, 100000), ' kappa=E/4=', float(E/4),
      ' 0.45125 >= kappa?', F(45125, 100000) >= E/4, ' 0.4513 >= kappa?', F(4513, 10000) >= E/4)

# heuristic ceilings (analytic scan): level-2 cost 0
import numpy as np
def Emin(a0, kind):
    best=9
    for A in np.linspace(a0, a0+2, 4001):
        B = max(0,(8-2*A)/3) if kind==1 else max(0,3-A)
        if kind==1 and a0+A < 2*B: continue
        best=min(best, 2*(A-a0)+B)
    return best
for kind in (1,2):
    als=[al for al in np.arange(2,3,0.0005) if Emin(4/al,kind)/4 < 1/al]
    print('ceiling kind',kind, max(als))
