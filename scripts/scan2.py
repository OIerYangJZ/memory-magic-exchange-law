import numpy as np, time, sys
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
src = open(_os.path.join(_HERE, 'conj1_scan.py')).read().split('# ---- 1. reproduce')[0].replace("TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14", "TMAX = 14")
exec(src)
EPS = [2.0**-k for k in (3,4,5,6,7)]
# all distinct axes of words with t<=3, plus axes of a sample of t=4,5 words, plus Haar, identity, Clifford, perturbations
allframes = {}
allframes['identity'] = (I2, I2); allframes['Clifford pair'] = (CL[5], CL[17])
for i in range(5): allframes[f'Haar {i}'] = (haar(), haar())
for k, (t, U) in ax_list:
    D = eigframe(U); allframes[f'axis t={t} {np.round(k,3)}'] = (D, D.conj().T)
ax45 = {}
for t in (4, 5):
    for U in SH[t]:
        n = axis_of(U)
        if n is None: continue
        kk = axis_key(n)                     # canonical key, defined in conj1_scan.py
        if kk not in axes: ax45.setdefault(kk, (t, U))
ax45l = sorted(ax45.items(), key=lambda kv: kv[0]); idx = rng.choice(len(ax45l), size=40, replace=False)
for i in idx:
    kk, (t, U) = ax45l[i]; D = eigframe(U); allframes[f'axis t={t} {np.round(kk,3)}'] = (D, D.conj().T)
for ang in (0.003, 0.01, 0.03):
    Dp = rot(rng.normal(size=3), ang) @ D_HT; allframes[f'HT axis perturbed {ang}'] = (Dp, Dp.conj().T)
# mixed frames: HT axis on the left, a cheap word on the right
for t in (1, 2, 3):
    V = SH[t][rng.integers(len(SH[t]))]; allframes[f'HT axis, G2 = D^dag * (t={t} word)'] = (D_HT, D_HT.conj().T @ V)
print(f"frames: {len(allframes)}, eps: {EPS}, tau<=14", flush=True)

def analyse(g, eps):
    """returns (c_post, tau_post), (c_all, tau_all), first tau violating the linear term by >1 word when 2^tau eps^2 < 1/8"""
    cp = (0.0, None); ca = (0.0, None); lin = None
    for tau in range(TMAX + 1):
        N = g[:tau+1].sum(); den = 2**tau*eps**2
        if N > 8 + 2*tau:
            r = (N - 8 - 2*tau)/den
            if r > ca[0]: ca = (r, tau)
            if den >= 1 and r > cp[0]: cp = (r, tau)
            if den < 1/8 and lin is None: lin = (tau, int(N))
    return cp, ca, lin

import json, os
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 270
done = set()
if os.path.exists(_d('scan2.jsonl')):
    for line in open(_d('scan2.jsonl')): done.add(json.loads(line)[0])
t0 = time.time(); res = []
for name, (G1, G2) in allframes.items():
    if name in done: continue
    if time.time() - t0 > BUDGET: print("BUDGET reached"); break
    for eps in EPS:
        Q = Qeps(eps); g, tb = counts(G1, G2, eps, Q)
        cp, ca, lin = analyse(g, eps)
        _, ca_t, _ = analyse(tb, eps)
        res.append((name, eps, cp, ca, lin, ca_t, g.tolist()))
    with open(_d('scan2.jsonl'), 'a') as f:
        for r in res[-len(EPS):]: f.write(json.dumps(r) + '\n')
    print(f"  {name:40s} done {time.time()-t0:.0f}s", flush=True)
res = [tuple(json.loads(l)) for l in open(_d('scan2.jsonl'))]
print(f"frames done: {len(res)//len(EPS)} / {len(allframes)}")
print("\n=== summary ===")
print("max post-onset c (2^tau eps^2 >= 1):")
for r in sorted(res, key=lambda r: -r[2][0])[:12]:
    print(f"  {r[2][0]:6.1f} @tau={r[2][1]:<3} eps=2^{int(np.log2(r[1])):<3} {r[0]}")
print("max all-tau c (any tau with N > 8+2tau):")
for r in sorted(res, key=lambda r: -r[3][0])[:12]:
    print(f"  {r[3][0]:6.1f} @tau={r[3][1]:<3} eps=2^{int(np.log2(r[1])):<3} {r[0]}   shells {r[6][:r[3][1]+1]}")
print("cells with N > 8+2tau while 2^tau eps^2 < 1/8 (linear-term pressure):")
for r in res:
    if r[4] is not None: print(f"  eps=2^{int(np.log2(r[1]))} tau={r[4][0]} N={r[4][1]} allowance={8+2*r[4][0]}  {r[0]}")
print("max tube c (any tau):")
for r in sorted(res, key=lambda r: -r[5][0])[:5]:
    print(f"  {r[5][0]:6.1f} @tau={r[5][1]} eps=2^{int(np.log2(r[1]))} {r[0]}")
