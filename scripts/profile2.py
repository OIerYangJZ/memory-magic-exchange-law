"""profile2.py -- same as ht_profile.py for four further axes (SHT and three axes flagged by scan2) at eps=2^-6, 2^-7,
up to t=17.  Usage: python3 profile2.py 0  and  python3 profile2.py 1  (two halves; outputs profile2_a.txt, profile2_b.txt).
Needs conj1_scan.py in the same directory."""
import numpy as np, time, sys, json
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
exec(open(_os.path.join(_HERE, 'conj1_scan.py')).read().split('# ---- 1. reproduce')[0].replace('TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14','TMAX = 17'))
targets = {}
want = {'A':(0.0,0.383,0.924), 'B':(-0.548,-0.548,0.632), 'C':(0.471,0.471,0.746)}
for t in range(1,6):
    for U in SH[t]:
        n = axis_of(U)
        if n is None: continue
        if n[np.argmax(np.abs(n))] < 0: n = -n
        for lab, w in want.items():
            if lab not in targets and np.allclose(n, w, atol=2e-3): targets[lab] = (t, U)
    if len(targets)==3: break
fr = {'SHT axis': frames['axis SHT (exact)']}
for lab,(t,U) in targets.items():
    D = eigframe(U); fr[f'axis {want[lab]} (t={t} word)'] = (D, D.conj().T)
which = sys.argv[1]
t0=time.time()
for name,(G1,G2) in list(fr.items())[int(which)::2]:
    for k in (6,7):
        eps=2.0**-k; Q=Qeps(eps); g,tb=counts(G1,G2,eps,Q)
        prof=[(tau, int(g[:tau+1].sum()), round((g[:tau+1].sum()-8-2*tau)/(2**tau*eps**2),1), round(tb[tau]/(36*2**tau*eps**2),2)) for tau in range(10,TMAX+1)]
        print(f"{name} eps=2^-{k} Q={Q}: (tau, N, c_req, shell_tube/Haar) = {prof}   [{time.time()-t0:.0f}s]", flush=True)
