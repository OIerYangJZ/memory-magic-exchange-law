"""ht_profile.py -- required constant of Conjecture 1 on the HT-axis frame as a function of tau, up to t=17,
at eps = 2^-5 and 2^-6, plus a Haar frame for comparison.  Needs conj1_scan.py in the same directory
(it re-uses its enumeration and counting code by exec).  Output: ht_profile_output.txt.  Runtime ~4 min, ~1 GB."""
import numpy as np, time, sys
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
exec(open(_os.path.join(_HERE, 'conj1_scan.py')).read().split('# ---- 1. reproduce')[0].replace('TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14','TMAX = 17'))
G1, G2 = frames['axis HT (exact)']
for k in (5, 6):
    eps = 2.0**-k; Q = Qeps(eps)
    g, tb = counts(G1, G2, eps, Q)
    print(f"\nHT-axis frame, eps=2^-{k}, Q={Q}")
    print("  tau : N_grid(<=tau)  c_req_grid   | N_tube(<=tau)  c_req_tube  | shell_tube/Haar")
    for tau in range(8, TMAX+1):
        Ng, Nt = g[:tau+1].sum(), tb[:tau+1].sum()
        den = 2**tau*eps**2
        print(f"  {tau:3d} : {Ng:8d} {(Ng-8-2*tau)/den:11.1f}   | {Nt:8d} {(Nt-8-2*tau)/den:11.1f}  | {tb[tau]/(36*den):6.2f}")
for name in ('Haar 0',):
    G1, G2 = frames[name]
    for k in (5, 6):
        eps = 2.0**-k; Q = Qeps(eps); g, tb = counts(G1, G2, eps, Q)
        print(f"\n{name}, eps=2^-{k}: c_req_grid by tau 10..{TMAX}: " + " ".join(f"{(g[:tau+1].sum()-8-2*tau)/(2**tau*eps**2):.1f}" for tau in range(10, TMAX+1)))
