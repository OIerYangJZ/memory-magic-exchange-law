import numpy as np, sys, time
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
exec(open(_os.path.join(_HERE, 'conj1_scan.py')).read().split('# ---- 1. reproduce')[0].replace('TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14','TMAX = 12'))
eps = 0.024532; Q = Qeps(eps); print("Q_eps =", Q)
g, tb = counts(I2, I2, eps, Q)
print("identity frame, eps=0.024532, Q=32: grid shells", g.tolist(), "cum@11 =", g[:12].sum(), "; tube shells", tb.tolist())
print("allowance c=48:", 8+22+48*2**11*eps**2, " c=96:", 8+22+96*2**11*eps**2, " required c:", (g[:12].sum()-30)/(2**11*eps**2))
# which grid points d are hit at t=10, 11, and how far
thr = 2 - eps**2; th = 2*np.pi*np.arange(Q)/Q; e1, e2 = np.exp(1j*th/2), np.exp(-1j*th/2)
for t in (10, 11):
    W = SH[t]; m0, m1 = W[:,0,0], W[:,1,1]
    v = np.abs(m0[:,None]*e1[None,:] + m1[:,None]*e2[None,:]); hit = v >= thr
    dd = np.where(hit.any(axis=1))[0]
    d_of = [int(np.argmax(v[i])) for i in dd]
    dist = [np.sqrt(2 - v[i].max()) for i in dd]
    print(f"t={t}: {len(dd)} words; grid points hit (d: count) = { {d: d_of.count(d) for d in sorted(set(d_of))} }")
    print(f"       distances: min {min(dist):.5f} max {max(dist):.5f}; distinct (rounded 1e-6): {sorted(set(np.round(dist,6)))}")
