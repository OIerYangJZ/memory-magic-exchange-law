# Max 3-volume |d1^d2^d3| of 4 points in a geodesic rectangle (long xL, short xS) on a 2-sphere of radius r.
import numpy as np
rng=np.random.default_rng(0)
def pts(r,xL,xS,n):
    s=rng.uniform(-xL/2,xL/2,n)/r; u=rng.uniform(-xS/2,xS/2,n)/r
    # geodesic coords: rotate along great circle (s) then offset (u)
    return r*np.stack([np.cos(u)*np.cos(s),np.cos(u)*np.sin(s),np.sin(u)],1)
def maxvol(r,xL,xS,trials=200000):
    P=pts(r,xL,xS,4*trials).reshape(trials,4,3)
    D=P[:,1:]-P[:,:1]
    return np.abs(np.linalg.det(D)).max()
r=1.0
for xL,xS in [(1.0,0.01),(0.5,0.01),(0.3,0.001),(0.1,0.01),(0.05,0.01)]:
    v=maxvol(r,xL,xS)
    print(f"xL={xL} xS={xS}: max vol={v:.3e}  claimed xL*xS*min(xS,xL^2/r)={xL*xS*min(xS,xL**2/r):.3e}  xL*xS*xL^2/r={xL*xS*xL**2/r:.3e}")
