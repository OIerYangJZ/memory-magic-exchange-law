"""Independent re-implementation of the Result-E exponent bookkeeping.
Units: log R'. Sizes (in R' units) are exponents s (<=1). Box: rho R' = R'^{1-a'}, lam R' = R'^{1-b'}.
Level1: 2a'+3b'>8 (valid when a'>=b', else radial extent is rho^2: 4a'+b'>8).
Sphere Y radius R'^{1-g}. Patch caps: S0 = 1-p0, U0 = 1-q0.
Level2 condition on rectangle with long l >= short sh (sizes):
  mode 'claimed': l+sh+min(sh, 2l-(1-g)) + 3 < 0
  mode 'correct': l+sh+(2l-(1-g)) + 3 < 0
cost = (S0-s)+(U0-u) >= 0.
"""
import numpy as np, sys
def lvl2_cost(ap,bp,g,mode,h=0.0025):
    p0 = max(bp,g,(g+ap)/2) if g<=ap else max(bp,g)
    q0 = max(ap,g)
    S0,U0 = 1-p0,1-q0
    s = S0-np.arange(0,4,h); u = U0-np.arange(0,4,h)
    Sg,Ug = np.meshgrid(s,u,indexing='ij')
    l = np.maximum(Sg,Ug); sh=np.minimum(Sg,Ug)
    if mode=='claimed': vol = l+sh+np.minimum(sh,2*l-(1-g))
    else: vol = l+sh+(2*l-(1-g))
    ok = vol+3 < 0
    cost = (S0-Sg)+(U0-Ug)
    return cost[ok].min()
def E(a0,mode,h=0.01,gs=np.linspace(0,2.0,81)):
    best=(np.inf,None)
    for ap in np.arange(a0,a0+1.6,h):
        for bp in np.arange(max((8-2*ap)/3,0)+1e-6, 3.0, h):
            E1=2*(ap-a0)+bp
            if E1>=best[0]: break
            if not (ap>=bp and 2*bp>=ap):  # geometric side conditions used in the notes
                continue
            E2=0.0
            for g in gs:
                E2=max(E2,lvl2_cost(ap,bp,g,mode,h=0.01))
                if E1+E2>=best[0]: break
            if E1+E2<best[0]: best=(E1+E2,(round(ap,3),round(bp,3),round(E2,3)))
    return best
mode=sys.argv[1]
for tau in [1.5,1.8,2.0,2.1,2.2,2.25,2.3]:
    e,arg=E(4/tau,mode)
    print(mode,f"tau={tau}L E={e:.3f} {arg} bits={e*tau/4:.3f}L tau/bits={tau/(e*tau/4):.3f}",flush=True)
