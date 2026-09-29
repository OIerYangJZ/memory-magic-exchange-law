"""Full-range scan of the Result-E profile with (i) notes' formulas, (ii) corrected formulas.
corrected: sagitta = xL^2/r_Y (not min(xS, .)); level-1 radial extent max(lam^2,rho^2);
patch cap at g=0 (r_Y ~ R') is e_s = lam R' (sphere may contain the core arc), which differs from
the notes' formula when 2b' < a'."""
import numpy as np
gs = np.linspace(0,2.0,201)
def cost_correct(S0,U0,g):
    K = -2-g
    M = np.maximum(S0,U0); m = np.minimum(S0,U0)
    ssum = np.where(K/4<=m, K/2, m+(K-m)/3)
    c = S0+U0-ssum
    c = np.where(3*M+m<K, 0.0, c)
    return np.maximum(c,0.0)
def cost_claimed(S0,U0,g,h=0.01):
    # brute force per g (vectorised over a small grid of reductions)
    out=np.zeros_like(g)
    d=np.arange(0,4,h)
    for i,(s0,u0,gg) in enumerate(zip(S0,U0,g)):
        Sg,Ug=np.meshgrid(s0-d,u0-d,indexing='ij')
        l=np.maximum(Sg,Ug); sh=np.minimum(Sg,Ug)
        ok = l+sh+np.minimum(sh,2*l-(1-gg))+3<0
        out[i]=((s0-Sg)+(u0-Ug))[ok].min()
    return out
def level2(ap,bp,mode):
    g=gs
    p0 = np.where(g<=ap, np.maximum.reduce([np.full_like(g,bp),g,(g+ap)/2]), np.maximum(g,bp))
    if mode=='correct':
        p0 = np.where(g==0, bp, p0)
    q0 = np.maximum(ap,g)
    S0,U0=1-p0,1-q0
    c = cost_correct(S0,U0,g) if mode=='correct' else cost_claimed(S0,U0,g)
    return c.max()
def lvl1_ok(ap,bp,mode):
    if mode=='correct':
        return (1-bp)+2*(1-ap)+(1-2*min(ap,bp))+4 < 0
    return 2*ap+3*bp>8
def E(a0,mode,h=0.01):
    best=(np.inf,None)
    for ap in np.arange(a0,a0+2.0,h):
        for bp in np.arange(0,3.0,h):
            if not lvl1_ok(ap,bp,mode): continue
            E1=2*(ap-a0)+bp
            if E1>=best[0]: break
            e=E1+level2(ap,bp,mode)
            if e<best[0]: best=(e,(round(ap,2),round(bp,2)))
            break  # smallest feasible bp for this ap is enough? (level2 not monotone in bp) -> handled below
    return best
def E_full(a0,mode,h=0.01):
    best=(np.inf,None)
    for ap in np.arange(a0,a0+2.0,h):
        for bp in np.arange(0,3.0,h):
            if not lvl1_ok(ap,bp,mode): continue
            E1=2*(ap-a0)+bp
            if E1>=best[0]: break
            e=E1+level2(ap,bp,mode)
            if e<best[0]: best=(e,(round(ap,2),round(bp,2)))
    return best
if __name__=="__main__":
    import sys
    mode=sys.argv[1]
    taus=[float(x) for x in sys.argv[2:]]
    for tau in taus:
        e,arg=E_full(4/tau,mode)
        b=e*tau/4
        print(f"{mode} tau={tau:.3f}L E={e:.3f} {arg} bits={b:.3f}L  min(bits,tau/2)={min(b,tau/2):.3f}L tau/f={tau/min(b,tau/2):.3f}",flush=True)
