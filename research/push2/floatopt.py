import numpy as np
def cost(a,b,g,zero=False):
    S0 = (1-b) if zero else min(1-b,(2-a-min(g,a))/2,1-g); U0=min(1-a,1-g); R0=1-g
    best=np.inf
    for XS in np.arange(min(U0,R0-1e-7),-4,-0.002):
        c1=min((XS+R0)/2,(R0-XS-3)/3-1e-7); cands=[c1]
        if -3-2*XS-1e-7>(XS+R0)/2: cands.append(-3-2*XS-1e-7)
        XL=max(min(max(cands),R0-1e-7,(1+XS)/2,S0),XS)
        if not (XS<=XL<R0 and 2*XL-1<=XS and XL+XS+min(XS,2*XL-R0)+3<0): continue
        c=max(0,S0-XL,U0-XS,S0+U0-XL-XS,2*(U0-XS)); best=min(best,c)
        if c==0: break
    return best
GS=np.linspace(0,2.01,202)
def E(a0,a,b):
    return 2*(a-a0)+b+max([cost(a,b,0,True)]+[cost(a,b,g) for g in GS[1:]])
def opt(a0):
    best=(9,None)
    for a in np.arange(a0,a0+0.08,0.005):
        bmin=max((8-2*a)/3+1e-4,0)
        for b in np.arange(bmin,min(a,(a0+a)/2)+1e-9,0.005):
            e=E(a0,a,b)
            if e<best[0]: best=(e,(round(a,4),round(b,4)))
    return best
for a0 in (1.80,1.795,1.79,1.785):
    e,arg=opt(a0); print(f"a0={a0} E={e:.4f} {arg}  4/E... alpha_needed<{4/a0:.4f}  ok_if_alpha<= {4/e:.4f} and <{4/a0:.4f}",flush=True)
