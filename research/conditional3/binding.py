import numpy as np
from cont_model import *
mu=0.002; E1=(8-2*(1.6+mu))/3; T=E1+10/33*mu
P=params(mu,T)
print("a0",P["a0"],"E1",E1)
for g in [0.85,1.0,1.2,1.4,1.55]:
    y=g-0.4+mu/11
    lde=max(1-g-P["a0"],1-2*P["a0"]); Dmax=1-2*g
    rows=[]
    for tang,D in [(True,None)]+([(False,lde+f*(Dmax-lde)) for f in np.linspace(0,1,11)] if Dmax>lde else []):
        a0,a,b=P["a0"],P["a"],P["b"]
        ldb=max(1-g-a,1-2*a); UT=min(1-a0,1-g); U0b=min(1-a,1-g)
        if tang: offc=lde+1; ST=min((1+lde)/2,1-g); S0b=min(1-b,(1+ldb)/2,1-g)
        else: offc=D+1; ST=min((1+lde)/2,lde+(1-D)/2,1-g); S0b=min(1-b,(1+ldb)/2,ldb+(1-D)/2,1-g)
        offs=max(0,y+offc); tube=max(0,(ST+2*UT-y)/3+1)
        boxes=max(0,ST-(1-b))+2*max(0,UT-(1-a)); per=min(tube,boxes+cost_upper(S0b,U0b,1-g,y))
        rows.append((offs+per,offs,per,tang,D))
    pm=max(rows)
    print(f"g={g:.2f} y={y:.3f}  worst piece: total={pm[0]:.3f} offsets={pm[1]:.3f} points/section={pm[2]:.3f} tangent={pm[3]} D={pm[4]}  high={high(P,g,y):.3f}")
