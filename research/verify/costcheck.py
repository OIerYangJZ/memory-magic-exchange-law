import numpy as np, fullscan as F
rng=np.random.default_rng(1); worst=0
d=np.arange(0,4,0.002)
for _ in range(300):
    g=rng.uniform(0,2); S0=rng.uniform(-2,1-g); U0=rng.uniform(-3,1-g)
    a=F.cost_correct(np.array([S0]),np.array([U0]),np.array([g]))[0]
    Sg,Ug=np.meshgrid(S0-d,U0-d,indexing='ij'); l=np.maximum(Sg,Ug); sh=np.minimum(Sg,Ug)
    ok=3*l+sh<-2-g; b=((S0-Sg)+(U0-Ug))[ok].min()
    worst=max(worst,abs(a-b))
print("max |analytic - brute| =",worst)
