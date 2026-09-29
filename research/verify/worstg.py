import numpy as np, fullscan as F
ap,bp=1.82,1.46
g=np.linspace(0,2,401)
p0=np.where(g<=ap,np.maximum.reduce([np.full_like(g,bp),g,(g+ap)/2]),np.maximum(g,bp)); q0=np.maximum(ap,g)
c1=F.cost_claimed(1-p0,1-q0,g); c2=F.cost_correct(1-p0,1-q0,g)
for name,c in [("claimed",c1),("correct",c2)]:
    top=g[c>=c.max()-0.005]
    print(name,"max",round(c.max(),3),"attained for g in",round(top.min(),3),"-",round(top.max(),3))
