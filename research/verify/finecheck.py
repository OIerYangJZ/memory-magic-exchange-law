import numpy as np, fullscan as F
F.gs=np.linspace(0,2.0,20001)
for tau in [2.20,2.22,2.23,2.24]:
    a0=4/tau; best=(9,None)
    for ap in np.arange(a0,a0+0.6,0.005):
        for bp in np.arange(1.2,1.8,0.005):
            if not F.lvl1_ok(ap,bp,'correct'): continue
            E1=2*(ap-a0)+bp
            if E1>=best[0]: break
            e=E1+F.level2(ap,bp,'correct')
            if e<best[0]: best=(e,(round(ap,3),round(bp,3)))
    e,(ap,bp)=best
    g=F.gs
    p0=np.where(g<=ap,np.maximum.reduce([np.full_like(g,bp),g,(g+ap)/2]),np.maximum(g,bp)); p0=np.where(g==0,bp,p0)
    c=F.cost_correct(1-p0,1-np.maximum(ap,g),g)
    print(f"tau={tau}: E={e:.4f} a'={ap} b'={bp} 2b'>=a' {2*bp>=ap} a'>=b' {ap>=bp} argmax g={g[c.argmax()]:.3f} E2={c.max():.4f} bits={e*tau/4:.4f}L")
