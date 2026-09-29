import numpy as np, time, pickle, os
_HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(_HERE, 'enum2.py')).read().split("tmin,tmax=")[0])   # reuse Cliffords, cores, count fn
epss=[2.0**-k for k in (4,5,6)]
taumin={}; nwords={}
for eps in epss:
    Q=int(np.floor(np.pi/(2*np.arcsin(2*eps)))); taumin[eps]=-np.ones(Q,dtype=int); nwords[eps]={}
t0=time.time()
for t in range(0,20):
    Wc=I2[None] if t==0 else np.concatenate([np.einsum('ij,njk->nik',T,cores(t-1)),cores(t)],0)
    for eps in epss:
        Q=int(np.floor(np.pi/(2*np.arcsin(2*eps))))
        hit=np.zeros(Q,dtype=int)
        for C in CL:
            W=np.einsum('nij,jk->nik',Wc,C)
            m00=W[:,0,0];m11=W[:,1,1]
            theta=np.angle(m11)-np.angle(m00); d=np.round(theta/(2*np.pi/Q)).astype(int)%Q
            th=2*np.pi*d/Q; tr=np.abs(np.exp(1j*th/2)*m00+np.exp(-1j*th/2)*m11)
            ok=np.sqrt(np.maximum(2-tr,0))<=eps
            np.add.at(hit,d[ok],1)
        tm=taumin[eps]; new=(tm<0)&(hit>0); tm[new]=t; nwords[eps][t]=hit.copy()
    print(t,f"{time.time()-t0:.0f}s",flush=True)
pickle.dump((taumin,nwords),open(os.path.join(_HERE, os.pardir, 'data', 'taumin.pkl'),'wb'))
for eps in epss:
    Q=len(taumin[eps]); L=-np.log2(eps); tm=taumin[eps]
    print(f"eps=2^-{int(L)} Q={Q}: covered {np.sum(tm>=0)}/{Q}; taumin by d:",list(tm))
    cov=tm[tm>=0]; print("   mean taumin (covered, excl d=0) = %.2f ; 3L=%d ; max=%d"%(cov[1:].mean() if len(cov)>1 else 0,3*L,cov.max()))
