import numpy as np, time, sys, pickle, os
H = np.array([[1,1],[1,-1]],dtype=complex)/np.sqrt(2)
S = np.array([[1,0],[0,1j]],dtype=complex)
T = np.array([[1,0],[0,np.exp(1j*np.pi/4)]],dtype=complex)
I2 = np.eye(2,dtype=complex)
def normphase(M):
    f=M.flatten(); idx=np.argmax(np.abs(f)>1e-9); return M/(f[idx]/abs(f[idx]))
def key(M): return tuple(np.round(normphase(M).flatten(),6))
cl={key(I2):I2}; fr=[I2]
while fr:
    new=[]
    for M in fr:
        for G in (H,S):
            N=G@M;k=key(N)
            if k not in cl: cl[k]=N;new.append(N)
    fr=new
CL=np.array(list(cl.values())); print('cliffords',len(CL)); assert len(CL)==24
A=H@T; B=S@H@T
def cores(k):
    if k==0: return I2[None]
    p=cores(k-1)
    return np.concatenate([np.einsum('ij,njk->nik',A,p),np.einsum('ij,njk->nik',B,p)],0)
def count_near(W,G1d,G2d,eps,Q):
    M=np.einsum('ij,njk,kl->nil',G1d,W,G2d)
    m00=M[:,0,0];m11=M[:,1,1]
    dcurve=np.sqrt(np.maximum(2-np.abs(m00)-np.abs(m11),0))
    theta=np.angle(m11)-np.angle(m00)
    d=np.round(theta/(2*np.pi/Q)).astype(int)%Q
    th=2*np.pi*d/Q
    tr=np.abs(np.exp(1j*th/2)*m00+np.exp(-1j*th/2)*m11)
    dg=np.sqrt(np.maximum(2-tr,0))
    return int((dg<=eps).sum()), int((dcurve<=eps).sum())
tmin,tmax=int(sys.argv[1]),int(sys.argv[2])
rng=np.random.default_rng(7)
def haar():
    z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));q,r=np.linalg.qr(z);return q*(np.diag(r)/np.abs(np.diag(r)))
Gs=[("I,I",I2,I2)]+[(f"haar{i}",haar(),haar()) for i in range(2)]
epss=[2.0**-k for k in (4,5,6,7,8)]
res=pickle.load(open('enum_res2.pkl','rb')) if os.path.exists('enum_res2.pkl') else {}
t0=time.time()
for t in range(tmin,tmax+1):
    Wc=I2[None] if t==0 else np.concatenate([np.einsum('ij,njk->nik',T,cores(t-1)),cores(t)],0)
    for name,G1,G2 in Gs:
        G1d=G1.conj().T;G2d=G2.conj().T
        for eps in epss:
            Q=int(np.floor(np.pi/(2*np.arcsin(2*eps))))
            tot=totc=0
            for C in CL:
                W=np.einsum('nij,jk->nik',Wc,C)
                a,b=count_near(W,G1d,G2d,eps,Q);tot+=a;totc+=b
            res[(t,name,eps)]=(tot,totc,len(Wc)*24)
    pickle.dump(res,open('enum_res2.pkl','wb'))
    print(f"t={t} words={len(Wc)*24} {time.time()-t0:.0f}s",flush=True)
