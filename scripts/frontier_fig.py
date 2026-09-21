# Regenerates fig1/fig2 from data/enum_res2.pkl and data/taumin.pkl (see chat transcript for the inline version).
import pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
taumin,_=pickle.load(open('data/taumin.pkl','rb')); res=pickle.load(open('data/enum_res2.pkl','rb'))
m=64
fig,axes=plt.subplots(1,2,figsize=(11,4.2))
for ax,eps in zip(axes,[2.0**-5,2.0**-6]):
    tm=taumin[eps]; Q=len(tm); K=Q-1; L=-np.log2(eps); lK=np.log2(K); lQ=np.log2(Q); Etau=tm.mean()
    qs=np.arange(0,m+1); S=qs*lQ; Tpre=(m-qs)*Etau; B=m*lK-S
    ax.plot(B,Tpre,'k.-',ms=3,label=f'fractional passthrough (E tau={Etau:.2f})')
    for a,ls,lab in [(1,':','alpha=1'),(2,'--','alpha=2'),(3,'-','alpha=3')]: ax.plot(B,a*B,ls,color='C3',lw=1,label=lab)
    for l in range(1,int(np.floor(lQ))+1):
        Tk=np.mean([tm[k] for k in range(2**l)]); ax.plot([m*l],[m*Tk],'x',color='C0')
    ax.set_xlabel('bits shed'); ax.set_ylabel('committed magic T_pre'); ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig2_frontier_smalleps.png',dpi=150)
fig,ax=plt.subplots(figsize=(5.5,4))
for k in (4,5,6,7,8):
    e=2.0**-k; ts=range(20); N=[res[(t,'haar0',e)][0] for t in ts]
    ax.semilogy(list(ts),[n if n>0 else np.nan for n in N],'o-',ms=3,label=f'eps=2^-{k}'); ax.semilogy(list(ts),[12*2**(t-2*k) for t in ts],'--',color='gray',lw=0.8)
ax.legend(fontsize=7); plt.tight_layout(); plt.savefig('fig1_tube_counts.png',dpi=150)
