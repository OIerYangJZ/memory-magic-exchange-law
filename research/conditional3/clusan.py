import numpy as np, sys
for f in sys.argv[1:]:
    d=np.loadtxt(f,comments='#')
    tc=d[:,0].astype(int); Q=d[:,3:7]
    print(f, 'points',len(d))
    for t in sorted(set(tc)):
        P=Q[tc==t]
        if len(P)<5:
            print(f'  T-count {t}: {len(P)} pts'); continue
        C=P-P.mean(0); u,s,vt=np.linalg.svd(C); nrm=vt[-1]
        print(f'  T-count {t}: {len(P)} pts, sv ratios {np.array2string(s/s[0],precision=2)}, normal {np.round(nrm/np.abs(nrm).max(),6)}')
