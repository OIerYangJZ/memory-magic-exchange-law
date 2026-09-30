import sys, glob, re, math
def touchard(m, E):
    # Poisson(E) raw moment: sum_j S(m,j) E^j
    S=[[0]*(m+1) for _ in range(m+1)]; S[0][0]=1
    for n in range(1,m+1):
        for k in range(1,n+1): S[n][k]=k*S[n-1][k]+S[n-1][k-1]
    return sum(S[m][j]*E**j for j in range(m+1))
files=sorted(glob.glob('mom_tau*.txt'), key=lambda f:int(re.findall(r'\d+',f)[0]))
rows={}
for f in files:
    tau=int(re.findall(r'\d+',f)[0])
    for line in open(f):
        m=re.match(r'(WORDS|CONTROL) E=(\S+) max=(\d+) moments: (.*)',line)
        if m:
            kind,E,mx,ms=m.group(1),float(m.group(2)),int(m.group(3)),[float(x) for x in m.group(4).split()]
            rows[(tau,E,kind)]=(mx,ms)
Es=sorted({k[1] for k in rows})
for E in Es:
    print(f"\nE={E:g}: ratio WORDS/Poisson (control/Poisson in brackets) for orders 2,4,6,8,10; max words (control)")
    for tau in sorted({k[0] for k in rows}):
        if (tau,E,'WORDS') not in rows: continue
        mw,w=rows[(tau,E,'WORDS')]; mc,c=rows[(tau,E,'CONTROL')]
        # normalise by the empirical mean to remove radius discretisation: use Poisson at the empirical mean
        Ew=w[0]; Ec=c[0]
        rw=[w[o-1]/touchard(o,Ew) for o in (2,4,6,8,10)]
        rc=[c[o-1]/touchard(o,Ec) for o in (2,4,6,8,10)]
        print(f" tau={tau:2d} " + " ".join(f"{a:6.2f}[{b:4.2f}]" for a,b in zip(rw,rc)) + f"   max {mw} ({mc})")
