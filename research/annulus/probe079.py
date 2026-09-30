import sys, os, numpy as np
al=float(sys.argv[1]); a0=4/al; b=(8-2*a0)/3+1e-4; T=a0-1e-4
sys.argv=['m',repr(a0),repr(a0),repr(b),repr(T)]
os.environ['DETA']='0.002'
import model2f as M
for g in [0.74,0.76,0.78,0.785,0.79,0.80,0.82,0.86]:
    v,eta=M.value(g)
    S0,U0=M.box_S0(g),M.box_U0(g)
    print(f"g={g} value={v:.5f} T={T:.5f} eta*={eta} HIGH0={M.E1+M.cost(S0,U0,g,0.0):.5f} S0={S0:.4f} U0={U0:.4f}")
    if eta:
        for y in np.arange(0.0, eta+1e-9, max(eta/6,0.002)):
            yj=y+0.002; yp=y
            gam=min(g,a0); pm=M.pieces_max(g,yj,yp)
            print(f"    y={yj:.3f} LOW={M.low_layer(g,yj,yp):.5f} pm={pm:.5f} HIGH={M.E1+M.cost(S0,U0,g,yj):.5f}")
