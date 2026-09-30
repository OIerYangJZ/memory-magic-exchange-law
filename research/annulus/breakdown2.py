import sys, numpy as np
al = float(sys.argv[1]); a0 = 4/al; b = (8-2*a0)/3+0.002; T = a0-0.002
sys.argv = ['m', repr(a0), repr(a0), repr(b), repr(T)]
import model2 as M
gam_of = lambda g: min(g, a0)
for g in [float(x) for x in sys.argv[2:]] if False else [1.25, 1.41, 1.5, 1.6]:
    S0, U0 = M.box_S0(g), M.box_U0(g)
    print(f"--- alpha={al} g={g} E1={M.E1:.4f} S0={S0:.3f} U0={U0:.3f}")
    for eta in np.arange(0.8, 1.21, 0.05):
        yj, yp = eta, eta-0.01
        gam = gam_of(g); pm = M.pieces_max(g, yj, yp)
        v01 = pm; v4 = max(0, 4*yj-2*gam)+pm
        n3 = np.maximum(0, 3*yj-gam-M.Xs)+pm; pp = M.E1+np.maximum(0, S0+U0+2+M.Xs-yp); v3 = np.max(np.minimum(n3, pp))
        hi = M.E1 + M.cost(S0, U0, g, eta)
        print(f"  eta={eta:.2f} HIGH={hi:.4f} pieces={pm:.4f} v01={v01:.4f} v3={v3:.4f} v4={v4:.4f} LOW={M.low_layer(g,yj,yp):.4f}")
