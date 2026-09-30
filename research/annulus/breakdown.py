import sys, numpy as np
sys.argv = ['x', '44/27', '44/27', '1.585', '4391/2700']
import indep_model_annulus as M
g = float(sys.argv[1]) if False else None
for g in [1.40, 1.50, 1.53, 1.56, 1.60]:
    S0, U0 = M.box_S0(g), M.box_U0(g)
    print(f"--- g={g}  E1={M.E1:.4f}  HIGH(eta=0)={M.E1+M.cost(S0,U0,g,0.0):.4f}")
    for eta in [0.9, 0.95, 1.0, 1.03, 1.06, 1.1]:
        gam = min(g, M.a0); G = max(0.0, gam + M.a0 + M.T - 4)
        terms = dict(d1=eta, d2=2*eta-G, d3=3*eta-gam, d4=4*eta-2*gam)
        hi = M.E1 + M.cost(S0, U0, g, eta)
        lo = M.low_at(g, eta)
        print(f"  eta={eta:.2f} HIGH={hi:.4f} LOW={lo:.4f} normals={M.normals(g,eta):.4f} " + ' '.join(f"{k}={v:.3f}" for k,v in terms.items()))
