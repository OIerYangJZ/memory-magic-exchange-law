import sys, numpy as np
al = float(sys.argv[1]); a0 = 4/al; b = (8-2*a0)/3 + float(sys.argv[2]); T = a0 - float(sys.argv[3])
sys.argv = ['m', repr(a0), repr(a0), repr(b), repr(T)]
import model2 as M
print(f"alpha={al} a0={a0:.5f} b={b:.5f} E1={M.E1:.5f} T={T:.5f}")
rows = []
for g in np.arange(0.01, 2.02, 0.02):
    v, eta = M.value(g)
    rows.append((v, g, eta))
rows.sort(reverse=True)
for v, g, eta in rows[:8]:
    S0, U0 = M.box_S0(g), M.box_U0(g)
    print(f"  g={g:.2f} value={v:.5f} (T-v={T-v:+.5f}) eta*={eta:.2f}  HIGH(0)={M.E1+M.cost(S0,U0,g,0.0):.4f} HIGH(eta*)={M.E1+M.cost(S0,U0,g,eta):.4f}")
