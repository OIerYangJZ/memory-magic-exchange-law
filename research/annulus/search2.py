import sys, importlib, numpy as np
from concurrent.futures import ProcessPoolExecutor
def run(args):
    a0, a, b, T, flags, gstep = args
    sys.argv = ['m', repr(a0), repr(a), repr(b), repr(T), flags]
    import model2; importlib.reload(model2)
    M = model2
    if not (a >= a0 and a >= b and a0 + a >= 2*b): return (a, b, None)
    if M.E1 > T: return (a, b, 9.0)
    worst = -1
    for g in np.arange(0.01, 2.02, gstep):
        v, eta = M.value(g)
        worst = max(worst, v)
        if worst > T + 0.05: break
    return (round(a,4), round(b,4), round(worst - a0, 5), round(M.kappa,4), round(M.E1,4))
if __name__ == '__main__':
    al = float(sys.argv[1]); flags = sys.argv[2]
    a0 = 4/al; T = a0 - 0.002
    grid = []
    for da in [0, 0.01, 0.02, 0.04]:
        for b in np.arange(1.30, 1.66, 0.02):
            grid.append((a0, a0 + da, float(b), T, flags, 0.04))
    with ProcessPoolExecutor(10) as ex:
        res = [r for r in ex.map(run, grid) if r[2] is not None]
    res.sort(key=lambda r: r[2])
    print(al, res[:6], flush=True)
