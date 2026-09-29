"""Check covol(V^perp cap O*) vs covol(V cap O*) (= 8 H_V) for random K-planes V."""
import numpy as np, math, random
exec(open('d0_check.py').read().split("# ---- D0")[0])
def to_OK(v):   # value of <x,y> as R^2 pair -> (p,q) ints
    v1, v2 = v
    p = (v1 + v2) / 2; q = (v1 - v2) / (2 * s2)
    assert abs(p - round(p)) < 1e-6 and abs(q - round(q)) < 1e-6, (p, q)
    return round(p), round(q)
def pair(x, y): return (float(np.dot(x[:4], y[:4])), float(np.dot(x[4:], y[4:])))
def kernel(basis, funcs):
    """Z-basis of {x in Z-span(basis): <x,f> = 0 for f in funcs} (values in O_K)"""
    n = len(basis)
    rows = []
    for i, x in enumerate(basis):
        r = []
        for f in funcs: r += list(to_OK(pair(x, f)))
        rows.append(r + [1 if j == i else 0 for j in range(n)])
    m = 2 * len(funcs)
    def reduce(rows, col, start):
        piv = start
        while True:
            nz = [r for r in range(piv, len(rows)) if rows[r][col] != 0]
            if not nz: return piv
            r0 = min(nz, key=lambda r: abs(rows[r][col]))
            rows[piv], rows[r0] = rows[r0], rows[piv]
            done = True
            for r in range(piv + 1, len(rows)):
                if rows[r][col]:
                    qq = rows[r][col] // rows[piv][col]
                    rows[r] = [a - qq * b for a, b in zip(rows[r], rows[piv])]
                    if rows[r][col]: done = False
            if done: return piv + 1
    p = 0
    for col in range(m): p = reduce(rows, col, p)
    ker = np.array([r[m:] for r in rows[p:]])
    return ker @ np.array(basis)
def covol(rows): return math.sqrt(abs(np.linalg.det(rows @ rows.T)))
random.seed(3)
for trial in range(10):
    v1 = sum(random.randint(-5, 5) * Ostar[i] for i in range(8))
    v2 = sum(random.randint(-5, 5) * Ostar[i] for i in range(8))
    Vp = kernel(list(Ostar), [v1, v2])          # V^perp cap O*, rank 4
    # pick two K-independent vectors of V^perp cap O*
    u = Vp
    V = kernel(list(Ostar), [u[0], u[1], u[2], u[3]])   # (V^perp)^perp cap O* = V cap O*
    cVp, cV = covol(Vp), covol(V)
    print("rank", len(Vp), len(V), " covol(V cap O*) = 8 H_V -> H_V = %.3f   covol(V^perp cap O*) = %.3f  ratio %.5f" % (cV / 8, cVp, cVp / (cV / 8)))
