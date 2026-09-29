import pickle, math, itertools, numpy as np, sys
S = math.sqrt(2)
N, R = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else "reps.pkl", "rb"))
P = np.array([[a + b * S for (a, b) in x] for x in R])
r = math.sqrt(N[0] + N[1] * S)
def mul(x, y): return (x[0]*y[0] + 2*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def sub(x, y): return (x[0]-y[0], x[1]-y[1])
def det3(u, v, w):
    t1 = mul(u[0], sub(mul(v[1], w[2]), mul(v[2], w[1])))
    t2 = mul(u[1], sub(mul(v[0], w[2]), mul(v[2], w[0])))
    t3 = mul(u[2], sub(mul(v[0], w[1]), mul(v[1], w[0])))
    return add(sub(t1, t2), t3)
def coplanar(idx):
    pts = [R[i] for i in idx]
    base = pts[0]; d = [tuple(sub(p[k], base[k]) for k in range(3)) for p in pts[1:]]
    for a, b, c in itertools.combinations(range(len(d)), 3):
        if det3(d[a], d[b], d[c]) != (0, 0): return False
    return True
n = len(P); spacing = r * math.sqrt(4 * math.pi / n)
D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
print(f"#points={n}  r={r:.4f}  mean spacing={spacing:.4f} ({spacing/r:.4f} r)")
for f in [2.0, 1.0, 0.5, 0.3, 0.2, 0.1, 0.05]:
    rho = f * spacing
    best, bestidx = 0, None
    for i in range(n):
        idx = np.where(D[i] <= 2 * rho)[0]   # points within a ball of radius rho around some center: use pairs within 2rho
        # greedy: points within rho of point i
        idx = np.where(D[i] <= rho)[0]
        if len(idx) > best: best, bestidx = len(idx), idx
    # count clusters with >=4 points within rho of a point, and how many are non-coplanar
    noncop = 0; tot = 0
    for i in range(n):
        idx = np.where(D[i] <= rho)[0]
        if len(idx) >= 4:
            tot += 1
            if not coplanar(list(idx)): noncop += 1
    print(f"rho = {f:>4} x spacing: max points in a ball = {best:3d};  balls with >=4 points: {tot:4d}, of which non-coplanar: {noncop}")
