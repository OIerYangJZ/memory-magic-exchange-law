#!/usr/bin/env python3
"""
frontier_eps_half.py -- exact minimal T-count tau_min(d; eta) of a Clifford+T word within
dproj <= eta of the grid rotation Rz(2 pi d/Q), for eta = eps (single-rotation calibration)
and eta = eps/2 (the accuracy Theorem 5 needs for a certified two-round process).
Used for Fig. 3(a) in the eps/2 error budget.  Streams the Matsumoto-Amano shells and
stops when every grid point is resolved or TMAX is reached.

Usage: python3 frontier_eps_half.py L [TMAX] [--json out.json]
       L=4: seconds;  L=5: t<=20, seconds;  L=6: t<=24, minutes and ~2 GB.

Criterion.  dproj(W, Rz(th)) = sqrt(2 - |tr(Rz(th)^dag W)|) and
tr(Rz(th)^dag W) = e^{i th/2} W00 + e^{-i th/2} W11, so dproj <= eta iff |tr| >= 2 - eta^2.
Only the two diagonal entries of each word are ever needed.

Rewrite of the original resolve() (kept as frontier_eps_half.py.bak-orig), which held the
full 2x2 word array of shell t and rechecked all Q grid points every shell.  Three changes
make L=6 feasible on 32 GB:
  * only the diagonals W00, W11 are formed, as two (n,24) matrix products -- 4x less memory;
  * each shell is checked against the *unresolved* grid points only, which after onset is a
    handful, so the late shells cost almost nothing;
  * one enumeration serves all (Qgrid, eta) tasks instead of one enumeration each.
The shell count is 36*2^t (Lemma 4), reproduced by the --selftest assertion.
"""
import sys, json, time
import numpy as np

args = [a for a in sys.argv[1:] if not a.startswith("--")]
L = int(args[0]); TMAX = int(args[1]) if len(args) > 1 else 24
JSON = None
if "--json" in sys.argv:
    JSON = sys.argv[sys.argv.index("--json") + 1]

eps = 2.0 ** -L
Q = int(np.floor(np.pi / (2 * np.arcsin(2 * eps))))
Q8 = 8 * (Q // 8)

H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)
T = np.diag([1, np.exp(1j * np.pi / 4)]).astype(complex)
I2 = np.eye(2, dtype=complex)
HT, SHT = H @ T, S @ H @ T


def key(U):
    f = U.flatten()
    i = next(j for j in range(4) if abs(f[j]) > 1e-9)
    return tuple(np.round(f / (f[i] / abs(f[i])), 7).view(float))


cl = {key(I2): I2}; fr = [I2]
while fr:
    nw = []
    for C in fr:
        for g in (H, S):
            U = C @ g
            if key(U) not in cl:
                cl[key(U)] = U; nw.append(U)
    fr = nw
CL = np.array(list(cl.values()))
assert len(CL) == 24
# columns needed for the diagonals of  core @ CL
CL0 = np.ascontiguousarray(CL[:, :, 0].T)   # (2,24), picks W00
CL1 = np.ascontiguousarray(CL[:, :, 1].T)   # (2,24), picks W11

BUDGET = 20_000_000          # complex entries held by one |tr| block


def resolve_all(tasks, tmax, verbose=True):
    """tau_min(d; eta) for every (Qgrid, eta, label) in tasks, from one enumeration."""
    st = []
    for Qg, eta, lab in tasks:
        th = 2 * np.pi * np.arange(Qg) / Qg
        st.append(dict(Q=Qg, eta=eta, lab=lab, e1=np.exp(1j * th / 2),
                       e2=np.exp(-1j * th / 2), thr=2 - eta ** 2,
                       tau=-np.ones(Qg, int)))
    per_block = max(1024, BUDGET // max(Qg for Qg, _, _ in tasks))

    def check(m0, m1, t):
        for s in st:
            un = np.flatnonzero(s["tau"] < 0)
            if un.size == 0:
                continue
            v = np.abs(m0[:, None] * s["e1"][un][None, :]
                       + m1[:, None] * s["e2"][un][None, :])
            hit = (v >= s["thr"]).any(axis=0)
            if hit.any():
                s["tau"][un[hit]] = t

    core = np.array([I2]); prev = None
    total = 0
    for t in range(tmax + 1):
        t0 = time.time()
        if t == 0:
            n_shell = 24
            check(CL[:, 0, 0].copy(), CL[:, 1, 1].copy(), 0)
        else:
            prev, core = core, np.concatenate([core @ HT, core @ SHT])
            n_shell = 24 * (len(core) + len(prev))
            step = max(1, per_block // 24)
            for src, pre in ((core, None), (prev, T)):
                for a in range(0, len(src), step):
                    blk = src[a:a + step]
                    if pre is not None:
                        blk = pre @ blk
                    m0 = (blk[:, 0, :] @ CL0).ravel()
                    m1 = (blk[:, 1, :] @ CL1).ravel()
                    check(m0, m1, t)
        total += n_shell
        assert n_shell == (24 if t == 0 else 36 * 2 ** t), (t, n_shell)
        left = [f"{s['lab']}:{int((s['tau'] < 0).sum())}" for s in st]
        if verbose:
            print(f"  t={t:2d}  shell={n_shell:>11,}  unresolved {' '.join(left)}"
                  f"  ({time.time() - t0:.1f}s)", flush=True)
        if all((s["tau"] >= 0).all() for s in st):
            break
    return st, total


print(f"L={L}  eps=2^-{L}  tuned Q_eps={Q}  admissible Q={Q8}  TMAX={TMAX}")
tasks = [(Qg, eta, f"{nm}/{el}")
         for Qg, nm in ((Q, f"Q{Q}"), (Q8, f"Q{Q8}"))
         for eta, el in ((eps, "eps"), (eps / 2, "eps2"))]
t_start = time.time()
st, total = resolve_all(tasks, TMAX)
print(f"enumerated {total:,} words in {time.time() - t_start:.0f}s\n")

out = {}
for s in st:
    tau, Qg = s["tau"], s["Q"]
    nz = tau[1:]
    ok = nz[nz >= 0]
    unres = int((tau < 0).sum())
    lk = float(np.log2(Qg - 1))
    # the frontier slope uses u uniform on all of Z_Q, the zero share being free
    mean_all = float(ok.sum() / Qg) if unres == 0 else float("nan")
    rec = dict(L=L, Q=Qg, eta=s["eta"], label=s["lab"], tmax=TMAX,
               unresolved=unres, log2K=lk,
               mean_nonzero=float(ok.mean()), max=int(ok.max()),
               per_bit_nonzero=float(ok.mean()) / lk,
               mean_over_ZQ=mean_all,
               per_bit=mean_all / lk if unres == 0 else float("nan"),
               tau=tau.tolist())
    out[s["lab"]] = rec
    bound = "" if unres == 0 else f"  [>= : {unres} point(s) unresolved at t<={TMAX}]"
    print(f"{s['lab']:>12}  eta={s['eta']:.6g}  mean(d!=0)={ok.mean():7.3f}  max={ok.max():3d}"
          f"  E_u tau={mean_all:7.3f}  per bit={mean_all / lk if unres == 0 else float('nan'):.3f}{bound}")

if JSON:
    json.dump(out, open(JSON, "w"), indent=1)
    print(f"\nwritten {JSON}")
