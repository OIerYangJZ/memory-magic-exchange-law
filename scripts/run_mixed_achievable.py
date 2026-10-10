#!/usr/bin/env python3
"""
run_mixed_achievable.py -- measure the achievable committed cost per share
under two-word mixing (Proposition S2, Supplemental Sec. S4); the result is the
"two-word gridsynth mixing" entry of Table S3.

What it does
  * takes the SAME angle set as the paper's Fig. 3(b) / Table I point
    (20 000 grid angles u in Z_Q, Q = Q_eps = 7 853 981 633 at eps = 1e-10;
    angle = 2*pi*u/Q), or any list of angles you give it;
  * runs gridsynth at the word accuracy of two-word mixing,
        d_proj <= sqrt(eps/r)  with r = 2  ->  7.0710678e-06
    (and, optionally, sqrt(eps/4) = 5e-06 as a spare point);
  * reports mean / s.d. / range of the T-count, the per-bit ratio
    tau / log2(K_eps), and the comparison numbers of Table S3.

Cost model (Prop. 3): the two-word Campbell/Hastings program uses two words
of the same accuracy with opposite error vectors; its expected T-count is
the mean of the two, i.e. the gridsynth T-count at that accuracy plus O(1)
for finding the partner word. We report the gridsynth number and say so.

Usage
  python3 scripts/run_mixed_achievable.py --u-file data/u_20000.txt   # grid indices, one per line
  python3 run_mixed_achievable.py --angles angles.txt             # radians, one per line
  python3 run_mixed_achievable.py --regen 20000 --seed 20260912   # fresh uniform u (only if the
                                                                  # original file is lost; then say so in the paper)
Options
  --eps 7.0710678e-06 5e-06     accuracies (d_proj, = gridsynth's -e)
  --backend cli|pygridsynth     default cli (the Haskell `gridsynth` binary on PATH)
  --jobs N                      parallel workers (default: all cores)
  --cache cache.json            content-addressed cache (angle string, eps) -> T-count
Runtime: gridsynth at 7e-6 is fast (well under 0.1 s per angle); 20 000 angles
with 10 workers is a few minutes on the M5.
"""
import argparse, json, math, os, subprocess, sys, time
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
from concurrent.futures import ProcessPoolExecutor, as_completed

EPS_PAPER = 1e-10
Q_EPS = 7853981633            # tuned modulus at eps = 1e-10 (paper, Sec. IV)
LOG2K = 32.871                # log2 K_eps (paper)

def gridsynth_cli(angle_str, eps):
    """T-count of the gridsynth word for angle (decimal string, radians) at accuracy eps."""
    out = subprocess.run(["gridsynth", "-e", repr(eps), angle_str],
                         capture_output=True, text=True, check=True).stdout
    # gridsynth prints the gate string (Matsumoto-Amano, right-to-left) on the last non-empty line
    line = [l for l in out.strip().splitlines() if l.strip()][-1]
    return line.count("T")

def gridsynth_py(angle_str, eps):
    """pygridsynth. up_to_phase=True is the d_proj convention of the paper
    (Paper 1, Rem. 6), matching scripts/gridsynth_run.py."""
    import mpmath
    from pygridsynth import GridsynthConfig, gridsynth_gates
    mpmath.mp.dps = 40
    gates = gridsynth_gates(mpmath.mpf(angle_str), mpmath.mpf(repr(eps)),
                            cfg=GridsynthConfig(up_to_phase=True))
    return str(gates).count("T")

def work(args):
    angle_str, eps, backend = args
    f = gridsynth_cli if backend == "cli" else gridsynth_py
    return angle_str, eps, f(angle_str, eps)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--u-file"); ap.add_argument("--angles")
    ap.add_argument("--regen", type=int); ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--eps", type=float, nargs="+", default=[math.sqrt(EPS_PAPER / 2), math.sqrt(EPS_PAPER / 4)])
    ap.add_argument("--backend", default="cli", choices=["cli", "pygridsynth"])
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--cache", default=_d("mixed_cache.json"))
    ap.add_argument("--out", default=_d("mixed_achievable_results.json"))
    a = ap.parse_args()

    # ---- angles as exact decimal strings (2*pi*u/Q to 40 digits, so the cache key is stable)
    import mpmath
    mpmath.mp.dps = 50
    if a.u_file:
        us = [int(x) for x in open(a.u_file) if x.strip()]
        angles = [mpmath.nstr(2 * mpmath.pi * u / Q_EPS, 40) for u in us]
    elif a.angles:
        angles = [mpmath.nstr(mpmath.mpf(x.strip()), 40) for x in open(a.angles) if x.strip()]
    elif a.regen:
        import random
        rng = random.Random(a.seed)
        us = [rng.randrange(Q_EPS) for _ in range(a.regen)]
        angles = [mpmath.nstr(2 * mpmath.pi * u / Q_EPS, 40) for u in us]
        print(f"[note] regenerated {a.regen} uniform u with seed {a.seed}; these are NOT the paper's angles")
    else:
        sys.exit("give --u-file, --angles or --regen")

    cache = json.load(open(a.cache)) if os.path.exists(a.cache) else {}
    todo = [(s, e, a.backend) for e in a.eps for s in angles if f"{s}|{e!r}" not in cache]
    print(f"{len(angles)} angles x {len(a.eps)} accuracies; {len(todo)} to compute, {len(cache)} cached")
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(work, t) for t in todo]
        for fu in as_completed(futs):
            s, e, tc = fu.result(); cache[f"{s}|{e!r}"] = tc; done += 1
            if done % 500 == 0:
                json.dump(cache, open(a.cache, "w"))
                print(f"  {done}/{len(todo)}  {time.time()-t0:.0f}s", flush=True)
    json.dump(cache, open(a.cache, "w"))

    results = {}
    for e in a.eps:
        tcs = [cache[f"{s}|{e!r}"] for s in angles]
        n = len(tcs); mean = sum(tcs) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in tcs) / (n - 1))
        L_word = math.log2(1 / e)
        results[repr(e)] = dict(n=n, mean=mean, sd=sd, sem=sd / math.sqrt(n), min=min(tcs), max=max(tcs),
                                per_bit=mean / LOG2K, word_bits=L_word, rs_pred=3 * L_word + 2.7)
        delta = 2 * e * e   # diamond budget of the two-word mixture: 2 d_proj^2
        print(f"\neps_word = {e:.4e}  (log2(1/eps_word) = {L_word:.2f}; two-word mixture diamond error <= {delta:.2e})")
        print(f"  gridsynth T-count: mean {mean:.2f} +- {sd/math.sqrt(n):.2f} (s.d. {sd:.2f}, range {min(tcs)}-{max(tcs)}, n={n})")
        print(f"  per bit of the share: {mean/LOG2K:.3f}   [Ross-Selinger scaling estimate 3 log2(1/eta) + 2.7: {3*L_word+2.7:.2f}]")
        print(f"  reference: mixed-diagonal formula of Kliuchnikov et al. (Quantum 7, 1208) at delta={delta:.1e}: {1.52*math.log2(1/delta)-0.01:.2f}")
        print(f"  -> Table S3 entry: 'two-word gridsynth mixing: {mean:.1f} T per share at eps=1e-10, r=2' (plus O(1) for the partner word)")
    json.dump(results, open(a.out, "w"), indent=2)
    print(f"\nwritten {os.path.relpath(a.out)}")

if __name__ == "__main__":
    main()
