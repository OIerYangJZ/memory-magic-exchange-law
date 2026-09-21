#!/usr/bin/env python3
"""scan2_report.py -- reads scan2.jsonl (written by scan2.py) and prints every number quoted in
M2_revision_0914.md Sec. 1.1b, with the totals stated explicitly.  Run: python3 scan2_report.py"""
import json, collections, numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_DATA = _os.path.join(_HERE, _os.pardir, 'data')
def _d(name): return _os.path.join(_DATA, name)
res = [json.loads(l) for l in open(_d('scan2.jsonl'))]
# record = [frame name, eps, (c_post, tau), (c_all, tau), lin, (c_tube, tau), shell counts g[0..TMAX]]
EPS = sorted({r[1] for r in res}, reverse=True)
TMAX = len(res[0][6]) - 1
frames = sorted({r[0] for r in res})
print(f"TOTALS: frames = {len(frames)}, accuracies = {len(EPS)} (eps = {[f'2^{int(round(np.log2(e)))}' for e in EPS]}), "
      f"cells (frame, eps) = {len(res)}, tau range 0..{TMAX} -> (frame, eps, tau) triples = {len(res)*(TMAX+1)}")
assert len(res) == len(frames) * len(EPS)

def cat(name):
    if name.startswith('axis'): return name.split(' [')[0]
    if name.startswith('Haar'): return 'Haar'
    return name
byf = collections.defaultdict(list)
for r in res: byf[r[0]].append(r)
cats = collections.defaultdict(list)
for name, rs in byf.items():
    cats[cat(name)].append((max(r[2][0] for r in rs), max(r[3][0] for r in rs), max(r[5][0] for r in rs), name))
print("\nPer-category maxima over the frame's five accuracies (c_req = (N_grid(<=tau) - 8 - 2 tau)/(2^tau eps^2)):")
print(f"{'category':36s} {'#frames':>7} {'post-onset c: max':>18} {'median':>7} {'min':>6} {'any-tau max':>12} {'tube max':>9}")
for c, v in cats.items():
    vals = np.array([x[0] for x in v])
    print(f"{c:36s} {len(v):7d} {vals.max():18.1f} {np.median(vals):7.1f} {vals.min():6.1f} {max(x[1] for x in v):12.1f} {max(x[2] for x in v):9.1f}")

allmax = max(res, key=lambda r: r[2][0])
print(f"\nGlobal maximum of post-onset c: {allmax[2][0]:.1f} at tau={allmax[2][1]}, eps=2^{int(round(np.log2(allmax[1])))}, frame '{allmax[0]}'")
print(f"Global maximum of any-tau c:    {max(r[3][0] for r in res):.1f}")
over48 = [r for r in res if r[2][0] > 48]
print(f"Cells with post-onset c > 48: {len(over48)}; their shell counts are all equal: "
      f"{len({tuple(r[6][:r[2][1]+1]) for r in over48}) == 1}; shells = {over48[0][6][:over48[0][2][1]+1]}")
over40 = [r for r in res if r[2][0] > 40]
print(f"Cells with post-onset c > 40: {len(over40)}:")
for r in sorted(over40, key=lambda r: -r[2][0]): print(f"   {r[2][0]:5.1f} eps=2^{int(round(np.log2(r[1])))} tau={r[2][1]:2d}  {r[0]}")
lin = [r for r in res if r[4] is not None]
print(f"\nLinear-term check: cells with N > 8 + 2 tau while 2^tau eps^2 < 1/8: {len(lin)} of {len(res)}")
print("\nPer-accuracy maximum of post-onset c:")
for e in EPS:
    m = max((r for r in res if r[1] == e), key=lambda r: r[2][0])
    print(f"   eps=2^{int(round(np.log2(e))):<3}: {m[2][0]:5.1f} at tau={m[2][1]}  ({m[0]})")
print(f"\nAxis frames with post-onset c > 48 / > 40 / > 30: "
      f"{sum(1 for c,v in cats.items() if c.startswith('axis') for x in v if x[0]>48)} / "
      f"{sum(1 for c,v in cats.items() if c.startswith('axis') for x in v if x[0]>40)} / "
      f"{sum(1 for c,v in cats.items() if c.startswith('axis') for x in v if x[0]>30)} of "
      f"{sum(len(v) for c,v in cats.items() if c.startswith('axis'))}")
