# A Memory–Magic Exchange Law in Streaming Clifford+T Compilation

Source, data and scripts for the manuscript. Every number, table and figure in the paper is
produced by a script in `scripts/` from an input in `data/`, and this file says which.

The manuscript is `main.tex` (REVTeX 4.2, PRX Quantum style). Build it with
`latexmk -pdf -outdir=out main.tex`; the PDF is written to `out/main.pdf`.

**Status.** Unpublished. All results are proven except Theorem 9 and Corollary 7, which are
conditional on Conjecture H, and the achievable rates of Corollary 3 and Proposition 3, which
are conditional on the synthesis hypotheses stated there (Eq. (31) for Proposition 3). The constant `c` of Conjecture H is deliberately left unspecified: the
enumeration bounds it from below (`c >= 163.9`, see `astra1.py`) and a finite scan cannot
bound it from above.

## Environment

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

`numpy` and `mpmath` carry everything except the figures (`matplotlib`) and the
Ross–Selinger synthesis. `gridsynth_run.py` needs `pygridsynth`; `run_mixed_achievable.py`
needs either `pygridsynth` (`--backend pygridsynth`) or the Haskell `gridsynth` binary on
`PATH` (`--backend cli`).

Run scripts from the repository root:

```bash
.venv/bin/python scripts/<name>.py
```

The one exception is `scripts/enum2.py`, which reads and writes its pickle relative to the
current directory and must be run from inside `data/`.

## What produces what

| Paper | Script | Input / output in `data/` |
|---|---|---|
| Table 1, Fig. 1 (coset counts to `t = 22`) | `enum_arith.py`, `arith_analysis.py`, `fig1_fig.py` | `enum_arith.pkl`, `enum_res2.pkl` |
| Table 2 (low-cost words per frame) | `lowcost_counts.py` | `lowcost_counts.json` |
| Table 3 (the bounds at `eps = 1e-10`) | `thm_numbers.py` | — |
| Fig. 2(a) (exact optima, both error budgets) | `frontier_eps_half.py`, `make_fig2.py` | `taumin_eps_half.json`, `frontier_L{5,6}_eps2.json` |
| Fig. 2(b) (Ross–Selinger at chemistry scale) | `gridsynth_run.py`, `make_fig2.py` | `gridsynth_tau.json`, `gridsynth_tau_epsover2.json`, `gridsynth_cache.json` |
| Sec. 9, exact `tau_min` for `L <= 6` | `taumin.py` | `taumin.pkl`, `taumin_exact.json` |
| Sec. 9, axis cosets: first scan (25 frames x 4 accuracies) | `conj1_scan.py` | `conj1_scan_output.txt` |
| Sec. 9, axis cosets: second scan (272 frames x 5 accuracies, 1360 cells) | `scan2.py`, then `scan2_report.py` | `scan2.jsonl`, `scan2_summary.txt` |
| Sec. 9, required-`c` profiles to `t = 17` | `ht_profile.py`, `profile2.py` | `ht_profile_output.txt`, `profile2_{a,b}.txt` |
| Sec. 9, the identity-coset step (232 words, `c >= 163.9`) | `astra1.py` | `astra1_output.txt` |
| Sec. 10.1, all constants (Thm. 14, Cor. 6/7, Prop. 3) | `mixing_bounds.py` | `mixing_bounds_output.txt` |
| Sec. 10.1, single-word `gridsynth` means at the mixing accuracies | `run_mixed_achievable.py` | `u_20000.txt`, `mixed_cache.json`, `mixed_achievable_results.json` |
| Lemma 13 and its `tau_0 >= 1` hypothesis | `check_lem_cost.py` | — (stdout) |
| Base enumeration (Matsumoto–Amano normal forms, resumable) | `enum2.py` | `enum_res2.pkl` |
| Thm. 4/5 (elementary tube bound, rate 11/5): exact rational certificate | `research/cert/certificate.py 9/5` | `research/cert/certificate_output.txt` |
| Rate 2.23 with optimised box shape (remark after Thm. 5) | `research/cert/certificate_223.py 179/100 179/100 737/500 400/223` | `research/cert/certificate_223_output.txt` |
| Thm. 6/7 (height dichotomy, rate 17/7; proof in App. F, research notes D0–D7 in `research/dioph/DIOPH_NOTES.md`): exact rational certificate | `research/dioph/cert_dioph.py 28/17 28/17 157/100 11183/6800 400 100 12` (needs numpy; use `.venv/bin/python`) | `research/dioph/cert_17_7_output.txt` |
| Rate 12/5, same method, larger margin (paragraph after Thm. 7; needs the relaxed level-1 bound D1 of the notes) | `research/dioph/cert_dioph.py 5/3 5/3 31/20 41/25 400 100` | `research/dioph/cert_12_5_output.txt` |

`scan2_report.py` recomputes, from the stored `scan2.jsonl`, every count quoted in the
axis-coset paragraph of Sec. 9 — the 272/1360/20400 totals, the 17 cells above `c = 48`,
the per-accuracy maxima, and the linear-term check.

## Quick reproduction

Cheap (seconds to a couple of minutes, all from cached data):

```bash
.venv/bin/python scripts/thm_numbers.py      # Table 3
.venv/bin/python scripts/mixing_bounds.py    # every constant of Sec. 10.1
.venv/bin/python scripts/scan2_report.py     # every count in the axis-coset paragraph
.venv/bin/python scripts/astra1.py           # the 232-word step forcing c >= 163.9
.venv/bin/python scripts/check_lem_cost.py   # Lemma 13 random-instance check
.venv/bin/python scripts/make_fig2.py        # Fig. 2
```

Expensive, and only needed to rebuild the cached data:

| | cost |
|---|---|
| `enum2.py` to `t = 22` | ~580 s, ~1.5 GB, single core. Each extra `t` doubles both. |
| `scan2.py` | hours; resumable, appends to `data/scan2.jsonl` |
| `profile2.py 0` / `profile2.py 1` | ~2 min each |
| `frontier_eps_half.py 6 24` | ~17 s, ~2 GB |
| `run_mixed_achievable.py` | ~2 min on 10 cores with a warm `mixed_cache.json` |

## Determinism

Every randomised choice is seeded, so the outputs above are reproducible bit for bit
(timing lines excepted):

- Haar and perturbed frames: `numpy.random.default_rng(914)` in `conj1_scan.py`; `scan2.py`
  inherits the same generator by `exec`.
- The 20 000 grid angles of Fig. 2(b) and of the mixing measurement:
  `random.Random(20260903).randrange(7853981633)`, stored as `data/u_20000.txt`.
- `gridsynth_run.py`: `--seed 20260903` by default.
- `check_lem_cost.py`: `random.seed(916)`.

Synthesis is cached by content — `(angle, accuracy) -> T-count` in `gridsynth_cache.json`
and `mixed_cache.json` — so reruns cost nothing and do not depend on the synthesiser's own
tie-breaking.

## Conventions worth knowing before reading the scripts

- Accuracy is the up-to-phase operator norm, `d_proj(A,B) = inf_phi ||A - e^{i phi} B||`,
  which is what `gridsynth` measures and what the paper calls `eps`. With `pygridsynth`
  this is `up_to_phase=True`; passing the default instead measures a different quantity.
- `L := log2(1/eps)`, and the tuned modulus is `Q_eps = floor(pi / (2 arcsin 2 eps))`.
- The frontier slope is `E_u tau / log2 K` with `u` uniform on all of `Z_Q`, the zero share
  costing nothing — not the mean over nonzero shares, which is larger by `Q/(Q-1)`.
- Fig. 2(a) reports two error budgets. The certified `r = 2` point of Theorem 8 synthesizes
  each committed share at `eps/2`; evaluating at `eps` instead is a single-rotation
  calibration and not a point of the family. The two sequences approach 3 from opposite
  sides (2.478, 2.721, 3.113 against 3.333, 3.334, 3.205).

`data/enum_res2.pkl.bak-t19` is the `t <= 19` enumeration, kept for comparison.
