# Ancillary files: A Memory–Magic Exchange Law in Streaming Clifford+T Compilation

These files are the code, data and exact certificates behind every number, table and figure of
the paper, frozen with this arXiv version. The same material is maintained at
https://github.com/OIerYangJZ/memory-magic-exchange-law.

## Environment

Python 3 with the packages in `requirements.txt`. The data were produced with Python 3.14,
numpy 2.5.2, mpmath 1.4.1, matplotlib 3.11.1 and pygridsynth 2.0.0. Run every command from
this directory (the one containing `scripts/` and `data/`):

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/thm_numbers.py
```

## The exact certificates (parts of the proofs)

Pure Python with exact rational arithmetic (`fractions.Fraction`); `cert_dioph.py` also uses
numpy, but only to propose witnesses that are then re-checked exactly. Each runs in seconds.

| Paper | Command | Stored output |
|---|---|---|
| Theorems 4–5 (elementary tube bound, rate 11/5), App. E | `python research/cert/certificate.py 9/5` | `research/cert/certificate_output.txt` |
| Rate 2.23 (paragraph after Theorem 5) | `python research/cert/certificate_223.py 179/100 179/100 737/500 400/223` | `research/cert/certificate_223_output.txt` |
| Theorems 6–7 (height dichotomy, rate 17/7), App. F | `python research/dioph/cert_dioph.py 28/17 28/17 157/100 11183/6800 400 100 12` | `research/dioph/cert_17_7_output.txt` |
| Rate 12/5 (paragraph after Theorem 7; conditional on a relaxed Lemma 14) | `python research/dioph/cert_dioph.py 5/3 5/3 31/20 41/25 400 100` | `research/dioph/cert_12_5_output.txt` |
| Float re-implementation quoted in App. F (largest exponent 1.6347) | `python research/verify4/indep_model.py` | stdout |

`research/dioph/DIOPH_NOTES.md` contains the working notes (D0–D7) from which Appendix F was
written; `research/verify4/` also has the exact checks of the covolume and duality facts used
in Lemmas 19 and 23.

## What produces what

| Paper | Script (in `scripts/`) | Input / output in `data/` |
|---|---|---|
| Table I, Fig. 1 (coset counts to T-count 22) | `enum_arith.py`, `arith_analysis.py`, `fig1_fig.py` | `enum_arith.pkl`, `enum_res2.pkl` |
| Table II (low-cost words per frame) | `lowcost_counts.py` | `lowcost_counts.json` |
| Table III (the bounds at eps = 1e-10, c = 48 and 256) | `thm_numbers.py` | — |
| Fig. 2(a) (exact optima, both error budgets) | `frontier_eps_half.py`, `make_fig2.py` | `taumin_eps_half.json`, `frontier_L{5,6}_eps2.json` |
| Fig. 2(b) (Ross–Selinger at eps = 1e-10) | `gridsynth_run.py`, `make_fig2.py` | `gridsynth_tau.json`, `gridsynth_tau_epsover2.json`, `gridsynth_cache.json` |
| Sec. IX, exact tau_min for L <= 6 | `taumin.py` | `taumin.pkl`, `taumin_exact.json` |
| Sec. IX, axis cosets, first scan (25 frames x 4 accuracies) | `conj1_scan.py` | `conj1_scan_output.txt` |
| Sec. IX, axis cosets, second scan (272 frames x 5 accuracies) | `scan2.py`, then `scan2_report.py` | `scan2.jsonl`, `scan2_summary.txt` |
| Sec. IX, required-c profiles to t = 17 | `ht_profile.py`, `profile2.py` | `ht_profile_output.txt`, `profile2_{a,b}.txt` |
| Sec. IX, identity-coset step (232 words, c >= 163.9) | `astra1.py` | `astra1_output.txt` |
| Sec. X A, all constants (Theorem 14, Corollaries 6–7, Proposition 3) | `mixing_bounds.py` | `mixing_bounds_output.txt` |
| Sec. X A, single-word gridsynth means at the mixing accuracies | `run_mixed_achievable.py` | `u_20000.txt`, `mixed_cache.json`, `mixed_achievable_results.json` |
| Lemma 13 and its hypothesis tau_0 >= 1 | `check_lem_cost.py` | stdout |
| Base enumeration (Matsumoto–Amano normal forms) | `enum2.py` (run from inside `data/`) | `enum_res2.pkl` |

Cheap (seconds to minutes, from the stored data): `thm_numbers.py`, `mixing_bounds.py`,
`scan2_report.py`, `astra1.py`, `check_lem_cost.py`, `make_fig2.py`, `fig1_fig.py`.
Expensive, only needed to rebuild the stored data: `enum2.py` to t = 22 (~10 min, ~1.5 GB),
`scan2.py` (hours, resumable), `profile2.py`, `ht_profile.py`, `frontier_eps_half.py`,
`gridsynth_run.py`, `run_mixed_achievable.py`.

Every random choice is seeded (`numpy.random.default_rng(914)` for frames,
`random.Random(20260903)` for the 20 000 grid angles in `data/u_20000.txt`), and synthesis
results are cached by (angle, accuracy) in `data/gridsynth_cache.json` and `data/mixed_cache.json`,
so reruns reproduce the stored outputs (timing lines excepted).

Conventions: accuracy is the up-to-phase operator norm d_proj (gridsynth's `up_to_phase=True`),
L = log2(1/eps), and the tuned modulus is Q_eps = floor(pi / (2 arcsin 2 eps)).
