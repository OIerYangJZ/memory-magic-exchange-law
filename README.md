# A Memory–Magic Exchange Law in Streaming Clifford+T Compilation

Source, data and scripts for the manuscript. Every number, table and figure in the paper is
produced by a script in `scripts/` from an input in `data/`, and this file says which.

The manuscript is `main.tex` (REVTeX 4.2, PRX Quantum style): main text, appendices with all
proofs of the main-text theorems, and references. The Supplemental Material is `supplement.tex`
(extended numerics, side information, model dependence and batching, the law under mixing).
The two files cross-reference each other through `xr-hyper` (labels of the other file carry
the prefix `S:` or `M:`), so build them in the same directory in the order
`pdflatex main && pdflatex supplement && pdflatex main && pdflatex supplement && pdflatex main`.

**arXiv.** [arXiv:2609.37368](https://arxiv.org/abs/2609.37368). Version 1 (tag `arxiv-v1`) proves the
unconditional rate 17/7; version 2, the present sources, proves every rate below 5/2.

**Status.** Not yet peer reviewed. All results are proven except Theorem 8 and Corollary S3, which are
conditional on Conjecture 1 (H), and the achievable rates of Corollary 4 and Proposition S2, which
are conditional on the synthesis hypotheses stated there (Eq. (S11) for Proposition S2). Numbers that
need a value of `c` assume `(c0, c1, c) = (8, 2, 256)`. The constant `c` of Conjecture H is deliberately left unspecified: the
enumeration bounds it from below (`c >= 163.9`, see `astra1.py`) and a finite scan cannot
bound it from above.

## Environment

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

`numpy` and `mpmath` carry everything except the figures (`matplotlib`), the
Ross–Selinger synthesis and the SMT check. `gridsynth_run.py` needs `pygridsynth`;
`run_mixed_achievable.py` needs either `pygridsynth` (`--backend pygridsynth`) or the Haskell
`gridsynth` binary on `PATH` (`--backend cli`, the default); `research/uniform52/verify_z3.py`
needs `z3-solver`. All three are in `requirements.txt`. The C++ example of Remark S2,
`research/conditional3/ballpts.cpp`, needs a C++17 compiler; the build lines (with optional
OpenMP, including macOS with Homebrew `libomp`) are at the top of the file.

`research/` holds the working notes, exact certificates and independent checks behind the
appendices; `research/README.md` says which of them the paper uses and which are historical.

Run scripts from the repository root:

```bash
.venv/bin/python scripts/<name>.py
```

The one exception is `scripts/enum2.py`, which reads and writes its pickle relative to the
current directory and must be run from inside `data/`.

## What produces what

| Paper | Script | Input / output in `data/` |
|---|---|---|
| Fig. 1 (schematic) | — (drawn by hand, `fig0_schematic.pdf`) | — |
| Table S1, Fig. 2 (coset counts to `t = 22`) | `enum_arith.py`, `arith_analysis.py`, `fig1_fig.py` | `enum_arith.pkl`, `enum_res2.pkl` |
| Table S2 (low-cost words per frame) | `lowcost_counts.py` | `lowcost_counts.json` |
| Table I (the bounds at `eps = 1e-10`) | `thm_numbers.py` | — |
| Fig. 3(a) (exact optima, both error budgets) | `frontier_eps_half.py`, `make_fig2.py` | `taumin_eps_half.json`, `frontier_L{5,6}_eps2.json` |
| Fig. 3(b), Table I (Ross–Selinger at chemistry scale) | `gridsynth_run.py` (20 000 angles at `eps`), then `gridsynth_run.py --eps 5e-11 --Q 7853981633 --out gridsynth_tau_epsover2.json` (the same angles at `eps/2`), then `make_fig2.py` | `gridsynth_tau.json`, `gridsynth_tau_epsover2.json`, `gridsynth_cache.json` |
| Sec. S1, exact `tau_min` for `L <= 6` | `taumin.py` | `taumin.pkl`, `taumin_exact.json` |
| Sec. S1, axis cosets: first scan (25 frames x 4 accuracies, `tau <= 15`) | `conj1_scan.py 15` | `conj1_scan_output.txt` |
| Sec. S1, axis cosets: second scan (272 frames x 5 accuracies, 1360 cells) | `scan2.py`, then `scan2_report.py` | `scan2.jsonl`, `scan2_summary.txt` |
| Sec. S1, required-`c` profiles to `t = 17` | `ht_profile.py`, `profile2.py` | `ht_profile_output.txt`, `profile2_{a,b}.txt` |
| Sec. S1, the identity-coset step (232 words, `c >= 163.9`) | `astra1.py` | `astra1_output.txt` |
| Sec. S4, all constants (Thm. S3, Cor. S2/S3, Prop. S2) | `mixing_bounds.py` | `mixing_bounds_output.txt` |
| Sec. S4, single-word `gridsynth` means at the mixing accuracies | `run_mixed_achievable.py` | `u_20000.txt`, `mixed_cache.json`, `mixed_achievable_results.json` |
| Lemma 8 and its `tau_0 >= 1` hypothesis | `check_lem_cost.py` | — (stdout) |
| Base enumeration (Matsumoto–Amano normal forms, resumable) | `cd data && ../.venv/bin/python ../scripts/enum2.py 0 22` | `enum_res2.pkl` |
| Thm. 9/10 (elementary tube bound, rate 11/5; proof in App. G): exact rational certificate | `research/cert/certificate.py 9/5` | `research/cert/certificate_output.txt` |
| Rate 2.23 with optimised box shape (remark after Thm. 10) | `research/cert/certificate_223.py 179/100 179/100 737/500 400/223` | `research/cert/certificate_223_output.txt` |
| Thm. 11/12 (height dichotomy, rate 17/7; proof in App. H, research notes D0–D7 in `research/dioph/DIOPH_NOTES.md`): exact rational certificate | `research/dioph/cert_dioph.py 28/17 28/17 157/100 11183/6800 400 100 12` (needs numpy; use `.venv/bin/python`) | `research/dioph/cert_17_7_output.txt` |
| Rate 12/5, same method, larger margin (research only, superseded by Thm. 4; needs the relaxed level-1 bound D1 of the notes) | `research/dioph/cert_dioph.py 5/3 5/3 31/20 41/25 400 100` | `research/dioph/cert_12_5_output.txt` |
| App. I (rational projections): exact certificate of the finite instance, rate 249/100 (proofs and review in `research/annulus/`) | `research/annulus/cert_proj.py 400/249 400/249 1192747/747000 319751/199200 400 100 12` | `research/annulus/cert_249_100_output.txt` |
| Thm. 4 (every rate below 5/2; proof in App. J, notes and review in `research/uniform52/`): continuum model and exact SMT check for all `0 < mu <= 1/10` | `research/uniform52/verify_z3.py 1/10 10/33` (needs `z3-solver`); model `research/uniform52/cont_model.py` | `research/uniform52/z3_outputs.txt`, `research/uniform52/z3_version.txt` |
| Remark S2 (ball example: 126 words of T-count 25 within radius 0.0025 of g/\|g\|, g = ((3+2√2)+i+j−k)/2, all on one sphere section; notes and review in `research/conditional3/`) | `research/conditional3/ballpts.cpp` (`ballpts 25 5.828427 1 1 -1 0.0025`, C++17, build line in the file; prints the stored 126 words, possibly in another order) and `research/conditional3/review/check_sections.py` (needs numpy) | `research/conditional3/review/bp_25_0.0025.txt`, `research/conditional3/review/known_output.txt` |

`scan2_report.py` recomputes, from the stored `scan2.jsonl`, every count quoted in the
axis-coset paragraph of Sec. S1 — the 272/1360/20400 totals, the 17 cells above `c = 48`,
the per-accuracy maxima, and the linear-term check. `scan2.jsonl` was produced by an earlier
version of the frame keys in `conj1_scan.py`, whose sign convention for a rotation axis
depended on floating-point noise: its 272 frames contain the 213 distinct axes of words of
cost `<= 3`, six of them twice with opposite signs (identical counts), and a sample of forty
cost-4/5 axes. The keys are now canonical and the samples are taken from sorted lists, so a
fresh `scan2.py` run is deterministic across machines; it samples a different set of forty
cost-4/5 axes, and the numbers in the paper are those of `scan2_report.py` on the stored file.

## Quick reproduction

Cheap (seconds to a couple of minutes, all from cached data):

```bash
.venv/bin/python scripts/thm_numbers.py      # Table I
.venv/bin/python scripts/mixing_bounds.py    # every constant of Sec. S4
.venv/bin/python scripts/scan2_report.py     # every count in the axis-coset paragraph
.venv/bin/python scripts/astra1.py           # the 232-word step forcing c >= 163.9
.venv/bin/python scripts/check_lem_cost.py   # Lemma 8 random-instance check
.venv/bin/python scripts/make_fig2.py        # Fig. 3
```

Expensive, and only needed to rebuild the cached data:

| | cost |
|---|---|
| `enum2.py 0 22` (from inside `data/`) | ~20 min and ~2.6 GB peak, single core; the `t = 22` layer alone is ~10 min. Each extra `t` doubles both. The times it prints are cumulative. |
| `scan2.py` | hours; resumable, appends to `data/scan2.jsonl` (see above) |
| `conj1_scan.py 15` | ~1 min |
| `lowcost_counts.py` | ~20 s |
| `ht_profile.py`, `profile2.py 0`, `profile2.py 1` | ~20 s each |
| `frontier_eps_half.py 6 24` | ~17 s, ~2 GB |
| `gridsynth_run.py` (both runs) | seconds with the stored `gridsynth_cache.json`; without it, 40 000 Ross–Selinger syntheses |
| `run_mixed_achievable.py --u-file data/u_20000.txt --backend pygridsynth` | seconds with the stored `mixed_cache.json` |

Times measured on a 10-core Apple-silicon laptop (Python 3.14, numpy 2.5).

## Determinism

Every randomised choice is seeded, so the outputs above are reproducible bit for bit
(timing lines excepted; the figures are written without matplotlib's version string, so
they are byte-identical too):

- Haar and perturbed frames: `numpy.random.default_rng(914)` in `conj1_scan.py`; `scan2.py`
  inherits the same generator by `exec`.
- The 20 000 grid angles of Fig. 3(b) and of the mixing measurement:
  `random.Random(20260903).randrange(7853981633)`, stored as `data/u_20000.txt`.
- `gridsynth_run.py`: `--seed 20260903` by default.
- `check_lem_cost.py`: `random.seed(916)`.

Two places where floating point used to decide a discrete outcome are now made exact:
`lowcost_counts.py` re-evaluates every word whose distance is within `1e-9` of `eps` with
50-digit `mpmath` arithmetic from its gate sequence and counts exact ties as within `eps`
(this changed one stored count, at `t = 16`, that no table uses, from 11772 to 11773); and
`conj1_scan.py`/`scan2.py` key rotation axes canonically and sample from sorted lists.

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
- Fig. 3(a) reports two error budgets. The certified `r = 2` point of Theorem 5 synthesizes
  each committed share at `eps/2`; evaluating at `eps` instead is a single-rotation
  calibration and not a point of the family. The two sequences approach 3 from opposite
  sides (2.478, 2.721, 3.113 against 3.333, 3.334, 3.205).

`data/enum_res2.pkl.bak-t19` is the `t <= 19` enumeration, kept for comparison.
