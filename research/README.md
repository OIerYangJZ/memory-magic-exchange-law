# research/ — working material behind the appendices

The paper (`main.tex`, with the Supplemental Material `supplement.tex`) is the authoritative
statement of every result. This directory keeps the working notes, exact certificates and
independent checks from which the appendices were written. The notes are dated working
documents: their theorem, section and equation numbers refer to the `main.tex` of their date,
and some of their intermediate claims were later superseded. For example `NOTES.md`
(2026-09-28) still says that the unconditional rate in the general CW model stops at 2; the
paper now proves every rate below 5/2 (Theorem 4). Review notes occasionally mention scratch
files on the authors' machines (`school-server:/tmp/...`, `/tmp/main_before_rate52.tex`,
`qip/extended_abstract.tex`) that are not part of the repository.

## Used by the paper

The exact commands and expected outputs are in the table of `../README.md`.

| Path | Paper | Content |
|---|---|---|
| `cert/` | Thms. 9–10, App. G | `certificate.py`, `certificate_223.py`: exact rational certificates for the rates 11/5 and 2.23 |
| `dioph/` | Thms. 11–12, App. H | `cert_dioph.py`: exact certificate for 17/7; `DIOPH_NOTES.md`: the lemmas D0–D7; `fast.py`: float optimiser used to find the parameters |
| `annulus/` | App. I | `cert_proj.py`: exact certificate of the finite instance (rate 249/100); model and probe scripts; `review/` |
| `uniform52/` | Thm. 4, App. J | `cont_model.py`: continuum form of the case analysis; `verify_z3.py`: exact SMT check; `review/` |
| `verify4/indep_model.py` | App. H | independent float re-implementation of the D1–D7 model (the value 1.6347 quoted there) |
| `conditional3/` | Remark S2 | `ballpts.cpp`: the ball example (build line at the top of the file); notes on the conditional rate-three programme |

## Independent checks

`verify/`, `verify2/`, `verify3/`, `verify4/`, `verify5/`: re-implementations and spot checks
of the determinant-method bookkeeping, of the certificates against the lemmas as written in
the paper, of Theorem 6 and Corollary 2, and of the numbers of Proposition S1. Each script's
docstring says what it checks.

## Historical explorations

- `NOTES.md` (2026-09-28): results A, A1, B and D there became Theorem 6, Corollary 2,
  Theorem 7 and Proposition S1; result C (one-sided cosets) is not used in the paper; the
  summary's "the unconditional rate stays at 2" is superseded by Theorem 4.
  `zomega.py` (exact Z[ω] arithmetic) and `check_identity.py`, `check_axis.py`, `check_lde.py`
  are the checks of those notes; `batch_qrom_cost.py` is the cost model of result D.
- `push/`, `push2/`: attempts to go beyond the rates 2.22 and 2.23 with the first lemma set,
  superseded by Theorems 11 and 4.
- `ternary/`: small-cap equidistribution on ternary Z[√2]-spheres (not used in the paper).
- `rigid_scheme.py`, `rigid_scheme_ab.py`, `rigid_plane.py`, `rigid_AB.py`, `certify_plane.py`:
  early versions of the determinant-method bookkeeping, superseded by `cert/` and `dioph/`.
