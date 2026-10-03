# research/rate3: handoff (2026-10-02)

The goal is to raise the paper's unconditional exchange rate (arXiv:2609.37368; v2 proves every α < 5/2) to α = 3,
or to anything above 5/2.

**Status: not reached. The unconditional rate is still 5/2.** All reasoning is in `NOTES.md` §1–§15 and the numerics
are in `numerics/RESULTS.md`. This file is a summary and an index. Read it first; then go to NOTES for proofs and
exact statements. Per-round work log (from 2026-10-02 on): `PROGRESS.md`.

## 1. Setting (one paragraph)

K = Q(√2), and O is the maximal order of (−1,−1)_K. Γ is Clifford+T modulo phase, acting on the 3-regular tree at
P = (√2); the T-count is the tree distance.

Single-qubit states Γ|0⟩ ↔ points (X, Y, Z) ∈ O_K³ on X² + Y² + Z² = 4^k (sde k, T-count τ ≈ 2k), counted mod phase.

Rate α needs, at scale ε = 2^{−L}, a bound n(τ) ≤ 2^{τ/α+o} for all τ < αL (NOTES §1). At the 5/2 critical scale
ε = 2^{−2τ/5}:
- the volume is 2^{τ/5};
- the cell count is 2^{2τ/5};
- the determinant method (paper, App. E) gives exactly one point per cell.

## 2. Proven (rigorous) — pointers

| Result | Where |
|---|---|
| For two rounds only Hopf fibres matter: the r = 2 rate is about Clifford+T states in Bloch caps (SA) | §2 |
| Exact two-round characterisation (covering by translates): r = 2 rate = cell-occupancy exponent, so α_{r=2} = 3 ⟺ SA | §7 |
| Hence any proof of α = 3 must prove SA; no modelling loophole | §7, §12(a) |
| Rounds ≥ 2 reduce to Hecke-ball tubes Λ_τ·V\|0⟩. Arithmetic frames suffice, uniformly in height | §2, §8(a) |
| Band decomposition: z-bands obey the volume law (divisor bound) | §10(A) |
| Caps reduce to short shifted convolutions S_ℓ of CM Hilbert theta series; \|S_ℓ\| ≤ L^{1−θ} ⇒ rate 1 + 1/(1−θ) | §10(A) |
| Whole-sphere pair correlation of states is Poissonian with a power saving at all scales above the mean spacing (Ramanujan; it is a sum of squares) | §12(b) |
| Variance route: in a band B·vol = cells², so a power-saving in-band pair-count asymptotic ⇒ α > 5/2 (r = 2) | §12(b) |
| θ > 1/3 in mean square over ℓ already suffices | §12(b) |
| Moment criteria (r = 2): raw 2k-point correlation upper bounds (2^{o}·Poisson) ⇒ rate 3 − 1/k; a centred fourth moment ⇒ 3 | §13(a), §14(a) |
| Averaged over all norms ≍ 2^τ, the tube volume law holds on [2L, 3L] (Marshall Prop. 5.2 transplanted) | §9 (sketch only, not written up) |

Conditional results: Conjecture H ⇒ α = 3 (in the paper), and SA ⇒ α = 3 for r = 2.

## 3. Numerical evidence (`numerics/`; server copy in `~/mmx_claude/rate3num`)

| Test | Range | Result |
|---|---|---|
| A: exact S_ℓ, all bands and all ℓ ≤ 1/ε | sde k ≤ 16 | Z = S_ℓ/(8√S₀) is k-independent (99.99% quantile ≈ 4.5); θ = 1/2 up to logarithms |
| B: r = 2 caps at the critical scale | k ≤ 16 | Poisson. k = 16: max 166, Poisson max 168, volume 107, determinant level ≈ 1.1·10⁴ |
| C: Hecke-ball tubes (rounds ≥ 3), all words of T-count ≤ t applied to V\|0⟩ | t ≤ 30 | For source T-count h ≥ 8 the maxima equal the Poisson maxima (t = 30: 478–485 against 482) |
| Special points (`caps.c`) | k ≤ 16 | Divisor-type "rich circles" at \|+⟩ and the H-eigenstates: ≤ 5.5× volume at k = 12, ≤ 1.4× for k ≥ 14 |
| D: factorial moments of cell counts, orders 2–8, both scales (`moments.py`) | k ≤ 16 | Poisson to three decimals for k ≥ 12; centred 4th moment within 1–7% of 3λ²+λ |

There is no numerical sign that α = 3 is false.

## 4. Open targets, in order of apparent tractability

Each line gives the conclusion for r = 2. Rounds ≥ 3 need the same statement for Hecke-ball sources V|0⟩, uniformly
in the height of V.
1. 5-point correlation upper bound for states in ε-caps, at most 2^{o} × Poisson ⇒ **13/5** (§15(b)); 6-point ⇒ **8/3** (§13(a)). The 4-point bound (⇒ 5/2) is itself unproved.
2. Power-saving asymptotic for the in-band pair count ⇒ **> 5/2** (§12(b)).
3. θ > 1/3 in mean square over ℓ for S_ℓ ⇒ **> 5/2** (§10(A), §12(b)).
4. Centred fourth moment of cap counts ⇒ **3** (§14(a)).

All four are beyond-square-root statements for a single Hecke orbit at a single prime. Ramanujan gives the
two-point (L²) level and nothing more. The rational toy model, Zhang–Zhu's tube conjecture on λS³, is open at the
same exponent (§8(d)).

## 5. Routes tried and where each stops (index)

| Route | Stops at | Where |
|---|---|---|
| Dirichlet fibration / hybrid sums S_χ | Burgess / Pólya–Vinogradov endpoint (dual length q^{1/4}): 5/2 | §3(a), §5(iv), agent_burgess |
| Quotients W_jW_i⁻¹, multiplicative energy | rate 2 | §3(b), §9 |
| Hecke-orbit bootstrapping | square-root barrier | §3(c) |
| Degree-D auxiliary forms; one global auxiliary polynomial (Siegel) | worse than E1; infeasible for a₀ < 3.5 | §3(d), §9 |
| 2-adic determinant; complex 3-point determinant | same problem at lower level; weaker than the 5-point K-determinant | §3(e)(f) |
| Fourier in the fibre angle (single a) | needs cancellation over a | §3(g) |
| Short differences; gap principle; Plücker for six points | no counting contradiction (Cayley–Menger) | §5(i)(ii) |
| Bilinear factorisation W = W₁W₂; Duke-type quotients | 8/3 only under an unproved fourth moment; else 2 | §5(iii) |
| Rigid 1-dimensional analogues (Cilleruelo–Córdoba, three-gap) | no non-abelian analogue | §5(v), §7, §10(C) |
| Sup-norm literature (Marshall, Blomer–Michel, Khayutin–Nelson–Steiner) | single-norm bounds are geometry of numbers only; one prime means no averaging | §6 |
| Process-side loopholes; subspace theorem; expander mixing; Hecke-stable profiles | none / a₀ > 4 / worse than Ramanujan / impossible | §7 |
| GRH, Lindelöf, Ramanujan, twisted Linnik, density hypothesis | rate ≤ 2 for tube upper bounds | §8(b), §9 |
| Decoupling and polynomial methods | window barrier: Archimedean methods see a norm window, one point per cell | §8(f), agent_decoupling |
| Linnik basic lemma, entropy, flattening, Bourgain–Gamburd | L²/ball inputs are blind to one tube (rigorous no-go) | §8(f), agent_linnik, §13(c) |
| Auxiliary-prime amplification; van der Corput (non-amenable); products of tube words; local spectral measures; random axes | as stated in §9 | §9 |
| Arithmetic frames, 2-adic depth (Postnikov) | PV size 2^L = cells reappears; Postnikov does not linearise a_χ | §9 |
| Excess invariance / transport; exact σ₂-sphere; Pancharatnam; double caps; sum-product; approximate APs; meet-in-the-middle; CM tori | as stated in §10(B)–(D) | §10 |
| Variance in a band (spectral) | Ramanujan only, Var_R ≲ N: localisation fails | §12(b) |
| Difference-vector family (averages over ~2^{2.4k} norms) | class-group characters: the diagonal h_M dominates; needs family cancellation | §13(b) |
| Residue classes mod P^τ (midpoints), Hecke-orbit unions, recursive splitting | localisation barrier | §13(c) |
| Norm-averaged volume law made precise; sub-family {P^τl}; positivity; moments over norms | weight 2^{−τ}; orbit averaging is neutral; loses 2^τ | §15(a) |
| 3-point correlation: prefix/slope/tree-distance/operator-norm forms; joint σ₁-cap × ∂tree-ball form | all return to depth-τ caps at depth-τ arithmetic points | §15(c) |
| Rarity + propagation; determinant + L² interpolation; two-embedding determinant; additive energy; circle method; genus-2 Siegel / Saito–Kurokawa; sieves; incidences; composite-norm lift; additive perturbation; random targets | as stated in §14(c) | §14(c) |

Spectral form of the obstruction (§14(b)):
- The family of forms has size F ≈ ε⁻². At the rate-α scale the twist λ_f(P^{2k}) has length F^{α/2}, so α = 2
  corresponds to F, 5/2 to F^{5/4} and 3 to F^{3/2}.
- Petersson over K has no decay at σ₂, where the weight is 2.

## 5b. Cloud session 2026-10-02 (rounds 4–7; NOTES §16, PROGRESS.md §4–§7)

- Positivity transfers from the §15(d) sum of squares all return the Ramanujan sup bound (§16(a)).
- P-adic Fourier form of the class count; exact Hecke recursion of deviation fields, Poisson growth saturates Ramanujan (§16(b)).
- **Black-box insufficiency model (§16(c)):** every input used so far (Ramanujan Weyl sums, the 5/2 tube bound, divisor bounds,
  Liouville/L1, band law, Hecke recursion, σ₂ shadow) is satisfied by a configuration with a vol² spike. A proof needs a new input.
- SA ⟺ small-scale equidistribution of signed prime-Hecke-angle sums over the shifted pair (m, n−m) (§16(d)).
- Two claims made during the session were withdrawn in round 6 (PROGRESS §6 勘误): the 5/2 tube bound is flat below the critical
  cube, so prefix splitting loses 2^{3p/5} (confirming §3(e)).

## 5c. Cloud session 2026-10-02, rounds 9–12 (NOTES §16(f)–(i), PROGRESS.md §9–§12)

- The K-determinant constrains only frames within one fibre ε-arc; cross-arc 5-tuples (the content of the 5-point criterion) carry no
  algebraic relation at the critical scale (§16(f)).
- Eleven true properties of the orbit are all absorbed by the spike model; a killer must be single-norm, single-cap, sub-√N (§16(g)).
- The band count B_W has no provable asymptotic at the critical scale; numerically B_W = main(1 + O(|W|^{−1/2})) (Test E, §16(h)).
- **Minimal open problem (M₀)** (§16(i)): Σ_{j≤N^{2/5}}|Σ_x Y_{j0}(ẑ·x)|² ≤ N^{8/5−δ} — zonal Weyl sums of one orbit about the Clifford
  axis; implies the band asymptotic; spectrally it is vertical Sato–Tate at scale 1/log N. Use it as a test of technique: (M₀) is a diagnostic. It neither implies nor is implied by SA (SA needs only the band upper bound), and proving it would not move the rate (PROGRESS round 13).

## 6. Corrections inside NOTES (read the corrected versions)

- The §6 side remark on multi-prime gate sets is corrected in §9.
- In §8(d), the Bourgain–Rudnick bound is n^{3/5}, not n^{7/10}.
- §8(c) ("numerics cannot reach the regime") was too pessimistic; corrected in §12(c).
- The §11 verdict (pointwise θ > 1/3) is corrected in §12(b): mean square suffices.

## 7. Compute and files

- **Server.** `ssh school-server`, code in `~/mmx_claude/rate3num`.
  - Use `nice`; about 26 GB of disk is free.
  - k = 16 in `sl.c` needs about 24 GB of RAM. A t = 30 run of `orbit.c` takes about 9 min on 20 threads.
  - The laptop has no OpenMP, so compile on the server: `gcc -O3 -march=native -fopenmp`.
- **Code.**
  - `sl.c`: S_ℓ and cell counts for states (Tests A, B).
  - `caps.c`: caps around fixed centres.
  - `orbit.c`: Hecke-ball orbits (Test C).
  - `brute.py`: brute-force check.
  - `analyze.py`: statistics.
- **Earlier related work.**
  - `research/conditional3/`: Theorem 1 (H ⇒ 3), moments, Theorem U.
  - `research/uniform52/`: the 5/2 proof, continuum reduction plus z3.
  - `research/annulus/`: 249/100.
