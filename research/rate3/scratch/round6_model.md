# Round 6 scratch: black-box insufficiency model (critical scale, §10 normalisation)
N = 4^k points, ε = N^{-2/5}, vol = ε²N = N^{1/5}, cells 1/ε² = N^{4/5} (per band 1/ε = N^{2/5}), B = εN = N^{3/5}.
Determinant scale ε_det = N^{-3/5} = ε^{3/2}; √N = vol^{5/2}; 1/ε = vol².

## The five black boxes (conclusions actually used in §3–§16)
(R) Ramanujan Weyl sums: |Σ_x f(x) − N∫f| ≤ √N·k^C·‖f‖_{H^s} for all f ∈ C^∞(S²) (any frequency; Sobolev constant).
(D) App. E: the points of any cell of radius ε_det lie on one arithmetic plane section (circle) carrying ≤ 2^{o(k)} points.
(C) Divisor bounds: any arithmetic circle (plane with O_K coefficients of height ≤ 2^{Ck}) carries ≤ 2^{o(k)} points.
(L) Liouville: distinct points are ≥ N^{-1} apart (Bloch vectors in O_K³, norm 4^k, two embeddings).
(V) Volume law for bands about arithmetic axes of height H: count ≤ 2^{o(k)}(H²·ε·N + 1) for width ε.

## Model M_X
Take a Poisson-like configuration Π of N − M points (any realisation of N − M i.i.d. uniform points conditioned on (R),
which holds a.s. up to log factors), and add a spike: M = X·vol points placed in one ε-cell c₀, one per ε_det-subcell
(there are (ε/ε_det)² = 1/ε = vol² subcells, so X ≤ vol is possible), each on a distinct "arithmetic" circle.
Check at X = vol (spike M = vol² = 1/ε = N^{2/5}):
(R): spike contributes ≤ M·sup|f| ≤ N^{2/5}‖f‖_∞ ≤ √N‖f‖_∞. ✓ (Ramanujan never excludes spikes below √N.)
(D),(C): one point per ε_det-subcell, each its own circle. ✓
(L): spacing ≥ ε_det = N^{-3/5} ≥ N^{-1}. ✓
(V): any band of width ε containing c₀ has ≥ εN H² ≥ N^{3/5} ≥ M allowance. ✓ (M = N^{2/5} ≪ B = N^{3/5}.)
Also: in-band pair count P(R) − B·vol ≥ M² = N^{4/5} = B·vol — exactly the main term: the model saturates the trivial
pair bound, consistent with the fact that no pair asymptotic is known.
⇒ The five black boxes, combined by any inequalities (positivity, C–S, Hölder, Hecke recursion, interpolation),
   cannot give N_c ≤ vol^{2−δ}. Every proof of α > 5/2 must use a property FALSE in M_vol.

## What kills M_X (true-or-conjectural statements the model violates)
- k-point moment criterion: spike gives M^k vs allowance N^{4/5}·vol^k: M^k/(N^{4/5}vol^k) = X^k N^{-4/5}. With X = vol = N^{1/5}:
  k=4 ⇒ 1 (borderline: 4-point ⇔ 5/2), k=5 ⇒ N^{1/5} (excluded: 5-point ⇒ 13/5). Matches §15(b).
- Centred 4th moment: M⁴ vs N^{4/5}·3vol²: excluded for X > vol^{1/2}. Matches §14(a) (⇒ 3).
- In-band pair asymptotic with power saving: excluded for X > vol^{δ}. Matches §12(b).
- Nothing weaker (no Weyl-sum bound at any frequency, no plane/circle structure, no divisor bound, no band law) kills it.

## Sharpening: can M_X be made a *Hecke orbit-like* configuration? Properties of the true orbit NOT in the list:
(H1) exact Hecke recursion T₂F_k = F_{k+1} + F_k + 4F_{k−1} between levels (ties the configuration at level k to levels k±1);
(H2) exact unique factorisation of words (tree structure: each point has a unique depth-p prefix, 2^p classes, and
     ε-close points share ≤ 3t/5 of their t syllables (L1));
(H3) two-embedding joint equidistribution (σ₂ shadow) at the Ramanujan rate;
(H4) exact band counts B_W = Σ_{m∈W} r(m)r(2^k−m) and exact circle multiplicities r(·)r(·)/64 (arithmetic of Z[ζ₈]).
Does the spike survive (H1)? Propagation forward: at level k+s the spike sits in 1.5·4^s cells with ≥ M points,
volume 4^s vol: excess while 4^s < X. The L² allowance at level k+s for cells with excess ≥ M/2 is ≈ 4^s vol²·(stuff)/... :
same order as 1.5·4^s. Backward: computed in round 5 (1.137^s loss). So (H1) alone does not kill M_X: the model can be
extended to a consistent family {M_X^{(k+s)}} (spikes in all translates). (H1) is a consequence of (R)+tree structure.
(H3): the spike's M points have σ₂ shadows: joint equidistribution at the Ramanujan rate only says #(c₀ × c') ≤ ε²η²N + √N,
which is vacuous for the spike (M ≤ √N). Survives.
(H2): L1 is a statement about FRAMES: same-parity words with shared prefix c are ≥ c₀2^{-(t−c)/2} apart in SU(2). Frames in a
common ε-cube of the tube share ≤ t/5 syllables; the spike has one frame per cube, so (H2) imposes nothing on the spike.

## CORRECTION to (D) (after reading main.tex Thm rate52 / Rem. what52 / App. E)
The 5/2 machinery proves: #{words of T-count ≤ t′ in T_ε(C)} ≤ 2^{(2/5 − μ/11)t′ + o} for ε ≥ 2^{-(8/5+μ)t′/4}, μ ≤ 1/10.
This is "≤ 2^{o} words per δ-cube along the core", δ = 2^{-2t′/5}: a constraint on the FIBRE coordinate of the frames, flat in
the cap radius below δ. App. E (11/5) has an S² "cells on a circle" structure (Lemmas E3–E5) but with weaker exponents.
Consequences: (i) round-5 §5.4/5.2 is withdrawn (splitting loses 2^{3p/5}, as §3(e) says); (ii) the model M_X must satisfy the
corrected (D): spike frames spread ≤ 2^{o} per ε-arc of the fibre — M = 1/ε frames, one per arc: ✓ allowed. (iii) a rich cap at
the determinant level is exactly "every ε-arc of the Hopf fibre over z₀ carries one frame of T-count t".
## Verdict (round 6)
M_vol satisfies (R),(D corrected),(C),(L),(V),(H1),(H2),(H3). The five black boxes plus the Hecke recursion and the tree
structure, combined by any inequalities, cannot give N_c ≤ vol^{2−δ}. A proof must use a property false in M_vol: the only
candidates in sight are the exact arithmetic (H4) (band counts Σ r(m)r(n−m), circle multiplicities) at a level beyond
divisor bounds, or genuinely new cancellation (5-point correlation, in-band pair asymptotic, centred fourth moment).
