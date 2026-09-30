# High moments over frames (2026-09-30; revised after review/REVIEW_moments.md)

Aim: replace the worst-case frame statement needed for rate 3 (Conjecture H, or (★_{1/3})) by an *averaged*
statement over Haar-random cosets, and see how far it can be pushed.

Notation:
- Λ_τ is the set of words of T-count ≤ τ, with N_τ := |Λ_τ| = 72·2^τ − 48.
- R(W) ∈ SO(3) is the rotation of W.
- σ is normalised area on S², and s(ψ) := (1−cos ψ)/2 is the normalised area of a cap of angular radius ψ.

## 0. Hopf form of the tube count (exact)

Write x1 := G1·ẑ and x2 := G2^{-1}·ẑ, and φ(ε) := 4 arcsin(ε/2). Then
#(Λ_τ ∩ T_ε(G1R_zG2)) = N_φ(x1,x2) := #{W ∈ Λ_τ : ∠(R(W)x2, x1) ≤ φ(ε)}.

- Haar-random cosets correspond to independent uniform (x1, x2) ∈ S² × S².
- Mean: E_τ(φ) := ∫∫N_φ = N_τ s(φ) ≍ 2^τ ε².

## 1. Moments ⇒ worst case (rigorous)

**Lemma 1 (local constancy).** If ∠(x_i, x_i′) ≤ φ/2 for i = 1,2, then N_φ(x1,x2) ≤ N_{2φ}(x1′,x2′) (triangle
inequality, since R(W) is an isometry). Hence, for every k ≥ 1,

    s(φ/2)² · sup_{x1,x2} N_φ(x1,x2)^{2k} ≤ M_{2k}(2φ) := ∫∫ N_{2φ}^{2k} dσ dσ.

**Lemma 2 (tiny tubes, from App. E–F).** Suppose ε ≤ R'^{−4−η}, i.e. τ ≤ (1−η′)L. Then
n(τ) ≤ 2^{O(τ/log τ)}(1 + C R'^{2−a0/3}) ≤ 2^{max(0, τ/2 − L/3) + o(τ)}, where R'^4 ≍ 2^τ and ε = R'^{−a0}.
Lemmas E1, E5 and F2 are stated in App. F for a0 ∈ (0,2), but their proofs hold verbatim for a0 ≥ 4.

*Proof.*
1. **One hyperplane per arc.** Cover the core by 13 arcs of length 1/2 and run Lemma E1 with b = 0. The determinant
   has size ≤ C ε² R'^8 = C R'^{8−2a0}, which tends to 0. So the points of X_t on each arc lie in one K-hyperplane H.
2. **Volume.** The arc piece of the tube lies in a 4-box with sides R'/2, R'/32 + O(ε²R'), 4εR' and 4εR'. So
   vol₃ conv(Y∩piece) ≤ εR'³/8.
3. **Circles.** Lemma F2 with y = 0 gives ≤ 1 + C R'^{1+(3−a0)/3} circles.
4. **Points.** Lemma E5 gives ≤ 2^{O(τ/log τ)} points per circle. If the points of an arc span only a plane, they
   already lie on one circle and step 4 applies directly. ∎

For τ ≤ (1−η′)L this is ≤ 2^{τ/6+o(τ)} < 2^{τ/3}, so the tiny range needs no hypothesis.

**Hypothesis (M_k).** There is η(τ) = o(τ) such that M_{2k}(τ, φ) ≤ 2^{η(τ)}·(E + E^{2k}) for all τ and φ,
with E = E_τ(φ).

For comparison, the Poisson(E) value of the 2k-th moment is Σ_j S(2k, j) E^j ≍_k E + E^{2k}.

**Theorem M.** Assume (M_k) with k ≥ 4. Then every process satisfies T̄_t ≥ α_k[(1−2δ)m log₂K − 2m h₂(δ) − S − mA],
with α_k := 3 − 2/k and A = o(L).

*Proof.* Apply Lemma gibbs with λ = 1/α_k. We need n(τ) ≤ 2^{τ/α_k + o(τ)} for τ ≤ τ* = ⌈α_k log₂K⌉.
1. **τ ≤ (1−η′)L.** Lemma 2 applies.
2. **τ ≥ (1−η′)L.** By Lemma 1, and s(φ/2) = ε²/4 exactly, n(τ) ≤ (16 M_{2k}(2φ)/ε⁴)^{1/2k}. With E = E_τ(2φ) ≍ 2^{τ−2L} this gives
   n(τ) ≤ 2^{o(τ)}·ε^{−2/k}·max(E, E^{1/2k}).
   - **E ≥ 1.** We need 2L/k + τ − 2L ≤ τ/α. The worst case is τ = τ* ≤ αL + O(1), where the condition reads
     (α−1) ≤ 2(1 − 1/k), i.e. α ≤ 3 − 2/k.
   - **E < 1.** We need (τ + 2L)/(2k) ≤ τ/α for τ ≥ (1−η′)L. This holds iff (1−η′)(2k−α) ≥ 2α. At α = α_k it
     holds iff 2k² − 9k + 6 ≥ 0. The roots are ≈ 0.81 and 3.69, so this means k ≥ 4, and at k = 4 it allows
     η′ ≤ 1/11. ∎

Values: k = 4 gives 5/2, k = 5 gives 13/5, k = 10 gives 14/5, k = 20 gives 29/10, and k → ∞ gives 3.
- The main.tex rate is 17/7.
- Unconditionally, research/uniform52 gives every α < 5/2; this is not in the paper.
- The tube form of Conjecture H, via the four-grid covering, implies (M_k) and (U_k) for every k.

## 2. Reduction to the unit scale: submultiplicativity (rigorous)

**Proposition 3.** For τ = τ1 + τ2 and every φ,

    ‖N_{τ,φ}‖_{L^{2k}}/E_τ(φ) ≤ c_Λ · ‖N_{τ2,φ}‖_{L^{2k}}/E_{τ2}(φ),

with c_Λ := N_{τ1}N_{τ2}/N_τ ≤ 72.

*Proof.*
1. **Factorisation.** Every W ∈ Λ_τ factors as UV with U ∈ Λ_{τ1} and V ∈ Λ_{τ2}. Split any word of T-count
   ≤ τ; no normal form is needed. Moreover UV ∈ T_ε(C) iff V ∈ T_ε(U^{−1}C), by left invariance of d_proj. Hence
   N_τ(C) ≤ Σ_{U∈Λ_{τ1}} N_{τ2}(U^{−1}C).
2. **Minkowski.** Haar measure on cosets is invariant under C ↦ U^{−1}C. Minkowski's inequality then gives
   ‖N_τ‖_{2k} ≤ N_{τ1}‖N_{τ2}‖_{2k}. Finally divide by E_τ = N_τ s(φ). ∎

**Corollary (unit-scale hypothesis).** (M_k) at every level where E ≥ 1 follows from

    (U_k):  M_{2k}(τ, φ) ≤ 2^{η(τ)} E_τ(φ)  whenever E_τ(φ) ≤ 1.

*Proof.* Given (τ, φ) with E_τ(φ) ≥ 1, take τ2 := max{t ≤ τ : E_t(φ) ≤ 1} and apply Proposition 3.
- Since N_{t+1}/N_t ≤ 4, this gives E_{τ2} ∈ (1/4, 1].
- If E_0(φ) > 1, then φ is macroscopic and the bound is trivial.

So Theorem M holds under (U_k) for k ≥ 4. Its sub-unit levels use (U_k) directly, since M ≤ 2^{o}E is exactly what
case E < 1 needs.

**Reformulation of (U_k).** Write F_j := Σ_{distinct W1..Wj} σ×σ{(x1,x2) : all W_i in the tube} for the j-point
correlation. Since M_{2k} = Σ_j S(2k,j) F_j and F_1 = E, (U_k) is equivalent to F_j ≤ 2^{o(τ)} E for 2 ≤ j ≤ 2k whenever
E ≤ 1. For Poisson, F_j = E^j.

In words: in a random tube of expected count E ≤ 1, the expected number of ordered j-tuples of distinct words is
at most 2^{o(τ)}·E.

**What Proposition 3 cannot do by itself.** Propagation only goes *up* in τ at fixed radius. The regime where
Ramanujan gives κ := ‖N‖_{2k}/E ≤ 2 is the fat regime φ ≥ 2^{−τ2/4}, which lies above the range τ ≤ 3L. The
determinant bound is scale-covariant, so propagating it gives nothing new. Some independent input at the unit scale
is needed.

## 3. Spectral form (exact)

**Theorem S.** Let B_ℓ^U be an L²(dA)-orthonormal basis of the U-invariant (octahedral) spherical harmonics of
degree ℓ. Take it to be a common eigenbasis of the operators S_τ F(x) := Σ_{W∈Λ_τ} F(R(W)x).
- They are self-adjoint, because Λ_τ is inversion-invariant.
- They vanish off the U-invariants, because Λ_τ is U×U-invariant.
- They commute because of the tree recursion 1_{Γ1} * 1_{Γj} = 24·1_{Γj+1} + 48·1_{Γj−1}, where Γ_j is the set of
  words of T-count exactly j. This makes every S_τ a polynomial in S_1.
Write S_τ F = λ_F(τ) F. Then

    N_φ(x1,x2) = Σ_{ℓ≥0} ĉ_ℓ(φ) Σ_{F∈B_ℓ^U} λ_F(τ) F(x1) F(x2),   ĉ_ℓ(φ) = 2π ∫_{cos φ}^1 P_ℓ(t) dt.

*Proof.*
1. Funk–Hecke expansion: 1[x·y ≥ cos φ] = Σ_ℓ ĉ_ℓ Z_ℓ(x,y), with Z_ℓ the reproducing kernel of H_ℓ.
2. Σ_W Z_ℓ(R(W)x2, x1) is the kernel of S_τ restricted to H_ℓ.
3. S_τ vanishes off H_ℓ^U, because Λ_τ is invariant under right multiplication by U.
4. The ℓ = 0 term is N_τ s(φ) = E. ∎

**Consequence.**

    ∫∫(N − E)^{2k} = Σ_{(ℓ_i,F_i), ℓ_i ≥ 1} Π_i ĉ_{ℓ_i} λ_{F_i}(τ) · (∫_{S²} F_1⋯F_{2k} dσ)².

- The series needs a smooth majorant of the cap, or a summation method, to converge absolutely.
- The weights I² are squares. Sign cancellation can come from the products of Hecke eigenvalues and from the signs
  of ĉ_ℓ.
- **Ramanujan with no cancellation.** Using |λ| ≤ (τ+1)2^{τ/2} gives only sup N ≲ 2^{τ/2}. This is the
  square-root barrier again.
- **k = 1** is Plancherel. With (R) and the missing factor (4π)^{−2} restored, it gives Var ≤ Cτ²E. That is
  (U_1) under (R). A variance ≍ E is not shown.
- **k ≥ 2** involves correlations Π λ_{F_i} weighted by squared 2k-fold product integrals of Hecke-eigen octahedral
  harmonics. After Hecke multiplicativity and the trace formula, this is the j-point correlation count again.

## 4. Status of (U_k)

- **(U_1) holds under (R)**, via Plancherel (§3). The earlier ball-count sketch was incomplete: it needs ball counts
  up to radius 2^{−τ/3}. (U_1) gives only the trivial rate 1.
- **Along this route, rate > 5/2 needs (U_5)**, i.e. F_j for j ≤ 10.
- **What F_j counts.** For j = 3 it is Σ over pairs (W1, W2) of the number of third words near the *lattice coset*
  through W1 and W2. That coset is W1·T_V with V = W1^{−1}W2, and T_V is the centraliser torus of V, a CM torus
  K(V)^×.
  - So F_j counts (j−1)-tuples of quotients V_i ∈ W1^{−1}Λ_τ that nearly commute (nearly share an axis).
  - Exact commuting gives the linear clusters (P-units of K(V)), which are harmless.
  - Near-commuting at the unit scale is the open part.
  - Integrality forces exact structure only at much smaller scales:
    - v′×w′ gives exact parallelism only at about R'^{−8};
    - the numerators give R'^{−6} for an exact arithmetic coset (3-wedge);
    - they give R'^{−4} for a common hyperplane (4-determinant).

## 5. Numerics (`moments.cpp`; outputs `mom_tau*.txt`)

**Setup.**
- For each τ: 100 uniform x2, and 2·10⁵ uniform x1 per x2.
- Exact N at radii with E ∈ {1, 4, 16, 64}.
- Power sums of order ≤ 10.
- Control: N_τ i.i.d. uniform points, same code.

**Normalised raw moments.** Each entry is E[N^m]/Poisson_m(Ê), with Ê the empirical mean. In brackets: the
i.i.d. control. The last column gives the maximum count for words (control).

E = 1:
| τ | m=2 | m=4 | m=6 | m=8 | m=10 | max |
|---|---|---|---|---|---|---|
| 10 | 0.99 | 0.98 | 1.05 | 1.55 | 4.00 | 13 (10) |
| 12 | 1.00 | 1.01 | 1.03 | 1.11 | 1.37 | 12 (10) |
| 14 | 1.00 | 1.00 | 1.01 | 1.04 | 1.17 | 12 (10) |
| 15 | 1.00 | 1.00 | 1.06 | 1.87 | 10.78 | 20 (10) |
| 16 | 1.00 | 1.00 | 1.01 | 1.04 | 1.11 | 11 (9) |
| 18 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 10 (10) |
| 20 | 1.00 | 1.00 | 1.00 | 0.99 | 0.99 | 10 (10) |

E = 64:
| τ | m=2 | m=4 | m=6 | m=8 | m=10 | max |
|---|---|---|---|---|---|---|
| 10 | 1.00 | 1.01 | 1.03 | 1.06 | 1.11 | 133 (108) |
| 13 | 1.00 | 1.01 | 1.04 | 1.10 | 1.27 | 182 (114) |
| 16 | 1.00 | 1.00 | 1.00 | 1.01 | 1.01 | 128 (113) |
| 18 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 126 (110) |
| 20 | 1.00 | 1.00 | 1.00 | 0.99 | 0.99 | 123 (113) |

(E = 4 and 16 behave the same way. Full output: `python3 momtable.py`.) All control ratios are 1.00 ± 0.02.

**Reading (corrected after review; the earlier reading is withdrawn).**
- **Mechanism at small τ.** The deviations are *linear clusters*, not Clifford 8-groups.
  - At τ = 10, E = 1, a single x2 carries 76% of the order-10 moment. Its cap holds 13 words with T-counts 4…10.
    All share one rotation axis to 10⁻⁸ rad, at irrational angles: the powers of one hyperbolic element, of size
    2(τ−m)+1. This is why the histogram excess sits at 9, 11 and 13.
  - At τ = 15 one x2 carries 91% of the moment (a 21-word cluster). An independent x1 draw gives an order-10 ratio of
    16.4 instead of 10.78, so these ratios are not stable estimates.
- **The Clifford-axis runs.** They are correct as a statement about x2 near an axis: locally m₁₀ ≈ 5·10⁹. But no
  sampled x2 came within 2.4φ (τ = 10) or 48φ (τ = 15) of a Clifford axis. Their share of the sampled excess is about
  1%, so they do not explain it.
- **The flat ratios at τ ≥ 18 are uninformative.**
  - The HT axis alone contributes at least 0.15 of the Poisson value at τ = 10 and 0.008 at τ = 16, computed exactly.
  - The chance that any of the 100 sampled x2 lands near it is 7·10⁻⁴ (τ = 10) and 10⁻⁵ (τ = 16).
  - So "convergence to Poisson" is not supported. The sequence is not even monotone.
- **The regime Theorem M needs (E ≪ 1) was not sampled.** The E ≥ 1 runs are redundant, given Proposition 3.
- **What a meaningful test needs** (review):
  - per-x2 output with a bootstrap;
  - factorial moments at E ≤ 1;
  - importance sampling near the axes of short words;
  - exact sums over the linear clusters.

## 6. Correction: the moments are not needed; the content is a one-scale reduction (added after the review started; not independently reviewed)

**Theorem U (rigorous).** Let θ ∈ [0, 1/3]. Suppose that for every τ, every tube whose Haar-expected count is ≤ 1
contains at most 2^{θτ+o(τ)} words of Λ_τ. Then the exchange rate is ≥ 3 − 2θ, with A = o(L).

*Proof.* Fix the tube radius ε, and let τ₂ := max{t : E_t(ε) ≤ 1}. Then E_{τ₂} ∈ (1/4, 1] and τ₂ = 2L + O(1).
1. **τ ≤ τ₂.** The ε-tube lies inside the unit-scale tube of level τ (radius ≍ 2^{−τ/2} ≥ ε). So
   n(τ) ≤ 2^{θτ+o(τ)} ≤ 2^{τ/α}, since θ ≤ 1/3 ≤ 1/α.
2. **τ > τ₂.** By the factorisation of Proposition 3, n(τ) ≤ N_{τ−τ₂}·2^{θτ₂+o(τ)} = 2^{τ−2L+2θL+o(τ)}. At
   τ = αL this is ≤ 2^{L} iff α ≤ 3 − 2θ, and the condition is weaker for smaller τ.
3. **Conclude.** Apply Lemma gibbs as in Theorem M. ∎

**Consequences.**
- (U_k) implies the hypothesis of Theorem U with θ = 1/k. By Lemma 1 at the unit scale, sup N ≤ (2^{o}E/s(φ/2)²)^{1/2k}
  ≍ 2^{τ/k}. So Theorem M is a corollary of Theorem U, and the moment hypothesis is *stronger* than what is used.
- Averaging over frames does not escape the worst case: at the unit scale, local constancy turns a 2k-th moment bound
  into exactly the sup bound 2^{τ/k} that Theorem U needs.
- **Rate 3 therefore needs only "H at one scale": tubes with O(1) expected words hold 2^{o(τ)} words.** Conjecture H
  predicts O(τ) there, from the T-powers and the linear clusters.
- **Rate 3 in turn needs** n(2L) ≤ 2^{2L/3}, i.e. θ_unit ≤ 1/3 at τ = 2L.

**Unconditional status at the unit scale.**
- The continuum lemma model at a0 = 2 (`cont_model.py 0.4 0`) gives T = E1 + 0.133 = 1.466. That is n ≤ 2^{0.367τ} at
  the unit scale, slightly above the necessary 2^{τ/3}.
- Through Theorem U this yields only rate 2.27, below the direct 5/2. The determinant method is scale-covariant, so
  propagating it from the unit scale loses.

**Numerics.** After the review, the sampled moments (§5) support nothing beyond "no gross deviation for random
frames". Linear clusters dominate the high moments whenever they are hit, and they are rarely hit.
