# Conditional rate three: what can and cannot be assumed (2026-09-30, revised after review)

Task (i): derive α = 3 from a recognised number-theoretic conjecture, in the way Browning–Kumaraswamy–Steiner (BKS,
arXiv 1609.06097) derive Sarnak's optimal covering exponent for S³ from a twisted Linnik conjecture.

The independent review is `review/REVIEW.md`. Its fixes are applied below.

## Verdict

1. **Twisted Linnik (§2, heuristic).** In the natural implementation, the BKS delta method over K = Q(√2) under the
   K-analogue of their Conj. 1.1 gives the tube count only up to N^{1/2+o(1)} (N = 2^τ), the same as (R). That is
   rate 2, below the unconditional 5/2. For balls it gives (B_α) only for α ≤ 2, which is already elementary (§1).
   - Beating N^{1/2} needs cancellation *between* the dual frequencies, about one per δ-cube.
   - The candidates we see are bilinear Kloosterman sums via Kuznetsov / the spectral large sieve. They lead back to
     sums of Hecke eigenvalues, i.e. the spectral side of Remark barrier.
   - We see no standard input that does this without being equivalent to the count. That is weaker than "none
     exists".
2. **Ball covering (§3–§4).** A frame-free ball bound (B_α) at radius 2^{−τ/α} gives rate α by a covering argument
   (§1, rigorous). But, heuristically, ball maxima grow like 2^{τ(2/3−5/(3α))} for α > 5/2, because of sphere-section
   clusters around Hecke points of low height. So no argument that covers a tube by balls can give more than 5/2.
   - The mechanism is confirmed by exact arithmetic (review census, τ ≤ 26) and by exact lower bounds up to τ = 60
     (§4). At the mean spacing (μ = 1) there are balls with ≥ 98 words at τ = 50 and ≥ 124 at τ = 60. The Poisson
     baseline is ≈ 16–19.
   - This does **not** exclude other frame-free hypotheses, e.g. few sections per ball, height-dependent section
     bounds, or averaged ball counts.
3. **What remains (§5).** The weakest sufficient hypothesis we can state is the tube-level (★_{1/3}), which is implied
   by Conjecture H.

Nothing here touches main.tex.

---

## 1. Covering principle (rigorous; checked by the review)

Notation: Λ_τ = {W : t_min(W) ≤ τ} with |Λ_τ| = 72·2^τ − 48, and d = dproj.
- (PU(2), d) is RP³ with d([x],[y]) = min(|x−y|, |x+y|).
- The Haar measure of a d-ball is 2(θ − sin θ cos θ)/π with θ = 2 arcsin(ρ/2). This lies between 0.354ρ³ and
  0.424ρ³ for ρ ≤ √2.

**(B_α)**: max_ξ #(Λ_τ ∩ B(ξ, 2^{−τ/α})) ≤ 2^{η(τ)} with η(τ) = o(τ), for τ ≤ τ*.

**Lemma 1.** Assume (B_α), let T_ε(C) be the ε-tube of a coset C = {G₁R_z(θ)G₂}, and let 2^{−τ/α} ≥ cε. Then
n(τ) := #(Λ_τ ∩ T_ε(C)) ≤ c′·2^{η(τ)}·2^{τ/α}.

*Proof.* Put r = 2^{−τ/α} and θ_i = i·r for 0 ≤ i < 2π/r.
- d(R_zθ, R_zθ′) = 2|sin((θ−θ′)/4)| ≤ |θ−θ′|/2, and d is bi-invariant. So T_ε(C) is covered by the balls
  B(G₁R_z(θ_i)G₂, ε + r/4).
- Each such ball is covered by O(1 + (ε/r)³) = O(1) balls of radius r (volume bounds). If r ≳ 0.8 the covering is
  trivial. ∎

**Theorem 1.** If (B_α) holds with α ≤ 3, then every process as in Sec. 2 satisfies
T̄_t ≥ α[(1−2δ)m log₂K − 2m h₂(δ) − S − mA],
with A = 1 + log₂Z = O(log L) + max_{τ≤τ*} η(τ) = o(L).

*Proof.* Use Lemma gibbs with λ = 1/α and τ* = ⌈α log₂K⌉.
- From ε < sin(π/2Q) < π/2Q and K < Q: 2^{−τ/α} ≥ (2/π)2^{−1/α}ε for τ ≤ τ*. So Lemma 1 applies, and
  Z ≤ c(τ*+1)2^{max η}.
- The rest is the rate-2 bookkeeping with ½ replaced by 1/α, as in Theorem rate3seg. ∎

**Theorem 1′ (tube form).** Let (★_θ) be: n(τ) ≤ 2^{o(τ)} max(2^{θτ}, ε²2^τ) for all cosets and all τ ≤ τ*. Then
(★_θ) implies rate 1/θ for θ ≥ 1/3, because ε²2^{(1−θ)τ*} = O(1).
- **Relation to Conjecture H.** H counts *grid-restricted* words (∃ d ∈ Z_Q with Q ≤ Q_ε), so it does not bound
  n(τ) directly. Cover T_ε by four shifted grids with Q = Q_{2ε} at accuracy 2ε. This gives
  n(τ) ≤ 4(c₀ + c₁τ + 4c·2^τε²), hence (★_{1/3}).
  - Alternatively, state (★) for the grid alphabet, which is all that the Gibbs step uses (as in Theorem main2).
  - (★_{1/3}) is weaker than H: it is implied by H and not known to imply it. It gives the rate but not the
    quantization.
- **Unconditional ladder, each on its own τ-range:**
  - (R) gives θ = 1/2 (rate 2).
  - Theorems elem and dioph give θ = 0.4513 for τ ≤ (20/9)L and θ = 0.4111 for τ ≤ (17/7)(L+2).
  - research/uniform52 gives θ → 2/5 for τ ≤ (5/2 − μ)L.

**Calibration.**
- **Minimal spacing.** If W₁ ≠ W₂ ∈ Λ_τ, then d(W₁,W₂) ≥ 2^{−τ}/4.
  - Put y = ½trd(w₁w̄₂) ∈ ½O_K and x = n₁n₂ − y² ∈ ¼O_K.
  - x is totally ≥ 0, and x = 0 iff W₁ = W₂.
  - So σ₁(x) ≥ 1/(16σ₂(n₁n₂)), and d² ≥ 1 − c² = σ₁(x)/σ₁(n₁n₂) ≥ 1/(16 N(n₁n₂)).
  - This gives (B_1).
- **(B_α) for every α < 2 is unconditional** (review, fix 5).
  - After rescaling by powers of √2, one T-count parity lies on the sphere nrd = n_τ.
  - The 3-wedge of four points in a ball of radius ρ has σ₁·σ₂-size ≤ Cρ³R'⁶, with coordinates in ⅛O_K. So it
    vanishes for ρ < c·2^{−τ/2}.
  - The points then lie on a circle, and Lemma E5 gives 2^{O(τ/log τ)}.
  - So Theorem 1 yields every rate α < 2 without (R), and the twisted-Linnik ball statement (B_2) of §2 adds
    nothing.
- **Coplanarity.** The 5-point ball version of Lemma E1 puts all points of a ball of radius < c·2^{−2τ/5} = c·R'^{−8/5}
  on one K-hyperplane.

## 2. Twisted Linnik over Q(√2) (heuristic scaling)

**Setting.**
- w ∈ O, nrd(w) = n, N(n) = 2^τ, balanced by a unit so that σ₁(n) ≍ σ₂(n) ≍ R² with R⁴ ≍ N.
- σ₁(w) lies in a region of S³(R) and σ₂(w) is free.

The K-analogue of Heath-Brown's delta method (Browning–Vishe) gives
Σ_q Σ_c (Q₁Q₂)^{−2} N(q)^{−4} S_q(c) I_q(c).
- S_q(c) ≈ N(q)² S_K(n, F*(c); q) (BKS §3).
- I_q(c) is, up to the factor (Q₁Q₂)², the Fourier transform at c/q of the restricted surface measure.

**Error bookkeeping.**
- Under Linnik-strength cancellation in q (the K-analogue of BKS Conj. 1.1, with partial summation as in their
  (5.5)): E ≈ M·#𝒞_w/N(Q). This mnemonic reproduces BKS's εN^{1/2} (Lemma 4.1, #𝒞 = O(ε^{−1}N^{4δ}),
  Q = ε√N).
- The σ₂ weight is the whole shell, so Q₂ = R and the dual size is O(1). Thinner shells make it worse.
- **Ball of radius ρ.** Q₁ = ρR, so N(Q) = ρN^{1/2} and #𝒞 ≈ 1/ρ. Hence E ≈ ρN^{1/2}.
- **Tube of radius ε.** The torus weight has ⊥-width εR and in-plane radial width ε²R, so Q₁ = εR.
  - The annulus Fourier decay is (R|c_Π|/q)^{−1/2}. The weighted dual count Σ_{|c_Π|≤1/ε}(ε/|c_Π|)^{1/2} ≈ 1/ε
    comes from ≈ 1/ε² vectors of mean weight ε.
  - Hence E ≈ ε²N·ε^{−1}/(εN^{1/2}) = N^{1/2}, independent of ε. This agrees with Theorem tube.
  - The Bessel twist parameter |c_Π|/√F(c) ≤ 1 lies inside the range of Conj. 1.1.
- **Ternary reformulation (GLH).** It needs an arithmetic frame G₂ ∈ Γ of height h. Then v = w k′ w̄ has nrd of norm
  2^{2τ+h}, so there are about 2^{τ+h/2} points. Square-root cancellation gives 2^{τ/2+h/4} ≥ 2^{τ/2}.

**Consequences (heuristic).**
- In this framework twisted Linnik gives (★_{1/2}), i.e. rate 2, which (R) already gives.
- For balls it gives (B_2), which is elementary anyway (§1).

## 3. Sphere-section clusters (heuristic mechanism, exact examples)

Fix l ∈ O_K⁺ with N(l) ≈ H, g ∈ O primitive with nrd(g) = l, and s ∈ ½O_K.

**The section.** Y_g = {w ∈ O : nrd(w) = n, ⟨w,g⟩ = s}.
- Put m = nl − s². Then Y_g lies at d ≤ √(σ₁(m)/σ₁(nl)) from ξ_g = g/|g|.
- It is inside B(ξ_g, ρ) as soon as σ₁(m) ≤ ρ²σ₁(nl): "an event (l,s)".

**Frequency.** Events are lattice points of ½O_K in a box of area ≈ ρ²√N(nl).
- Heuristically they occur at rate ≈ 1.41ρ²2^{τ/2}N(l)^{1/2} per l.
- The expected number with N(l) ≤ H is ≈ 0.47ρ²2^{τ/2}H^{3/2}.
- The review's census confirms this: 36 events against 39 expected (μ = 16, τ = 20–26), 52 against 61 (μ = 4), and
  16 against 24 (μ = 1).
- Heuristically l and 4l (g ↦ 2g) give the same sections, so count primitive g only.
- Existence at H ≈ H_min is a random-model assumption: equidistribution of s² ≈ nl at σ₁ as l varies.

**Size.** z = w ḡ gives a bijection between pairs (w, g) on sections and the P-primitive z ∈ O with Re z = s and
nrd z = nl.
- The bijection is exact up to the 48 norm-one units when N(l) is odd. It is off by a factor 2–4 when P | l.
- So Σ_g |Y_g| = 48·#Z, and for prime l the ≈ 48σ(l) centres share Z uniformly. Here r_O(l) = 48Σ_{d|l}N(d),
  checked on 12 values of l.
- There are no local or spinor obstructions, since D is unramified at every finite place and h(O) = 1. Every
  census event had #Z > 0.
- Empirically #Z ≈ (18–42)·N(m)^{1/2} (median 29), and the per-centre size is ≈ #Z/σ(l).
- By Siegel (ineffective), #Z ≫ N(m)^{1/2−o(1)}.

**Optimisation.** Take the smallest H with an event, H_min ≈ (ρ^{−2}2^{−τ/2})^{2/3}. The largest cluster is then
≈ ρ2^{τ/2}/√H_min ≍ ρ^{5/3}2^{2τ/3}.
- With ρ = 2^{−τ/α} this is 2^{τ(2/3−5/(3α))} words, at height H ≍ 2^{τ(4/α−1)/3}.
- The exponent is positive iff α > 5/2. At α = 3 it is 2^{τ/9}.
- The threshold ρ ≥ 2^{−2τ/5} = R'^{−8/5} coincides with the 5-point coplanarity scale of §1. That is an arithmetic
  coincidence of the two thresholds, not an explanation.

**What this shows (heuristic).**
- M(ρ) := max_ξ #(Λ_τ ∩ B(ξ,ρ)) ≳ max(1, ρ³2^τ, ρ^{5/3}2^{2τ/3}).
- Hence every ball-covering bound n(τ) ≤ M(ρ)/ρ, for any ρ ≥ ε, is ≥ 2^{2τ/5}. Lemma-1-type arguments cannot beat 5/2.
- Rate > 5/2 needs most δ-cubes *along one tube* to be empty. A few heavy cubes do not contradict that, so this is
  a barrier for ball covering, not for the problem.

**Compatibility with Conjecture H (heuristic).**
- For tubes with ε ≥ r (r the section radius): the existence condition gives s/(2^τr²) ≤ r^{−1/3}2^{−τ/3} ≤ 1 for
  r ≥ 2^{−τ}.
- For ε < r: assuming Duke-type equidistribution on Y, a tube meets ≲ s(ε/r)^{3/2} ≤ c·2^τε² cluster words.
- For fixed l the centres are Hecke points of level l. A tube meets ≤ Hε² + O(√H) of them, so the dominant clusters
  contribute ≲ ρ2^{τ/2} = 2^{τ/6} ≪ 2^{τ/3}.
- Summing that worst case over all successful l is not enough, so a proof would need more (§5).

## 4. Numerics

**(a) Word-centred ball maxima (`balls.cpp`, `balls_words.txt`).** These are maxima over balls centred at words, not
over all centres.
- The code is correct (review).
- The first control (`--random`) used only the right-Clifford action. Λ_τ is invariant under a group of order 1152
  (left and right Cliffords, inversion). With a matched control (`review/balls_sym.cpp`, `review/symrandom.txt`) the
  word maxima exceed the control only at τ = 21 (μ = 16) and τ = 26 (μ = 4, 16).
- **So the table-level statistics are inconclusive.** The earlier readings ("exceed at every μ and grow",
  "repulsion") are withdrawn.

**(b) Exact sections (review census; `review/known_output.txt`, `review/census_*.txt`).** Word-centred balls see
only caps of sections, because the centres ξ_g are Hecke points and not words. Measured at ξ_g:
- **τ = 25, μ = 16 (ρ = 0.0025, Haar mean 16).** Take g = ((3+2√2)+i+j−k)/2, with nrd g = 5+3√2 and N = 7. The
  ball B(ξ_g, ρ) holds **126** words of T-count 25, all on ⟨w,g⟩ = (357+256√2)/2. `ballpts` confirms this
  independently.
- **Census for τ = 20…26 (H ≤ 40, μ = 16).** B(ξ_g, ρ) ≥ 64, 168, 41, 48, 48, 126, 48.
- **τ = 21, 168 words.** This is the exceptional l = 1 case (n nearly a square, N(m) = 42).
- **The earlier "τ = 26 cluster" was misidentified.**
  - The word-centred ball there caps two sections at 1.13ρ and 1.28ρ from ξ_g, of sizes 126 (T-count 25) and 219
    (T-count 26), inherited from n_τ·l coinciding across τ.
  - B(ξ_g, ρ) itself is empty at τ = 26.
  - The claimed "height 7.4 matched" was not a test.

**(c) Large τ, certified lower bounds (`clusters_arith.py`, `arith_tau*.txt`).**
- **Method.** For even τ, list every l with N(l) ≤ H_max and N(l) odd. Find every s with the whole section inside
  B(ξ_g, ρ_μ). Count the P-primitive z exactly. Then max_g |Y_g| ≥ 48#Z/r_O(l) = #Z/σ(l).
- These are exact lower bounds for max_ξ #(Λ_τ ∩ B(ξ, ρ_μ)), counting words of T-count exactly τ.
- **Baseline.** The i.i.d. Poisson maximum over M ≈ 72·2^τ/1152 symmetry-independent balls at μ = 1 (the quantile
  P(Pois(1) ≥ k) ≈ 1/M) is ≈ 11, 13, 16 and 18.5 at τ = 30, 40, 50 and 60.
- **Status.** τ = 30 used H_max 400, τ = 40 used 800, τ = 50 used 1500, all complete. τ = 60 (H_max 3000) is partial:
  6 of 53 events counted.

| τ  | μ=1     | μ=4   | μ=16  | l at μ=1: N(l), #Z/N(m)^{1/2} |
|----|---------|-------|-------|-------------------------------|
| 30 | 22.5    | 36    | 78    | 63, 39                        |
| 40 | 61.5    | 80.9  | 203   | 425, 49                       |
| 50 | 97.6    | 200.9 | 271.2 | 887, 32                       |
| 60 | ≥ 124.5 | ≥ 124.5 | ≥ 228.4 | 199, 10.5               |

**Reading.** At the mean spacing (μ = 1), balls around Hecke points ξ_g carry 22 → 62 → 98 → ≥ 124 words as τ goes
30 → 40 → 50 → 60.
- That is 2, 5, 6 and ≥ 6.7 times the Poisson baseline, growing roughly like 2^{τ/9}.
- The counts are exact. Only their interpretation as a growth law is heuristic.

## 5. What this means for α > 5/2

- In the natural circle and spectral implementations, square-root inputs (Ramanujan, twisted Linnik, Lindelöf via
  Waldspurger) give rate 2. Geometry of numbers already gives 5/2.
- Covering a tube by balls cannot beat 5/2 (heuristic §3, numerics §4b–c).
- A proof of α > 5/2 must be global along the coset. It has to show both:
  - most δ-cubes of one tube are empty;
  - the tube meets few cluster-bearing Hecke points.
  Ramanujan at level l, summed over all successful l, is not enough.
- The weakest hypothesis we can state that gives rate 3 is (★_{1/3}). It is tube-level and implied by Conjecture H
  (via the four-grid covering).

**Safe for a paper remark (review §5):**
- Theorems 1 and 1′ with the fixes, the minimal-spacing lemma, the elementary (B_{2−}), and 5-point coplanarity below
  c·2^{−2τ/5};
- the labelled heuristic "bounds that cover a tube by balls cannot give more than 5/2; this is why Conjecture H is
  stated for tubes";
- one exact example: 126 words of T-count 25 in the ball of radius 0.0025 (Haar mean 16) around g/|g|,
  g = ((3+2√2)+i+j−k)/2, all on ⟨w,g⟩ = (357+256√2)/2;
- the certified large-τ lower bounds of §4c.

**Not safe:**
- "no frame-free hypothesis can give more than 5/2";
- "(B_α) is false" stated without a heuristic label;
- "twisted Linnik cannot give more than rate 2" stated without "in this framework".

Files: `balls.cpp`, `ballpts.cpp`, `clusan.py`, `clusters_arith.py`; `balls_words.txt`, `balls_random.txt`,
`clus21.txt`, `clus26.txt`, `arith_tau{30,40,50,60}.txt`; `review/`.
