# Task (a): can counting cluster-bearing points along a tube push the rate past 5/2? (2026-09-30)

**Short answer: no.** Clusters are harmless for tubes. The configurations that bind the rate at 5/2 are not
clusters. They are about R'^{E1} sections, or boxes, carrying O(1) points each, where E1 = (8−2a0)/3. So controlling
clusters cannot move the rate.

Every rate above 5/2 is equivalent, up to R'^{o(1)}, to the statement that most boxes of length R'^{−E1} along
*every* tube are empty. That is a sparsity statement for isolated points at the δ-scale. No known method (geometry
of numbers, circle method, spectral method) reaches it, and each reduction tried below returns to it. The analysis
follows; nothing here touches main.tex.

## 1. The binding configuration at 5/2 (continuum model, `binding.py`)

In the uniform52 continuum model (a0 = 8/5 + μ, μ = 0.002, T = E1 + (10/33)μ), at the binding height η* = g − 2/5 + μ/11:

| class g | y = η* | worst piece: offsets | points per section | HIGH |
|---|---|---|---|---|
| 0.85 | 0.45 | 0.75 | **0.000** | E1 |
| 1.00 | 0.60 | 0.60 | **0.000** | E1 |
| 1.20 | 0.80 | 0.40 | **0.000** | E1 |
| 1.40 | 1.00 | 0.20 | **0.000** | E1 |
| 1.55 | 1.15 | 0.05 | **0.000** | E1 |

**The LOW side.** Every binding piece is a deep crossing class with points-per-section exponent 0. The total is
normals (route N/P) × offsets ≈ R'^{E1}.

**The HIGH side.** It is R'^b = R'^{E1} boxes at per-box cost R'^0.

**Conclusion.** At the limit the adversary is R'^{E1} sections, or boxes, with one point each. A cluster
(a section with R'^c points, c > 0) never binds. The model already charges clusters only once each, through the
global LOW count.

## 2. Clusters along a tube (heuristic accounting, consistent with §1)

The section clusters of NOTES.md §3 sit at Hecke points ξ_g of level l. Let s(l) = ρ2^{τ/2}/√H be the per-centre
size (H = N(l)). There are ≈ ρ²2^{τ/2}H^{3/2} successful l per dyadic H.

**At most O(1) cluster centres per successful l in the tube.**
- The total is Σ_H ρ²2^{τ/2}H^{3/2}·s = ρ³2^τ Σ H up to H ≤ ρ²2^τ.
- That is ≈ ρ⁵2^{2τ}, which is ε^{5−2α} at ρ = ε, τ = αL.
- This is ≤ ε^{−1} for all α ≤ 3.

**At most Hε² + O(√H) centres per l (Ramanujan at level l).**
- The total is ≈ ρ⁶2^{5τ/2} = ε^{6−5α/2}.
- This is ≤ ε^{−1} for α ≤ 14/5.

So even the Ramanujan-level worst case for clusters is harmless up to 14/5, beyond 5/2. Clusters are not the
obstruction. (The large-H end of these sums counts generic points once per section, so it overcounts.)

## 3. Why every "slot count" gives exactly E1

A slot is a potential location of one point. Each method below bounds the tube by (#slots) × (points per slot),
and each time #slots = R'^{E1}.

**(i) Boxes (Lemma E1).** Length R'^{−b} with 2a0 + 3b > 8, i.e. R'^{E1} boxes, each with ≤ R'^{o(1)} points at
the limit.

**(ii) Short differences.** This is a new reformulation.
- *Pigeonhole.* If the tube has N points, at least N/2 consecutive pairs have core gap ≤ 4π/N.
- *The differences.* Each d = w' − w lies in the body |P_Πσ₁d| ≲ R'/N, |P_⊥σ₁d| ≤ 2εR', |σ₂d| ≲ R'. By
  successive minima and wedge integrality (as in Lemma F5), this body holds ≲ R'^{8−2a0}/N² lattice vectors.
- *The count.* If every d carried ≤ R'^{o(1)} pairs, then N/2 ≤ R'^{8−2a0}/N², i.e. **N ≤ R'^{E1+o(1)}**.
- *Scope.* This is not a finished proof. Near-tangent d and low-height d need Lemma F2/F4-type care, and the crude
  F2 bound loses R'^{(2/3)(2−a0)}.
- *Where it stops.* It reaches E1 at best, again because every available d is charged.
- *Heuristic.* With the expected per-d count (points of Y_d in its tube window) in place of 1, it would give up to
  α ≈ 2.89. But Σ_d #(pairs on d) is by definition the number of close pairs, so this is circular.

**(iii) Sections (Lemmas F4, F5).** Normals × offsets meeting the tube at the binding depth give R'^{E1}.

**(iv) Degree-2 level 1 (quadric boxes, 12a + 16b > 44).** The ceiling is 28/11 if level 2 is free.
- *Exact thin tori are harmless, rigorously.* If some member of the pencil (nrd − n, Q) has rank-2 quadratic part,
  then its eigenvalue is a root of gcd(χ, χ′). So the eigenspace V₁ is defined over K or over a quadratic
  extension K′.
  - In the K′ case, w ↦ P_{V₁}w is a K-linear isomorphism K⁴ ≅ V₁ ⊗ K′, since V₁ ∩ K⁴ = 0.
  - Hence the torus {nrd = n, Q = 0} carries ≤ #{w₁ : nrd(w₁) = c₁} ≤ R'^{o(1)} points, by the divisor bound
    in a CM extension of K′.
- *Near-degenerate pencils.* An integrality/Liouville argument on disc(χ) should make low-height near-tori exactly
  degenerate. This is not worked out.
- *The loss.* Curved quadric pieces at level 2 lose; the earlier model (research/push) estimates about 2.28. So this
  route needs a full quadric analogue of App. F, with no guarantee of passing 5/2.

**(v) Hecke factorisation.** W = UV with U ∈ Λ_{τ₁} gives N(τ; C) = Σ_U N(τ−τ₁; U^{−1}C).
- Averaging over the Hecke orbit of frames and Ramanujan at level τ₁ gives an error ≥ 2^{τ₁/2} × (worst case at
  level τ−τ₁).
- That is ≥ 2^{2τ/5}. No gain.

**(vi) Arithmetic frames.** With G₂ ∈ Γ of height h, the tube count is a binary problem with a congruence condition
modulo an ideal of norm 2^h.
- The needed equidistribution is a sum of twisted CM theta coefficients.
- Deligne gives square-root cancellation only, i.e. rate 2.

## 4. What an argument beyond 5/2 must contain

- **The required statement.** For a0 < 8/5 (α > 5/2): along every tube, all but R'^{a0+o(1)} of the R'^{E1} boxes
  of length R'^{−E1} are empty. Equivalently, the number of close pairs (core gap ≤ R'^{−a0}) in every tube is
  ≤ R'^{a0+o(1)}.
- **What it is not.** It is not implied by any square-root input. Ramanujan, twisted Linnik and Lindelöf all stop
  at rate 2 (NOTES.md §2).
- **Why local counting fails.** Some boxes genuinely hold many points (NOTES.md §3–4), so the statement must be
  global along the coset.
- **Structure.** The problem splits as "clusters", which are harmless (§2), plus "isolated points", which are the
  whole difficulty. For isolated points one needs something like a *Poisson upper bound for close pairs along a
  great circle* for the points of O on the product of spheres.
- **Sharpness.** No example is known that makes 5/2 sharp. Proposition S would need dense tubes at T-count 5L/2, and
  the heuristics predict volume plus sparse clusters.
- **Numerics.** Enumeration cannot decide the question. The asymptotic regime needs 2^{L/2} ≫ 72·c ≈ 400, i.e.
  L ≳ 17 and τ ≳ 43 at T-count 5L/2.

## 5. Status

| statement | status |
|---|---|
| α ≥ 5/2 − o(1) | unconditional (uniform52) |
| clusters along tubes harmless up to 14/5 | heuristic accounting (§2) |
| binding configuration at 5/2 = one point per slot | verified in the continuum model (§1) |
| difference-pair argument reaches E1 | sketch (§3 ii) |
| exact thin quadric tori carry R'^{o(1)} points | proof sketch (§3 iv) |
| rate > 5/2 | open; needs sparsity of isolated points along tubes (§4) |
