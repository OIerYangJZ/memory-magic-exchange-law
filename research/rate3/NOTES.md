# Attempt at α = 3 (2026-10-01)

**Verdict (updated 2026-10-02, through §14): α = 3 is not reached. No exponent beyond 5/2 is obtained.** This note
records the new reductions and the routes tried, and says where each one stops. Start with `HANDOFF.md` for a
summary and an index. Status and open targets are at the end of §14. Nothing here touches main.tex. Earlier work is in
`research/conditional3/` (NOTES, NOTES_a, NOTES_moments) and `research/uniform52/`.

## 1. What rate α > 5/2 needs (recap)

The Gibbs step only needs the number of grid cells along the core that are reachable at cost ≤ τ, because A ↦ W is
injective into cells given the context. That count is trivially ≤ 2^L. So rate α needs
n(τ) ≤ 2^{τ/α + o(L)} only for τ < αL.

- **Where it binds.** For α = 5/2 + δ, the binding T-count is τ ≈ 5L/2. There the determinant method gives exactly
  2^L, one point per ε-cell.
- **What is needed there.** n(5L/2) ≤ 2^{(1 − 2δ/5 + O(δ²))L}, i.e. a power-fraction of the ε-cells along *every*
  tube must be empty. Volume predicts 2^{L/2}.
- **Rate 3.** At τ = 2.9L the requirement is 2^{0.967L}, against a volume of 2^{0.9L}. As α → 3, only a vanishing
  power of slack over the volume law remains.

## 2. Two rounds: only Hopf fibres matter (new reduction)

In round 1 there is no prefix, so G₁ = I. For r = 2 this is the only committed round. So the two-round exchange
rate depends only on the cosets {R_z(θ)G₂}.

**State form.**
- W ≈ R_z(θ)G₂ iff W⁻¹|0⟩ ≈ e^{iθ/2}G₂⁻¹|0⟩. So these tubes are the ε-neighbourhoods of Hopf fibres, i.e. of right
  complex lines of the Z[ω]-module O.
- W⁻¹ is determined by W⁻¹|0⟩ up to right multiplication by the 8 diagonal elements T^j of Γ.
- Hence #(Λ_τ ∩ tube) ≤ 8·#{Clifford+T states ψ of T-count ≤ τ whose Bloch vector lies in a cap of radius cε}.
- The core position is the phase of ψ, taken mod π/4.

**Conjecture SA (state approximation).** Every single-qubit state has at most 2^{o(τ)}(1 + 2^τ ε²)
ε-approximations among Clifford+T states of T-count ≤ τ.

**Consequences.**
- SA implies α = 3 for r = 2, by Theorem 1 of conditional3/NOTES.md.
- Conversely, if SA fails with phases spread along the fibre, the process of Prop. S beats 3 for r = 2. The
  phase-spread caveat is the one in Remark what52.

**Rounds t ≥ 2.** These have an arithmetic prefix P of T-count h. The relevant set is Λ_τ·φ with φ = P⁻¹|0⟩, which
is state-to-state approximation, i.e. Conjecture H in Hopf form.
- SA for |0⟩ alone does not suffice: applying it to PW ∈ Λ_{τ+h} loses a factor 2^h.
- So for general r, nothing weaker than H (tube form) is visible.

**Classical analogue.** SA is the S-arithmetic analogue (S = {√2}) of the Bourgain–Rudnick cap conjecture for
x² + y² + z². The literature status for Z³:
- *Individual caps:* no bound better than the Jarník-type R^ε(1 + r²/R^{1/2}) of Bourgain–Rudnick (2012) is known.
- *Almost all caps:* Humphries–Radziwiłł (arXiv 1910.01360) and Lutsko (arXiv 2402.12822, averaged over n).
- *Under GRH:* Shubin (arXiv 2108.00726) gets the variance of the conjectured order.

The adversary picks the frame G₂ freely, because the final round is free. So almost-all and variance results do not
help: local constancy turns them back into the worst case (NOTES_moments §6).

## 3. Routes tried; where each stops

**(a) Dirichlet fibration on a Hopf fibre.**
- *Construction.* Take t ≈ z₀u. Choose p/q ∈ Q(ω) by Minkowski with |q|, |p| ≤ X in both embeddings and
  |σ₁(qz₀ − p)| ≤ cX^{−3}. Then the tube points fibre over g = qt − pu ∈ Z[ω].
- *Fibres.* Each fibre is a circle parallel to the core. It carries ≤ 2^{o(τ)} points (norm equation in Q(ω)).
- *Count.* #g ≍ R'^4(X²ε + X^{−2})² ≥ 4εR'^4. This gives rate 2 alone.
- *Where the difficulty went.* In the coordinates (⟨w,v⟩, Δ(w,v)) this is the Clifford-frame count on a sublattice
  of index N(H)². The g-projection is onto, and the whole difficulty sits in a congruence mod H in the fibre: a
  shifted-convolution sum, which is square-root at best.

**(b) Quotients and energy.**
- W_jW_i⁻¹ lies in the tube of the *standard* torus, where the volume law is proved (Thm. clifford), but at T-count
  2τ.
- Even with minimal multiplicative energy this gives N² ≤ 2^{o}ε²4^τ, i.e. rate 2.

**(c) Hecke-orbit bootstrapping.**
- Use n_τ(x) ≤ Σ_{σ∈Λ_j} n_{τ−j}(R(σ)x). Split g = n_{τ−j}^{(2φ)} at level λ, and use Ramanujan on min(g, λ) with
  a φ-ball average.
- The error is j·λ^{1/2}·2^{τ/2}, independent of j. This is the square-root barrier again.

**(d) Degree-D auxiliary forms on the quadric** (leading-order count).
- Monomials c^a t₁^b t₂^e r^f with f ≤ 1, where the quadric removes r².
- The leading term requires b > 6 − 2a₀, and σ₂ patches do not change this. That is worse than (8 − 2a₀)/3 for
  a₀ < 5/2.
- The reason: the curvature gain enters only through monomials linear in the normal coordinate. This is consistent
  with the last paragraph of Remark what52.

**(e) 2-adic determinant.**
- Points sharing a suffix σ of length s have their 4×4 determinant in nrd(σ)²O_K.
- This is the same as factoring σ off, so it is the same problem at T-count τ − s for 2^s frames, and it loses
  2^{s/3}.

**(f) Complex determinant for Hopf fibres.**
- The 3-point Δ over Z[ω] vanishes only for ℓ < 1/(8εR'^4).
- Its size is ℓ²ε²R'^8, against ℓ³ε²R'^8 for the 5-point K-determinant, so it is weaker.

**(g) Fourier analysis in the fibre angle.**
- *Expansion.* The count is Σ_{|a|≤1/ε} ĉ_a e^{−ia·arg z₀} S_a(2^k), where S_a are shifted convolution sums of the
  coefficients of the CM theta series of the Hecke characters (u/|u|)^a (a ≡ 0 mod 8) of Q(ω).
- *Single a.* Each S_a with a ≠ 0 is at least 2^{τ/2} under square-root cancellation.
- *Over a.* Beating that needs cancellation over the weight a, which is the count itself.

## 4. Status

- **Unconditional:** α ≥ 5/2 − o(1) is unchanged.
- **Necessary for any α > 5/2** (already the conclusion of NOTES_a): sparsity of isolated points along every tube
  at T-count 5L/2. For r = 2 this is a worst-case statement about Clifford+T *states* in Bloch caps.
- **Classical analogue:** in the Z³ analogue, no worst-case cap bound of this kind is known beyond the elementary
  Jarník range.
- **Conditional:**
  - Conjecture H gives α = 3 for all r (paper).
  - SA (states only) gives α = 3 for r = 2 (§2).

## 5. Second attempt: sparsity of isolated points along tubes (2026-10-01, later)

**Target.** It is enough to beat the uniform52 exponent at the critical point by any fixed amount. Precisely: there
should be κ, η > 0 such that for all a₀ ∈ [8/5 − η, 8/5 + η] and all tubes,
#(tube points) ≤ R'^{8/5 − κ + o(1)}. That gives α = 5/2 + c(κ, η) > 5/2. **Not achieved.**

What was found, in order.

**(i) Short differences fill the lattice exactly.**
- Consecutive points of a fully occupied tube differ by d ∈ O with |σ₁d| ≲ δ in *all* four directions and
  |σ₂d| ≲ R'. The body is a product of balls, so #(O ∩ B) ≍ Σ_{N(m) ≲ R'^{4/5}} r_O(m) ≍ R'^{8/5}, which equals the
  number of cells.
- Each d is realised by O(1) pairs. The direction of P_Π σ₁d fixes the cell up to the radial wobble ε²R', which is
  sub-cell.
- So "one point per cell" is consistent with the lattice at every gap J: pairs NJ against J²R'^{8/5} available
  differences. No counting contradiction exists at the level of differences.

**(ii) No gap principle at cell scale.**
- Up to four points impose nothing.
- Five points within J cells have a determinant in O_K/16 of norm O(J³). It is forced to vanish only for J < c.
- With six points, the Plücker relation Σ(−1)^j D_j w_j = 0, with all D_j in a finite set, makes the sixth point
  one of O(J^{c}) affine combinations of the other five. The sphere condition then holds identically (Cayley–Menger),
  so there is no contradiction.
- Fixing σ₂-patches does not help either: pigeonhole needs J > J^{9/5}.

**(iii) Fourier/spectral variants stop at 2.**
- *Bilinear (Cauchy–Schwarz over factorisations W = W₁W₂).* Two tubes meet only if they are Clifford-parallel.
  This gives N² ≲ 2^τ · (return count at length 2τ₂) · M(τ₁).
  - With volume laws inserted it gives 8/3, but only under an unproved fourth moment at the unit scale.
  - With the unconditional sup bounds it gives nothing.
- *Quotients near the identity.* Duke-type equidistribution does apply, because the axis cap is macroscopic.
  But the pair count still gives N ≲ ε2^τ, i.e. rate 2.

**(iv) Where 5/2 sits analytically (heuristic, Hopf fibre).** This explains why each route stops where it does.
- *Fibration.* Use the Dirichlet fibration (§3a) with Bezout qa − pb = 1. Then exactly
  |Hλ + κg|² + |g|² = Hn, with H = |p|² + |q|² and κ = q̄b + p̄a. So
  N ≤ #{(λ', g) : N(λ') + N(g) = Hn, λ' ≡ κg (mod H), g ∈ Box}. This is a Clifford tube, where the divisor bound
  holds, plus one congruence on an isotropic line of the quadric mod H.
- *Sizes.* At X⁴ = 1/ε the box has T₀ ≈ εR'^4 points and N(H) ≈ 1/ε. The expected count is T₀/N(H) = E, the
  volume.
- *The dictionary.* An error T₀/N(H)^θ gives rate α = 2 + θ:
  - θ = 0 (congruence ignored) gives 2;
  - θ = 1/2 gives 5/2;
  - θ = 1 gives 3.
- *What θ = 1/2 is.* The error N(H) ≈ 1/ε is the square root of the number q = N(H)² of residue classes of g. That
  is the Pólya–Vinogradov completion error, and it equals the cell count.
- *Shape of the box.* After rebalancing by units of O_K, the g-box is a balanced box of side q^{3/16} against a
  period q^{1/4}: an incomplete sum of length q^{3/4}.
- **So, heuristically, α > 5/2 ⟺ beating Pólya–Vinogradov-type completion for this incomplete sum.** For signed
  1-dimensional character sums of length q^{3/4}, PV is essentially sharp in the worst case. Our quantity is a count,
  not a signed sum, so this is a warning sign rather than an obstruction.
- *Two complementary pictures.* The geometric 5/2 (determinant boxes = cells) and the analytic 5/2 (PV completion)
  coincide.

**(v) Rigid 1-dimensional analogues.**
- For lattice points on circles in short arcs, unique factorisation beats Jarník (Cilleruelo–Córdoba). For n = p^k
  the points form an AP in angle, which gives three-gap rigidity.
- The non-abelian analogue here is words of the free product near a torus coset. No three-gap-type rigidity is known
  for it, and §3d shows that higher-degree auxiliary forms lose on S³ × S³.

**Verdict.** Not solved. Every available input gives exactly 2 (spectral, quotients, square-root inputs) or exactly
5/2 (determinants, PV-type completion). The missing ingredient is cancellation beyond completion for an incomplete
sum, or an equivalent non-abelian gap principle. Both are open problems in their own right.
- *Toy model.* The cleanest place to try is Z⁴: four squares of a prime power near a great circle of S³. There the
  same barrier sits at tube width n^{1/10}: n^{2/5} cells against volume n^{1/5}. Any power saving in that model
  would already be new.

## 6. Third pass: the same count in the sup-norm literature (2026-10-01)

Our tube count is a **Hecke-return count near a torus at a single norm**:
#{γ ∈ O, nrd γ = π^τ, γ within ε of a torus coset}. The automorphic-forms literature needs exactly this count for
sup-norm and geodesic-restriction bounds.

**Single norm: only elementary bounds are known.**
- Marshall (arXiv 1204.0781, Duke 2016), Lemma 3.2: M(ℓ, n, κ) ≪ (κ² + κ^{1/2}) n^{1+ε} + n^ε, uniformly in the
  geodesic, by geometry of numbers. In our normalisation this is ≈ 2^τ ε^{1/2}, weaker than our 5/2 analysis.
- Blomer–Michel (definite case) and Khayutin–Nelson–Steiner (arXiv 2207.12351), Type I/II: also geometry of numbers
  for a fixed lattice R(ℓ)⁰.

**Every improvement averages over the norm.**
- Marshall's Prop. 5.2: Σ_{M/2<m<M} M(ℓ, m, δ)/√m ≪ δ² M^{3/2} (volume law on average) once M ≥ δ^{−2−ε}.
- The proof expands f = (smoothed tube) spectrally and uses Σ_m g(m/M) λ_i(m) ≪ M^{−A} (his Lemma 5.3). That comes
  from the analytic continuation of L(s, ψ_i): shift the contour to Re s = −A.
- KNS bound fourth moments through theta kernels together with a second-moment count over the determinant.

**Why this is unavailable to us.**
- In Clifford+T every norm is a power of the single prime P = (√2).
- The analogue of Σ_m g(m/M)λ(m) is Σ_k g(2^k/M) λ(P^k). Its generating function is the local factor
  L_P(ψ, X) = 1/(1 − a_ψ X + 2X²), whose poles lie at |X| = 2^{−1/2} under Ramanujan. So there is no contour to
  shift: the sum has size 2^{K/2}, with no cancellation.
- This is the square-root barrier again (rate 2), now seen as "one prime".
- Equivalently: the Gibbs sum Σ_τ n(τ) 2^{−λτ} of our rate argument is spectrally Σ_π L_P(π, 2^{−λ}) w_π. For
  λ < 1/2 the evaluation point lies beyond the radius of convergence of every non-trivial L_P.

**Verdict.**
- Pushing to α > 5/2, and a fortiori to 3, needs a *single-norm* bound for Hecke returns near a generic torus
  beyond geometry of numbers.
- That is not known even in the settings (Marshall, Blomer–Michel, KNS) where the community needed it; they avoided
  it by averaging over the norm.
- *Side remark (corrected in §9).* A gate set with generators at many primes, with cost = log norm, does **not**
  escape. Its words of cost ≤ τ number ~4^τ, so the rate-relevant range shrinks to τ ≤ 1.5L, while averaging
  over norms only reaches 2^τ ≥ ε^{−2}, i.e. τ ≥ 2L. For Clifford+T itself, averaging over all norms of size 2^τ
  would cover the relevant range [2L, 3L]. The obstruction is that only the single norm (√2)^τ occurs.

## 7. Fourth pass: no loophole on the process side; the two-round rate is exactly the occupancy exponent

This pass looked for slack in the information-theoretic reduction, so that some averaging or spread condition might
replace the worst-case tube bound. There is none. The following is rigorous.

**Proposition (covering by translates; sharpens Prop. S, and settles "lies between" in Remark what52).**
*Hypothesis.* For some ε, with the tuned Q = Q_ε, some τ and some frame G, let U ⊂ Z_Q be a set such that each
u ∈ U has a word W_u with t(W_u) ≤ τ and d(W_u, R_z(Δu)G) ≤ ε/2.
*Conclusion.* For r = 2 there is a deterministic zero-error CW process with T_pre ≤ mτ and
S ≤ m(log₂(Q/|U|) + log₂ ln Q + 2). No condition on the gaps of U is needed.

*Proof.*
1. **Shifts.** Choose a public set of shifts K ⊂ Z_Q with ∪_{k∈K}(U − k) = Z_Q and |K| ≤ ⌈(Q/|U|) ln Q⌉ + 1.
   It exists by the random choice plus union bound (P[a uncovered] ≤ e^{−|K||U|/Q} < 1/Q).
2. **Round 1.** For share a, pick k ∈ K with a + k ∈ U. Emit W_{a+k} ≈ R_z(Δ(a+k))G, and store k.
3. **Round 2.** Emit V within ε/2 of G⁻¹R_z(Δ(a₂ − k)). This costs nothing in T_pre.
4. **Correctness.** W V is within ε of R_z(Δ(a + a₂)), by unitary invariance and the triangle inequality. ∎

**Corollary (exact two-round rate).**
- Let 𝒰_τ(G) be the set of grid cells reachable along the tube of G at cost ≤ τ.
- The Gibbs lower bound (Cor. profile, cell version) gives α ≥ α* whenever max_G |𝒰_τ(G)| ≤ 2^{τ/α* + o(L)} for all
  τ ≤ α*L.
- The Proposition gives α ≤ τ/(log₂|𝒰_τ(G)| − log₂ L)·(1 + o(1)) for every τ and G.
- Hence the r = 2 exchange rate equals sup{α : max_G |𝒰_τ(G)| ≤ 2^{τ/α + o(L)} for all τ < αL}.
- **So α = 3 for r = 2 is equivalent to the cell form of SA.** Any proof of 3 must prove SA; no modelling or
  averaging loophole exists. The frame is free because the final round is free, and the translates reuse a single
  tube.
- For r ≥ 3, rounds t ≥ 2 need SA for the arithmetic sources φ = P⁻¹|0⟩ (Hecke-ball version, i.e. Conjecture H in
  Hopf form).

**Status of the attempt to prove SA (cumulative).** No route found.
- Every input available is listed in §§3, 5 and 6.
- In the language of §6, the needed statement is a single-norm Hecke-return bound near a generic torus beyond
  geometry of numbers. It is not known in any setting.

Further checks made in this pass, each a dead end:
- *Subspace theorem.* The tube points are Dirichlet-level approximations: the product of the linear forms is
  ≈ R'^{8−2a₀} ≫ 1. The theorem needs a₀ > 4.
- *Cilleruelo–Córdoba products.* These need Liouville in a single complex embedding. Z[ω] has two complex places, and
  σ₂ is free. Pigeonholing σ₂ forces equality only when εη2^τ < 1/4, i.e. N > 16ε²2^{2τ}.
- *Expander mixing on an ε-net.* The error is 2^{τ/2}ε^{−1/2}, worse than Ramanujan.
- *Hecke-stable profiles.* The target 2^{τ/3} grows by 2^{1/3} per level against the tree's factor 2, so no
  inductive argument over levels can produce it.

## 8. Fifth pass (2026-10-01)

**(a) Reduction to arithmetic frames (rigorous).**
- Every great circle lies within ε of an arithmetic one, G₁, G₂ ∈ Γ, of height about 2L, because cosets meet Λ_h
  once h ≥ 2L.
- Hence Q2 holds for all frames ⟺ Q2 holds for arithmetic frames, uniformly in their height.
- *Equivalent tree form.* For V in the Clifford tube (off-diagonal entry small), the count of V with
  d(y, Vx) ≤ τ, for vertices x, y of the tree, should be ≤ 2^{o}(1 + 2^τ ε²). That is, Clifford-tube elements must be
  equidistributed in the 2-adic tree. Theorem clifford pays exactly the factor 2^{d(x,v₀) + d(y,v₀)} for ignoring
  this.
- Splitting off prefixes is circular: Σ over prefixes P of length h of n_τ(P⁻¹T) is the proven Clifford count, and
  each term is the arithmetic-frame count again.

**(b) Beyond GRH.**
- For each (character, harmonic), Deligne/Ramanujan/Lindelöf applied to the theta series bounds the hybrid sums S_χ
  only by roughly the square root of the *complete* count. This is larger than T₀, so it gives at best rate 2.
- So α > 5/2 is beyond what GRH-type inputs give by the standard routes. The unconditional 5/2 is already past them.

**(c) Numerics cannot reach the asymptotic regime.**
- Meet-in-the-middle (W = AB, A ∈ prefixes of length m, B = suffix·Clifford) makes one tube cost O(2^{τ/2}). Close
  pairs Bloch(A^†G₁|0⟩) ≈ Bloch(B G₂⁻¹|0⟩) are found with a spatial hash at resolution 2 arcsin ε.
- But at τ = 5L/2 the ratio volume/cells is ≈ 92·2^{−L/2}, because of the Haar constants (72·2^{τ}ε² against
  π/(4ε)).
- It drops below 1/100 only for L ≳ 26, i.e. τ ≳ 65. At feasible τ ≤ 48 it is ≈ 0.1, so no search can separate
  "volume × poly" from "cells^{1−κ}".
- Not pursued.

**(d) Side observation (corrected after agent_decoupling).**
- Lifting through the Hopf map, the S³ determinant method bounds the points of x² + y² + z² = n² in caps of radius
  r ≤ n^{3/5} by n^{2/5+o(1)}. Monotonicity in r gives this for every r ≤ n^{3/5}.
- The general Bourgain–Rudnick bounds give min(r²n^{−1/2}, r)·n^{o(1)}. Their Lemma 2.2, the "(1 + λ)" bound, gives
  n^{3/5} at r = n^{3/5}, not n^{7/10} as first written here.
- So the lift improves on Bourgain–Rudnick for r ∈ (n^{9/20}, n^{3/5}]. No prior statement of this lift was found.
- Related work: Zhang–Zhu (arXiv:2606.08650) prove ≤ λ^{3/4} lattice points in unit-width tubes around circles on
  λS³ and conjecture λ^ε. Their conjecture would imply the toy-model volume law. The determinant method already gives
  ≈ λ^{2/3} at unit width.

**(e) Delegated directions, reports in this folder:** Burgess/bilinear forms on the hybrid sums S_χ
(`agent_burgess.md`); decoupling and polynomial methods (`agent_decoupling.md`); Linnik basic lemma and entropy
(`agent_linnik.md`).

**(f) Reports received.**
- *Burgess (`agent_burgess.md`): fails.*
  - With the weight dropped, the critical scale is exactly the Burgess endpoint: the dual sum has length q^{1/4}.
    Beating it needs uniform subconvexity C^{1/8−η}, which is beyond Weyl and beyond Milićević's sub-Weyl bound for
    prime-power moduli.
  - The weight a_χ(Hn − N(g)) breaks Burgess's collection step: the summation set is a definite quadric with no
    translation invariance.
  - Type II bounds are trivial: the kernel is a non-negative incidence and all oscillation sits in multiplicative
    coefficients.
  - Dispersion needs an average over the shift or over the modulus, and neither exists.
- *Decoupling (`agent_decoupling.md`): fails. Main structural finding, the window barrier.*
  - The determinant input holds for the union of tube points over all norms in a window |nrd − n| ≲ w², and that union
    has ≍ one point per cell on average (proved in the toy model, sketched for O).
  - So every Archimedean method stops exactly at the cell count. Q1 needs arithmetic that singles out the exact norm.
  - Congruences mod M only reduce the union to cells/M. A power saving needs tree depth ≳ κ log(1/ε).
  - Decoupling caps at the critical width are exactly the cells. Even the conjectured endpoint gives n^{3/5}, against
    n^{2/5} cells.
  - Degree-D determinants: on one-cell boxes the log-norm of the determinant is (6/5)·C(D+2,3)·(D − 1)/2 > 0 for
    D ≥ 2.
- *Linnik / entropy (`agent_linnik.md`): fails.*
  - Basic lemmas, proved:
    - L1: shared prefix c ⟹ separation ≥ c₀2^{−(t−c)/2}, and this is sharp numerically.
    - L2: the joint σ₁/σ₂ pair count gives the volume law at every scale.
    - L3: a set of diameter r contains at most C(1 + r²2^t) words.
  - On tubes these only give N ≤ 2^tε, i.e. rate 2.
  - Tree-class splitting buys K^{1/2} (pairs) or K² (5-point determinant) where K and K³ are needed. Prefix and suffix
    sharing is no better than a longer prefix.
  - **Rigorous no-go.** Adding ε^{−1} points along one tube preserves every ball non-concentration bound with s ≤ 5/2,
    and every L² bound. So no argument whose only inputs are ball or L² non-concentration of μ_t can decide Q1.
    Flattening is already optimal.
  - The Bourgain–Gamburd commutator relation W₂W₁⁻¹W₃ = W₃W₁⁻¹W₂ for points on a torus coset forces exact structure only
    for ε ≤ R'^{−6}, and gives a gap principle only for ε ≲ R'^{−3}.
  - Translates of a rich tube cannot be used either. One-sided translates are parallel cosets, and Ramanujan scales the
    signal and the error equally. Two-sided translates carry too few incidences for any Furstenberg or
    Szemerédi–Trotter contradiction.

**Synthesis (2026-10-01).** Three independent barriers sit exactly at 5/2.
1. *Analytic:* the dual sum has length q^{1/4}, the Burgess endpoint.
2. *Geometric:* Archimedean methods see only a norm window, whose union has one point per cell.
3. *Dynamical:* non-concentration inputs are blind to a single full tube.

Any proof of α > 5/2 needs arithmetic of the exact single norm against a generic torus, of "beyond Burgess"
strength. Nothing in the literature supplies it.


## 9. Sixth pass (2026-10-01, after the user asked for further approaches)

Each of the following fails, at the place indicated.

- **Amplification with auxiliary odd primes.** Each tube point W gives Wγ (γ of norm p^j, near the torus T₂) at
  norm 2^τ p^j. Spectral averaging over p^j (Marshall) gives the composite counts. But lower-bounding by the
  T₂-returns loses exactly their density ε², so one gets only Markov-type bounds: F(C) ≤ mean(F)/ε² = 2^τ.
  Sup bounds need regularity of F at a scale larger than ε, which a rich tube does not have.
- **One global auxiliary polynomial.**
  - *Setup.* Choose, by Siegel's lemma, a degree-D polynomial with K-coefficients that is small on the whole
    σ₁-tube, then use Liouville in both embeddings.
  - *Bookkeeping.* Unknowns ≈ D³/3; conditions ≈ D·m² (jets of order m = μD along the core). This needs
    μ·a₀ ≥ 2(1 + c)/(1 − c) with c = μ²/(1/3 − μ²).
  - *Verdict.* Infeasible unless a₀ ≳ 3.5. So no single hypersurface captures a one-point-per-cell configuration,
    consistent with the window barrier.
- **2-adic (arithmetic-frame) formulation, depth aspect.**
  - Approximate the frame by an arithmetic one (height h ≈ 2L). The tree condition becomes t ≡ s·ū (mod P^h): the
    same congruence structure as §5(iv), now with modulus 2^h.
  - The core object is again the hybrid sum Σ_t ψ(t) a_χ(2^{k'} − N(t)), and the Pólya–Vinogradov size
    √(2^h) = 2^L = cells reappears.
  - Postnikov-type p-adic analyticity does not linearise a_χ along the quadric.
- **van der Corput / Weyl differencing on Λ_τ.** It fails because the free group is non-amenable: |Λ_τ g Δ Λ_τ| ≍ |Λ_τ|.
  Non-amenability is exactly what gives the spectral gap, which is capped at square root.
- **Products of tube words** (left quotients near the Clifford torus, where the volume law is proved; k-fold
  products; near-identity products). They land at T-count 2kτ, where the volumes grow faster than |U|^{2k}, so there
  is never a contradiction. Pair counting is intrinsically rate 2.
- **Fourier decay of local spectral measures.** F_τ = P_τ(A)F₀ on frame space, and the deviation from the main term is
  2^{τ/2} times the Fourier coefficient at "time" τ of the local spectral measure. Bounding that coefficient is the
  count itself.
- **Random-axis families** (target axes random, revealed at the start). Round 1 then needs only
  sup_{x₂} N(x₁, x₂) for almost every x₁. Via local constancy this costs the same high moments (U_k). It also changes
  the paper's model.
- **Standard conjectures.** GRH, Lindelöf, Ramanujan, twisted Linnik (BKS) and the density hypothesis each give at
  most rate 2 for tube *upper* bounds. Ball-covering with asymptotics fails below one expected point per ball. No
  named conjecture short of a Bourgain–Rudnick-type small-cap statement (Conjecture H / SA) implies α > 5/2.
- **Literature.**
  - Burrin–Gröbner (arXiv:2502.17678) treat x² + y² + z² = n², the classical analogue of SA. Individual caps reach
    only n^{−1/4+o(1)}, the square-root scale; almost all caps reach n^{−1/2+δ}.
  - Wang–Zhang (arXiv:2609.18290) is eigenfunction restriction, not lattice points.

**Sharpest formulation of the obstruction.**
- For Clifford+T the rate-relevant range is τ ∈ [2L, 3L].
- Averaged over all norms m with N(m) ≍ 2^τ, the tube count obeys the volume law there: Marshall's Prop. 5.2
  transplanted, since 2^τ ≥ ε^{−2−δ}.
- α = 3 is exactly the statement that the single norm (√2)^τ is not exceptional in this average, at sub-square-root
  scales.

## 10. Seventh pass: new ideas only, no literature (2026-10-01)

Setting: SA in point form. Count the r = t/u ∈ P¹(Q(ζ₈)) with Hermitian height |u|² + |t|² = 2^k exactly (in both
embeddings) and σ₁(r) ∈ D(z₀, ε). Critical scale: ε = 2^{−4k/5}, cells 2^{4k/5}, volume 2^{2k/5}.

**(A) Band decomposition (new and rigorous; the cleanest target found).**
- *Bands obey the volume law.* A latitude band about a Clifford axis is {N(u) ∈ W} (σ₁-window of relative width ε).
  Its count Σ_{m∈W} r(m) r(2^k − m) ≤ 2^{o(k)}·ε4^k is the volume law, proved by the divisor bound. More generally,
  bands about an arithmetic axis φ of height H obey it with loss H².
- *Caps reduce to azimuths.* A cap is the intersection of two bands about different axes. But a second band only
  restricts the lattice points of each circle to arcs, and the divisor bound cannot see arcs. So the cap count is
  exactly the azimuthal question in one band.
- *Fourier expansion.* Detect the azimuth window with Grössencharacters (u/|u|)^ℓ, 8 | ℓ, |ℓ| ≲ 1/ε. These are
  unramified Hecke characters of Q(ζ₈). Then cap count = main + Σ_ℓ ĉ_ℓ e^{−iℓα₀} S_ℓ, where
  S_ℓ = Σ_{m∈W} a_ℓ(m) ā_ℓ(2^k − m): a short shifted convolution of a CM Hilbert form with itself at shift 2^k, of
  length L ≈ ε4^k.
- *Dictionary.*
  - The trivial bound |S_ℓ| ≤ L gives rate 2.
  - A saving L^{1/3} reaches 5/2; a saving L^{1/3+δ} beats it.
  - Square-root cancellation gives 3.
  - Spectral methods give |S_ℓ| ≲ 2^k = √(complete length), i.e. rate 2.

**(B) Excess invariance: no transport proof exists.** Let E = count/volume.
- E is preserved by 2-adic splitting (children at lower level) and by transverse splitting (narrower tubes).
- E is diluted by merging, by left translates, and by right shifts along the torus. A shift by V near T₂ multiplies
  E by ε²η²/(ε + η)².
- At a₀ = 8/5 a saturated tube has E = 2^{L/2}. The window barrier allows exactly this at a₀ = 8/5, and allows more
  at a₀ > 8/5. So no sequence of splittings or translations reaches a scale where a known bound is violated.
- Theorem U's unit-scale route needs θ < 1/4, but the window barrier at the unit scale is 2^{τ/3}.

**(C) Cilleruelo–Córdoba transplanted fails for a structural reason.** On a circle of norm p^e, any m points share
Gaussian factors, of total size ≈ n^{m²/4}, because each prime contributes a 1-dimensional exponent: the tree is a
line. Words of norm P^τ are paths in a 3-regular tree. The minimal total sharing (prefixes and suffixes) is ≈ m²
syllables, i.e. a divisibility of about 2 per pair, against the δR' = R'^{2/5} needed. Non-commutative unique
factorisation forces no sharing.

**(D) Other dead ends.**
- *One global auxiliary polynomial.* Infeasible for a₀ < 3.5 (§9).
- *Exact σ₂-sphere.*
  - Circumcentres: for exact-norm 5-tuples the circumcentre is 0, but that only restates nrd = n.
  - Cayley–Menger and parallelogram relations are identities; the midpoints of exact points are window points.
- *Pancharatnam/Bargmann invariants and cross-ratios.* Their small parts factor as Δ_{24} Δ̄_{13}, i.e. no new
  information.
- *Double caps (σ₁ and σ₂).* Liouville gives separation only below η ≈ 4^{−k}/ε.
- *Sum-product.* The Hopf map is the sum-product structure (|u|² + |t|², ūt), but no small-scale theorem applies.
- *Approximate APs in phase and k-fold products of tube words.* They land at T-count 2kτ, where volumes grow faster.
- *Meet-in-the-middle.* The tube count equals the number of ε-close pairs between two Hecke orbits Λ_{τ/2}x₁ and
  Λ_{τ/2}x₂. Each orbit is spread (≤ 2^{L/6} points per cap), but bounding close pairs needs independence of the two
  orbits at scale ε, which is the problem itself.
- *Arithmetic CM tori and section clusters.* These give O(τ) exact points or fewer than volume, so no candidate
  counterexample exists either.

**Verdict.** No new working mechanism. The obstruction is now located three ways: the window barrier (B),
non-commutativity of factorisation (C), and the short shifted convolution S_ℓ (A).

## 11. Numerics for §10(A) (2026-10-01; details in numerics/RESULTS.md)

Computed S_ℓ exactly for every band and every ℓ ≤ 1/ε, for k ≤ 16 (up to 1.7·10^10 states), at β = 4/5 and β = 2/3.

- **Test A.** The normalised Z = S_ℓ/(8√S_0) has a k-independent distribution:
  - E Z² = 1.00;
  - kurtosis 3.45;
  - 99.99% quantile 4.5 for every k = 13, …, 16.
  Effective θ = 1/2 up to logarithms, far from the θ = 1/3 line.
- **Test B.** Direct cap counts at the critical scale are Poisson: var/mean = 1.00, and the maxima equal the Poisson
  maxima (k = 16: 166 against 168, volume 107). The determinant bound allows about 1.1·10⁴.
- **Exceptions at arithmetic points.** At |+⟩ and at the H-eigenstates there are rich circles (states sharing
  |u − t|²). The excess reaches 5.5× at k = 12, decays to ≤ 1.4× for k ≥ 14, and is divisor-bounded.

Verdict: the numerics support SA and θ = 1/2. They give no proof. *Corrected in §12(b):* θ > 1/3 in mean square over ℓ
already suffices. The whole-sphere pair count is provable, but the in-band one is not, and an in-band asymptotic
would not be blind.

## 12. Eighth pass (2026-10-01): necessity, a variance route, Hecke-ball numerics

**(a) Any proof of α = 3 must prove SA.**
- "α = 3 for every r" includes r = 2, and by §7 the two-round rate is 3 ⟺ SA. So there is no route to 3 that avoids
  SA, an S-arithmetic Bourgain–Rudnick small-cap theorem for a single Hecke orbit.
- Its rational toy (Zhang–Zhu's tube conjecture on λS³, §8(d)) is open. The best bound there is λ^{2/3}, from the
  determinant method, which is the same barrier as our 5/2.

**(b) Variance route (new reduction, rigorous).** Work in a z-band R of width ε, cut into cells of side ε.
- *The inequality.* max_c |N_c − vol| ≤ (Σ_c (N_c − vol)²)^{1/2} = (P(R) − B·vol)^{1/2}. Here P(R) is the
  number of pairs of states in R lying in the same cell, B = |R ∩ states| ≈ ε4^k, and vol = εB/2π.
- *Critical balance.* At the critical scale, B·vol = cells² exactly.
- *Reduction.* Q1 (r = 2) ⟸ P(R) = B·vol·(1 + O(2^{−δk})) for every band. This is a power-saving *asymptotic* for a
  two-point statistic. The no-go of §8(f) concerns L² *upper* bounds with constant or 2^{o(k)} slack, so it does
  not apply here.
- *The global version holds (rigorous).* The smoothed whole-sphere pair count is a sum of squares,
  Σ_φ k̂(φ)|λ_φ(P^{2k})φ(z)|². Ramanujan gives it with relative error O(k²/vol). But the whole sphere gives only
  sup ≤ vol + √N, the Ramanujan bound.
- *Localising to the band destroys positivity.* The band statistic is v*Gv with v_φ = λ_φ φ(z) and an
  off-diagonal Gram matrix G (‖G‖ = |R|). Ramanujan with Bessel's inequality gives only P(R) − B·vol ≲ k²N, i.e.
  Ramanujan again.
- *What would be needed.* The L² mass of the fluctuation F = Σ_{φ≠1} λ_φ(P^{2k})φ(z)φ (the smoothed state measure
  minus its mean) must be equidistributed in z-bands of measure ≲ vol⁻¹, with a power saving. Spectrally this is
  cancellation in Σ λ_φ λ_φ' φ(z) φ'(z) φ''(z) ⟨φφ'φ''⟩ over ~ε⁻⁵ triples, i.e. effective QUE for a combination of
  ε⁻² Hecke forms at a shrinking scale. That is not available: Brooks–Lindenstrauss QUE on S² is ineffective, and
  Watson–Ichino plus subconvexity gives rates only for fixed test functions.
- *Arithmetic form, possibly worth more thought.*
  - P(R) − B·vol = Σ_D (angular Weyl sums of the representations by the binary form Q_D on D^⊥).
  - The sum runs over differences D ∈ O_K³ with |σ₁D| ≤ ε2^k. The discriminants |D|² have arbitrary norms, so
    unlike the one-point count this *is* a family.
  - The band cuts each circle C_D to two short arcs. So one needs Hecke-angle equidistribution in short arcs for
    K(√−|D|²), on average over D. Not pursued.

**(c) Hecke-ball numerics (Test C, rounds ≥ 3; numerics/orbit.c, RESULTS.md).**
- *Method.* Take the orbit Λ_t·φ of a source φ = V|0⟩, with V of T-count h. All 72·2^t − 48 words are enumerated
  in Matsumoto–Amano form and binned at the critical scale s = 2^{−2t/5}.
- *Results at t = 28 (λ = 278 per cell; t = 30 agrees, see RESULTS.md).*
  - h ≥ 8: var/mean 0.99–1.02. The maximum is 375–382, against a Poisson maximum of ≈ 380.
  - h = 4: var/mean 1.42, from short stabiliser clusters.
  - h = 0: the T-multiplicity 8.
- So Conjecture H in Hopf form, with arithmetic sources, also looks Poisson.
- *Correction to §8(c).* Measured counts at T-count 28 are ≈ 1% of the one-per-fibre-cell level. So numerics at
  T-count 28–32 *do* separate volume from cells; §8(c)'s ratio formula was too pessimistic. They still cannot
  separate "volume × poly" from "cells^{1−κ}" asymptotically, of course.

## 13. Ninth pass (2026-10-02): the difference-vector family, moment criteria for r = 2

**(a) Moment criterion for two rounds (new, rigorous).** Caps on S² form a 2-parameter family, whereas tubes in S³
form a 4-parameter one. Hence, for r = 2:
- *Hypothesis.* For all scales, Σ_caps N_c^{2k} ≤ 2^{o(t)}·#caps·(vol + vol^{2k}). That is, 2k-tuples of states in
  a common ε-cap are at most 2^{o} times Poisson.
- *Bound.* Then max N ≤ (#caps·vol^{2k})^{1/(2k)} = ε^{−1/k}·vol, and the Gibbs bound gives rate ≥ **3 − 1/k**
  (worst case τ = αL).
- *Cases.*
  - k = 1 is provable (pair correlation, §12(b)) and gives rate 2.
  - k = 2 gives exactly 5/2.
  - k = 3, a 6-point upper bound, would give **8/3**.
- *Comparison.* This needs only *upper* bounds with 2^{o} slack, unlike the variance route of §12(b), which needs a
  power-saving asymptotic. For general tubes the same argument gives 3 − 2/k (conditional3).
- *Why it is still open.*
  - 6-tuples in a cap correspond to tuples of Clifford-tube elements with a common first half. Hence
    Σ_caps N⁶ ≈ Σ_ψ N(cap ψ)⁵, which is circular.
  - The Plücker parametrisation ψ_j = (f_{1j}ψ₂ − f_{2j}ψ₁)/f₁₂ only reorganises the count.
  - S² determinants do not vanish at the critical scale: the 6×6 quadratic-monomial determinant needs ε < 2^{−2k}.

**(b) The difference-vector family.** Write Ŵ = P + P′ and D = P′ − P. Then Ŵ ⊥ D, |Ŵ|² + |D|² = 4N and Ŵ ≡ D (2).
- *Rewriting.* P(R) = Σ_{Ŵ ∈ band, |Ŵ|² ∈ 4N − window} r_{S+S}(Ŵ). The Ŵ run over the spheres of all norms
  M = 4N − |D|², about 2^{2.4k} of them. This is a genuine average over norms.
- *Weights.* r_{S+S}(Ŵ) = r_{L_Ŵ}(4N − M), where L_Ŵ = Ŵ^⊥ ∩ O_K³ is a binary lattice of discriminant ~M.
- *Trivial genus character (provable).* Its contribution is a lattice-point count in a 6-dimensional region, thin
  in σ₁ but long in σ₂. No dual vector fits, since ε²4^k > 1.
- *Class-group characters χ of K(√−M) (about 4^k per M).* Their contribution is Σ_{M,χ} h_M⁻¹ a_χ(4N − M)·W_χ(M; band),
  with twisted band Weyl sums W_χ.
  - Trivially this is 2^{3.6k}, against the main term 2^{1.6k}. Square-root heuristics give 2^{0.8k}.
  - Cauchy–Schwarz over (M, χ) fails: Σ_χ|W_χ|² = h_M·#{same-class pairs in the band}, and the diagonal h_M·εr₃(M)
    dominates, since a band holds 2^{1.2k} ≪ h_M ≈ 4^k points.
  - Pointwise Lindelöf for L(½, π_f × χ) is also not enough. One needs cancellation over the family.
- *Verdict.* The family exists, but its cancellation is not reachable by large-sieve or Cauchy–Schwarz arguments.

**(c) The localisation barrier (summary of §12(b) and this pass).**
- Every second-moment route faces the same step: localise the provable global variance (≈ N) to a set of
  measure ε. Routes tried: in-band pairs, residue classes mod P^τ via midpoints of Clifford-tube words,
  Hecke-orbit unions of caps, recursive child splitting.
- Each time Cauchy–Schwarz loses exactly the missing factor (e.g. √(|G|·q·vol) against ε|A|).
- Bounds per eigenform are fine: latitude circles are curved, so ∫_band|φ|² ≲ ε·l^{1/3}. The obstruction is
  entirely in the cross terms between eigenforms.
- Bourgain–Gamburd non-concentration near cosets of SO(2) is far too weak (ε^κ against the needed ε^{3/2+δ}).

## 14. Two criteria not stated above, and the remaining dead ends of passes 8–9 (2026-10-02)

**(a) Centred fourth moment (r = 2).**
- *Statement.* Suppose Σ_caps (N_c − vol)⁴ ≤ 2^{o(τ)}·#caps·(vol + vol²) at every scale, with ε = 2^{−L} and
  vol = 2^{τ−2L}. Then max |N_c − vol| ≤ (#caps·vol²)^{1/4} = 2^{τ/2 − L/2 + o}. This is ≤ 2^{τ/α} for all
  τ ≤ αL exactly when α ≤ 3, so the centred fourth moment alone gives **α = 3 for r = 2**. (For vol < 1 the bound
  2^{τ/4} suffices.)
- *Why it is hard.* Unlike §13(a), this needs main terms to cancel. It needs the 2-, 3- and 4-point correlations
  at scale ε with relative errors 1, vol⁻¹ and vol⁻² respectively.
- *Equivalent forms.* It is L²-equidistribution, with Poisson variance, of the close-pair (midpoint) measure,
  i.e. the midpoints of Clifford-tube words of T-count 2τ. Spectrally it is ‖F‖₄ ≲ 2^{o}‖F‖₂, i.e.
  Σ_ψ |Σ_{φ,φ′} a_φ a_φ′ ⟨φφ′ψ⟩|² ≲ 2^{o}‖a‖⁴ with a_φ = k̂(φ)λ_φ(P^τ)φ(z).

**(b) Spectral dictionary of the obstruction.**
- *The spectral sum.* A cap count equals Σ_f c_f λ_f(P^{2k}). It runs over the Hecke forms of weight ≤ 1/ε on D,
  i.e. Hilbert forms of weight (d + 2, 2) and level 1 over Q(√2). The family has size F ≈ ε⁻².
- *Twist length.* At the rate-α scale the twist λ_f(P^{2k}) has length N(P^{2k}) = 4^k = F^{α/2}. So:
  - α = 2 corresponds to twist length F, the natural reach of Petersson and trace-formula methods;
  - 5/2 corresponds to F^{5/4};
  - 3 corresponds to F^{3/2}.
- *Petersson over K.* The σ₂-weight is 2, so the Bessel factor J₁ at σ₂ does not decay for small argument and the
  off-diagonal Kloosterman terms are not suppressed. The family varies in the σ₁-weight only.
- *Satake angles.* S_ℓ = 2^k Σ_f c_f U_{2k}(cos θ_f). Vertical Sato–Tate at the single prime P (Serre) is
  qualitative only, with discrepancy ≍ 1/log.

**(c) Further dead ends (each checked, none recorded earlier).**
- *Rarity plus propagation.* L² makes rich caps rare (≤ vol²ε^{−2δ} of them). A rich cap C at level τ gives rich
  caps gC (g ∈ Λ_s) at level τ + s. But the L² allowance at level τ + s is 4^s·vol², more than the 2^s produced.
- *Determinant sup plus L².* Interpolation gives the 2k-th moments only for ε ≥ N^{−1/3}, above the critical scale.
- *Double caps (σ₁-radius ε, σ₂-radius η) with a two-embedding determinant.* Four points are coplanar once
  εη < 2^{−3k/2}. Pigeonholing σ₂ then costs ε²8^k.
- *Additive energy of the state set.* The exact energy is minimal (E ≤ 2^{o}|S|²; each circle C_Ŵ carries 2^{o}
  points). But close quadruples are not additive quadruples: the σ₁-approximate energy at scale ε2^k exceeds the
  close-quadruple count by ε⁻³.
- *Circle method for the pair count* (12 integer variables, 4 quadratic equations, thin σ₁ conditions). Far outside
  Birch's range.
- *Genus-2 Siegel theta series.*
  - Orthogonal pairs (Ŵ, D) are representations of diag(4N − d, d) by I₃.
  - With harmonics, the Maass relations of the Saito–Kurokawa lift turn the pair statistic into half-integral weight
    coefficients c(4d(4N − d)), summed along a quadratic sequence over a thin window.
  - Not pursued: it is the same family-cancellation problem as §13(b).
- *Upper-bound sieves* (large sieve modulo auxiliary primes). They give only the dimension bound for points on a
  surface and are blind to the global size.
- *Incidences.* By the determinant method the points of a cell lie on one circle. But each circle carries ≤ 2^{o}
  states, and different cells give different circles, so there are no rich incidences.
- *Lift to composite norms* (w ↦ wg₀ with nrd g₀ = m). This gives N_n(C) ≤ N_{nm}(Cg₀). But averaged-norm results
  control all norms ≍ 2^τm, and {nm} has density 2^{−τ} among them.
- *Additive perturbation* (w ↦ w + δ, δ small and near the plane). This reproduces the window count R′².
- *Random target axes.*
  - Lower bounds for random targets would imply worst-case ones.
  - But the frame is chosen by the process, so the supremum over frames remains.
  - For almost every source one still needs 2k-th moments over a 4-parameter family (§9).

**Status after §14.**
- *Unconditional:* still 5/2.
- *Any proof of 3 must prove SA* (§12(a)).
- *Cleanest open targets, in order of apparent tractability. Each would give, for r = 2, at least:*
  - a 6-point correlation upper bound gives 8/3 (§13(a));
  - a power-saving in-band pair-count asymptotic gives > 5/2 (§12(b));
  - θ > 1/3 in mean square for S_ℓ gives > 5/2 (§10(A), §12(b));
  - a centred fourth moment gives 3 (§14(a)).
- *Rounds ≥ 3* need the same statements for Hecke-ball sources V|0⟩, uniformly in the height of V.
- All of these are beyond-square-root statements for a single Hecke orbit. No known technique reaches them.

## 15. Tenth pass (2026-10-02, new session): three precise statements, no new mechanism

Per-round log in `PROGRESS.md`. Nothing here moves the unconditional rate.

**(a) The norm-averaged volume law, stated precisely (routine; Marshall Prop. 5.2 + Jacquet–Langlands).**
X_m := {w ∈ O : nrd w = m}, m ∈ O_K⁺ mod (3+2√2)^Z; N_m(f) := Σ_{w∈X_m} f(σ₁w/√σ₁m) for f ∈ C^∞(SO(3)) a
smoothed ε-tube indicator with Sobolev norms ≲ ε^{−C}. For g smooth on [1/2,2]² and M = M₁M₂,
  Σ_m g(σ₁m/M₁, σ₂m/M₂) N_m(f) = (∫f)·Σ_m g(⋯) r_O(m) + O_A(ε^{−C}M^{−A})   whenever M ≥ ε^{−2−δ}.
Proof: N_m(f) = Σ_ψ ⟨f,ψ⟩ψ(x₀)λ_ψ(m) over Hecke eigenforms of weight (a,0), a ≲ 1/ε; non-trivial ψ are cuspidal
on GL₂/K by JL; detect the box with the unramified Grössencharacters λ_ν of K; L(s, π_ψ ⊗ λ_ν) is entire, shift to
Re s = −A, conductor ≍ a²(1+|ν|)² ≲ ε^{−2}.
- P^τ has weight ≍ 2^{−τ} in the average. Positivity loses the full 2^τ; the second moment over norms gives
  (2^τ·vol)^{1/2} = 2^τε, worse than Ramanujan.
- The sub-family {P^τl : (l,P) = 1, N(l) ≍ Λ} has negligible spectral error for Λ ≥ ε^{−2−δ}, by multiplicativity
  λ_ψ(P^τl) = λ_ψ(P^τ)λ_ψ(l). But X_{P^τl} = X_{P^τ}X_l/48, so the statement is Σ_{v∈X_l} N_{P^τ}(cap·v) = volume:
  Hecke-orbit averaging of caps, i.e. agent_linnik §3(d) ("amplification is exactly neutral"). No content.
So "α = 3 ⟺ P^τ is not exceptional in the average" has no operational content.

**(b) Odd moments.** max_c N_c ≤ (Σ_c N_c^k)^{1/k} needs no parity. A k-point correlation upper bound
Σ_c N_c^k ≤ 2^{o}·#caps·(vol + vol^k) at all scales gives rate 3 − 2/k (r = 2):
2 ⇒ 2, 3 ⇒ 7/3, 4 ⇒ 5/2, **5 ⇒ 13/5 > 5/2**, 6 ⇒ 8/3. The 4-point bound itself is not known: the determinant
sup bound interpolated with L² misses it by ε^{−1}/vol = ε^{−1/2} at the critical scale (§14(c)). Test D
(numerics/RESULTS.md) checks m ≤ 8 at k ≤ 16: Poisson to three decimals.

**(c) Three-point correlation, four equivalent forms (all exact).** With V ∈ T_ε(T_z) ∩ Λ_{2τ},
n(P) := #{V : prefix_τ(V) = P}, T = Σ_P n(P) ≈ ε²4^τ (Theorem A), vol = ε²2^τ:
  Σ_P n(P)² = #{(V,V′) ∈ tube² sharing a prefix of length ≥ τ}
            = #{(u,t),(u′,t′) ∈ tube² : ut′ ≡ tu′ (mod √2^τ) in Z[ζ₈]}
            = Σ_{D = A⁻¹B ∈ T_{2ε}(T_z), A,B ∈ Λ_τ} #{P ∈ Λ_τ : PA, PB ∈ tube}
            = Σ_j ⟨A′_{2j}1_tube, 1_tube⟩  (A′_{2j}: shell adjacency at tree distance 2j).
By L1 the j-th term vanishes for j < 2log₂(1/ε) − O(1) (= 4τ/5 at the critical scale); for j ∈ [4τ/5, τ] it is
Σ_{ε-close (A,A′) ∈ Λ_j²} #{P ∈ Λ_{2τ−j} : PA, PA′ ∈ tube}: close pairs at depth j (known) times cap counts at depth
2τ − j around arithmetic points of depth j (Theorem A loses 2^j). The j = τ term is Σ_A N_τ(cap_A)² itself. This
confirms §13(a): every rewriting returns to cap counts at depth τ around depth-τ arithmetic points.
- *Bloch form.* ν = Bloch(u,t) ∈ Z[√2]³, |ν|² = 4^τ in both embeddings. Prefix class P ⟺ ν ∧ ν_P ≡ 0 (mod √2^τ).
  P-adically |ν|² → 0, so [ν] lies near the isotropic conic C₀ ⊂ P²(K_P), C₀(K_P) ≅ P¹(K_P) = ∂(tree); prefix
  classes are depth-τ balls on the boundary. **SA for r = 2 is the joint distribution of the integer points of
  x²+y²+z² = 4^τ over Z[√2] in (an ε-cap at σ₁) × (a depth-τ ball on ∂tree), at the single norm 4^τ.** Joining-type
  results (Aka–Einsiedler–Shapira style) are qualitative and need a second prime.
- *Operator form.* M(P,A) := 1[PA ∈ tube] is a 2^τ×2^τ 0/1 matrix with T ones (one per tube element, midpoint
  split). Σ_P n(P)² = ‖M1‖². ‖M‖_op ≤ 2^{o}vol would give the three-point bound for all weights. ‖M‖⁴ ≤ #4-cycles
  = the j ≥ 1 part of the multiplicative energy of the tube. Peter–Weyl splits M into ≈ ε^{−3} rank-one pieces; the
  triangle inequality gives ‖M‖ ≤ vol·ε^{−3/2}; the needed cancellation between pieces is the problem.

**(d) In-band pair count is a sum of squares (exact).** With a positive-definite azimuth window,
P(R) − B·vol = Σ_{ℓ≠0} ĉ_ℓ|S_ℓ|² ≥ 0, ĉ_ℓ ≈ ε. So the power-saving asymptotic of §12(b) is exactly
Σ_{0<|ℓ|≤1/ε}|S_ℓ|² ≤ B²2^{−δk}, i.e. RMS_ℓ|S_ℓ| ≤ L^{2/3−δ′} at the critical scale: §12(b)'s mean-square θ > 1/3,
now with the sign information that the deviation is non-negative. The large sieve over ℓ returns the trivial bound.

**Status after §15.** Unchanged: 5/2 unconditional. The cleanest open targets are those of §14, with 5-point
correlation (⇒ 13/5) added below 6-point.

## 16. Eleventh pass (2026-10-02, cloud session; per-round log in PROGRESS.md §4–§7)

**(a) The sign in §15(d) gives nothing beyond Ramanujan (rigorous accounting; details in `scratch/round4_sign.md`).**
With D(W) := Σ_c(N_c − vol)² = Σ_{ℓ≠0}ĉ_ℓ|S_ℓ|² ≥ 0 and the critical-scale dictionary 1/ε = vol², B = vol³, N = vol⁵:
- The functionals of {S_ℓ} with usable bounds are: per-ℓ |S_ℓ| ≤ √N·2^{o} (the shifted-convolution spectral bound of §10(A), taken
  from there; S² Ramanujan on the ~1/ε harmonics of e^{iℓα}1_W gives only √N ε^{−1/2} = N^{7/10}, worse than the trivial B); the large sieve Σ_{|ℓ|≤Λ}|S_ℓ|² ≤ (Λ+B)B; the isotropic
  whole-sphere sum of squares Σ_W D(W) ≤ N k^C; and exact long-range Parseval Σ_{ℓ mod Λ}|S_ℓ|² = ΛB·2^{o} for Λ ≥ B.
  The gap is Σ_{|ℓ|≤Λ}|S_ℓ|² for Λ ∈ [vol, vol²]: per-ℓ Ramanujan reaches exactly B² at Λ = vol, the target is B^{2−δ} at Λ = vol².
- Positivity transfer D(W₀) ≤ Σ_{W∈F}D(W) costs |F| and gains the family's power saving. Window sums (|F| = 1/ε, gain vol):
  deficit vol^{1+δ}, i.e. the Ramanujan sup bound √N. Axis averages and Hecke translates of the axis (orbit-in-band counts need
  2^s ≥ ε^{−2−δ}) have deficit vol³. The Gram matrix over windows is PSD with total ≤ N k^C but signed off-diagonal.
- A rich cell N_c = X vol forces RMS_{|ℓ|≤2vol²/X}|S_ℓ| ≥ X vol/2; against per-ℓ Ramanujan this is X ≤ vol^{3/2}, against the
  large sieve X ≤ vol. Nothing new.
- *Withdrawn (round 6 correction).* The round-4 claim that azimuth cells of width η ≥ vol^{−1+δ} give a power-saving in-band
  Poisson asymptotic was wrong: it majorised the anisotropic (ε-band × η-cell) pair count by the isotropic η-pair count, whose
  main term is larger by η/ε. Correct accounting (Schur bound on the anisotropic kernel) gives relative error k²/vol only for
  the sum over all bands, i.e. nothing for a single band. The gap remains the whole range 0 < |ℓ| ≤ 1/ε.

**(b) P-adic Fourier form and the depth-split determinant bound (round 5; `scratch/round5_padic.md`).**
- Class P ⟺ ν ∈ L_P := O_Kν_P + P^τO_K³ (index 2^{2τ}). Exactly n(P) − vol = 2^{−2τ}Σ_{η∈L_P^⊥∖0}Ŵ(η) with
  Ŵ(η) = Σ_{ν∈X∩cap}e_{P^τ}(⟨η,ν⟩). Per-η square root with absolute values gives vol³; Parseval over L_P^⊥ counts cap pairs
  with difference in L_P and its diagonal N_cap = 2^{6τ/5} exceeds the target vol². Ŵ(η) is a theta coefficient with
  characteristic η/P^τ at index = level (depth aspect, §14(b)). Grouping η by valuation is the Hecke recursion.
- *Withdrawn (round 6 correction).* Round 5 claimed that for prefix length p ≥ t/3 each prefix class contributes ≤ 2^{o} points
  to an ε-cell, assuming the determinant bound for the suffix orbit improves below its critical scale. It does not: by Theorem
  rate52 / Remark what52 the tube bound is "≤ 2^{o} words per δ-cube along the core" with δ the critical cube for the T-count,
  i.e. 2^{2t′/5} for T-count t′ at every radius ≤ 2^{−2t′/5}. A class at depth p therefore contributes ≤ 2^{2(t−p)/5+o}, and
  summing over 2^{p} classes gives 2^{2t/5+3p/5}: a loss of 2^{3p/5}, confirming §3(e). The hitting-set reformulation H_p is
  void. What survives: the structure of a hypothetical rich cap at the critical scale is "one frame per ε-arc of the Hopf fibre"
  (the determinant method constrains the fibre coordinate, not sub-cells of the cap).
- *Exact recursion of deviation fields.* F_j(z) := N_{cap(z,ε)}(j) − main; T₂F_k = F_{k+1} + F_k + 4F_{k−1} with T₂ = Σ_{Λ₂}g^*.
  On L²₀ Ramanujan gives spec T₂ ⊂ [−3,5], and each component grows exactly like 2^k: the Poisson growth saturates Ramanujan.
  Backward propagation of a spike against the determinant bound at lower levels gains nothing (factor 1.137^s against it).
- Literature (search only): no individual-cap bound beyond the determinant method for x²+y²+z² = n (Bourgain–Rudnick
  F₃(R,λ) ≪ R^ε(1+λ²)); Humphries–Radziwiłł is variance; Burrin–Gröbner averages over heights.

**(c) Black-box insufficiency model (round 6; `scratch/round6_model.md`).** Let (R) be Ramanujan Weyl sums at every frequency,
(D) the 5/2 tube bound ("≤ 2^{o} words per critical δ-cube along the core", flat below δ), (C) the divisor bound on arithmetic
circles, (L) Liouville/L1 separation, (V) the band volume law about arithmetic axes, (H1) the exact Hecke recursion of deviation
fields, (H3) joint σ₁/σ₂ equidistribution at the Ramanujan rate. A configuration consisting of a Poisson-like background plus
M = 1/ε = vol² frames over one critical ε-cap, one per fibre ε-arc, satisfies all of (R),(D),(C),(L),(V),(H1),(H3) (checked
item by item) and saturates the trivial in-band pair bound. Hence no combination of these inputs by positivity, Cauchy–Schwarz,
Hölder, interpolation or Hecke recursion gives N_c ≤ vol^{2−δ}. It is killed exactly by the criteria of §12–§15 (5-point
correlation, in-band pair asymptotic, centred fourth moment), and by nothing weaker. The only input not absorbed by the model
is the exact arithmetic of Z[ζ₈] (band counts and circle multiplicities as identities rather than divisor upper bounds).

**(d) Exact arithmetic form of a critical rich cap (round 7; `scratch/round7_arith.md`).** With Θ(m) := {arg u : |u|² = m}
(a set of r(m) ≤ 2^{o} signed sums of prime Hecke angles of Z[ζ₈]/Z[√2]), N_c ≤ 2^{o}#{m ∈ W : (Θ(n−m) − Θ(m)) ∩ (α₀ ± ε) ≠ ∅}.
SA at the critical scale ⟺ the multiset {Σ_{𝔮|n−m}±θ_𝔮 − Σ_{𝔭|m}±θ_𝔭 : m ∈ W} of size ≈ B is 2^{o}-uniform at scale B^{−2/3};
its Fourier coefficients are the S_ℓ of §10(A). The exact identities (band count, Lagrange identity for quotients, Bloch vector)
all express the cap count as an angle-selected sub-sum of an exact divisor sum and give no second identity for the selected part.

**Status after §16.** Unchanged: 5/2 unconditional. New: the black-box insufficiency model of (c) shows that every input used
in §3–§16 is consistent with a vol² spike, so a proof needs a genuinely new input; the cleanest open targets remain those of §14
(5-point correlation ⇒ 13/5, in-band pair asymptotic, centred fourth moment).
