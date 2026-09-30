# Referee report: `NOTES_moments.md` (high moments over frames), §§0–4 and raw §5 data

Scope: §§0–4 of the note, `moments.cpp`, `momtable.py` and `mom_tau10..20.txt` on the server. I ran light
checks with scratch code in `school-server:/tmp/rev_mom/` (`permom2.cpp`, `clus_exact.cpp`). Those files are
not in the repository, and I did not touch `run_mom.sh`.

| Claim | Verdict |
|---|---|
| §0 Hopf form | VALID (exact) |
| (a) Lemma 1 and the sup bound | VALID |
| (b) Lemma 2 (tiny tubes) | VALID WITH FIXES (minor) |
| (c) Theorem M | VALID WITH FIXES (one wrong side remark) |
| (d) Prop. 3 and its corollary | VALID WITH FIXES (choice of τ2) |
| (e) Theorem S | VALID WITH FIXES (commutativity argument, summability); its "consequences" are OVERSTATED |
| §4 heuristics | OVERSTATED in three places |
| Numerics: code | SOUND |
| Numerics: author's reading | WRONG (mechanism) and UNSUPPORTED (convergence) |

## §0 and (a)

- **Hopf form.** With W′ = G1⁻¹WG2⁻¹ = [[a, −b̄], [b, ā]], the tube condition is |b|² ≤ ε²(1−ε²/4) (eq. btube). The polar
  angle ψ of R(W′)ẑ has sin²(ψ/2) = |b|², and sin(φ/2) = ε√(1−ε²/4) for φ = 4 arcsin(ε/2). So the equality in §0 is
  exact, and E = N_τ ε²(1−ε²/4) exactly.
- **Lemma 1.** The triangle inequality and the product-of-caps argument are correct. Moreover s(φ/2) = sin²(arcsin(ε/2)) =
  ε²/4 *exactly*, so n(τ) ≤ (16 M_2k(2φ)/ε⁴)^{1/2k} holds with no (1+o(1)). Drop that factor in Thm M step 2.

## (b) Lemma 2

The checks requested all go through.

- **E1 with λ = 1/2, b = 0.** The proof of E1 only uses the side lengths of the box:
  - core direction: 2R′ sin(1/4) ≤ R′/2;
  - e1 direction: R′(1−cos ¼) + O(ε²R′) = 0.0311R′ ≈ R′/32;
  - two transverse sides of 4εR′ each, which is generous since the transverse offset is ≤ εR′.

  The σ1-determinant is ≤ 24·Π sides ≤ Cε²R′⁴, and the σ2-determinant is ≤ CR′⁴. Their product is CR′^{8−2a0}, which is
  below 2⁻⁸ once a0 ≥ 4+η and R′ ≥ R0(η). The hypotheses 2a+3b > 8 and a0+a ≥ 2b become a = a0 > 4 and a trivial
  inequality.
- **Stated range.** App. F says "E1–E5 hold for any a0 ∈ (0,2)". The proofs of E1, E5 and F2 use no upper bound on a0.
  Say so explicitly.
- **Volume.** A hyperplane section of a box has vol₃ ≤ 2 × (product of the three largest sides), by projecting to the
  coordinate hyperplane of the largest normal component (as in F3(b)). Here that gives 2·(R′/2)(R′/32)(4εR′)(1+O(ε²)) =
  εR′³/8·(1+o(1)). Correct, and only the order matters.
- **F2 with y = 0.** Admissible: every primitive normal lies in 𝒪* ⊂ 𝒪_K⁴, so h(n) ≥ 1 by (H2). H is K-rational because it
  is the affine span of ≥ 4 points of 𝒪. Then v = 3−a0+o(1) gives 1 + CR′^{2−a0/3} circles, as claimed.
- **Fixes.**
  1. The final bound is 2^{max(0, τ/2−L/3)+o(τ)}. As written, "≤ 2^{τ/2−L/3+o(τ)}" is false for τ < 2L/3.
  2. "For τ ≤ L" should read "for τ ≤ (1−η′)L".
  3. Add the case where the span has dimension ≤ 2 (a circle, so E5 applies).
  4. R0 and C depend on η. Summing over T-counts t ≤ τ costs a factor τ+1, and a0(t) only grows as t decreases.
- **Optional simplification.** The *linear* 4×4 determinant of numerators is ≤ Cε²R′⁴·CR′⁴. So for ε < cR′⁻⁴ the whole
  tube lies in one K-rational 3-space through 0, and no arcs are needed. For ε < cR′⁻⁶ the linear 3-wedge
  (σ1 ≤ CεR′³) forces a K-plane, i.e. an arithmetic coset, which gives O(τ) words.

## (c) Theorem M

- **Case E ≥ 1.** The bound is log n ≤ 2L/k + τ − 2L + o. It is increasing in τ, and at τ* ≤ αL + O(1) (log2 K < L +
  log2(π/2)) the excess is O(1)(1−1/α), which goes into A. So α ≤ 3−2/k is right.
- **Case E < 1.** Here n ≤ 2^{(τ+2L)/2k+o}, and the condition (1−η′)(2k−α) ≥ 2α is right. At α_k it reduces to
  2k²−9k+6 ≥ 0. The roots are (9±√33)/4 ≈ 0.81 and 3.69; the value at k = 4 is 2 > 0; and k = 4 allows η′ ≤ 1/11.
  Lemma 2 covers τ ≤ (1−η′)L for every α ≤ 3, since there τ/2 − L/3 ≤ τ/6.
- **A = o(L).** This holds because sup_{τ≤3L} η(τ) = o(L) whenever η(τ) = o(τ).
- **Does (U_k) suffice?** Yes. For E < 1, (M_k) is used only as M_2k ≤ 2^{o}(E+E^{2k}) ≤ 2^{o+1}E, which is (U_k). For
  E ≥ 1, use (d) with the fix below. The regime actually needed is τ ∈ [(1−η′)L, 2L+O(1)] at radius 2φ(ε), including
  very small E (≈ 2^{−L} at τ ≈ L).
- **Fix.** "k = 4 gives 5/2 (already unconditional)" is wrong as stated. The paper has 17/7. The research notes have
  249/100 (annulus) and every α < 5/2 (uniform52). So k = 4 adds only the endpoint.
- **Add.** The tube form of Conjecture H (n(τ) ≤ 2^{o}(τ+1+2^τε²)) implies (U_k) for every k, because
  M_2k ≤ E·(sup N)^{2k−1}. This confirms that (U_k) is a genuine averaged weakening.

## (d) Proposition 3 and its corollary

- **Factorisation.** Correct. It needs no normal form: any word with t ≤ τ T-gates splits after its first
  max(0, t−τ2) T-gates. The rest of the argument is also correct:
  - left invariance d_proj(UV, C) = d_proj(V, U⁻¹C);
  - the coset (G1, G2) maps to (R(U)⁻¹x1, x2), which is σ×σ-preserving;
  - Minkowski's inequality.
- **c_Λ ≤ 72.** Correct. It is equivalent to 72(2^{τ1}+2^{τ2}) ≥ 120, and it equals 24 when τ1 = 0.
- **Choice of τ2 (fix).** The step ratio N_{t+1}/N_t is 4, 2.5 and 2.2 for t = 0, 1, 2, and E_0(φ) = 24s(φ) may already
  exceed 1. So an E_{τ2} ∈ [1/2, 1] need not exist.
  - Take τ2 := max{t ≤ τ : E_t(φ) ≤ 1}. Then E_{τ2} ∈ (1/4, 1], and the constant becomes 72·4^{1−1/2k}.
  - If E_0 > 1, use τ2 = 0 and the trivial κ ≤ 24/E_0 ≤ 24.

  The conclusion is unchanged.
- **Last paragraph of §2.** The fat regime φ ≳ √τ·2^{−τ2/4} forces τ2 ≳ 4L > τ*. Correct.

## (e) Theorem S

- **Normalisation.** Correct. With an L²(dA)-orthonormal basis, Z_ℓ = (2ℓ+1)P_ℓ/4π and ĉ_ℓ = 2π∫P_ℓ, so the ℓ = 0 term
  is N_τ s(φ) = E.
- **Self-adjointness.** Follows from Λ_τ = Λ_τ⁻¹.
- **Vanishing off H_ℓ^U.** Follows from Λ_τU = Λ_τ, and UΛ_τ = Λ_τ gives image ⊂ H_ℓ^U.
- **Commutativity (fix the argument).** "The Hecke algebra of the tree is commutative" is not the right reason: U has
  order 24 and does not act transitively on spheres of radius ≥ 5, so C[U\Γ/U] is larger than the radial algebra. What
  holds is the tree recursion
  - 1_{Γ1}*1_{Γj} = 24·1_{Γj+1} + 48·1_{Γj−1} for j ≥ 2;
  - 1_{Γ1}*1_{Γ1} = 24·1_{Γ2} + 72·1_U.

  So every S_τ is a polynomial in A_1 and 24P_U, and these commute.
- **Centred-moment formula.** Correct as a formal identity. The 2k-fold series of an indicator's harmonic expansion is
  not absolutely convergent, so state it for a smooth or Poisson-kernel majorant, or with a summation method.
- **Overstated consequences.**
  - "Weights are squares, so all cancellation must come from the eigenvalues" ignores the signs of ĉ_ℓ. It is true only
    for a positive-definite majorant.
  - "k = 1: variance ≍ E" is not shown. What (R) and Plancherel give is Var = (4π)⁻²Σĉ_ℓ²Σλ_F² ≤ Λ²s(φ) ≤ Cτ²E, with
    Λ ≤ 24Σ_j(j+1)2^{j/2}. Also restore the (4π)⁻² factor.
  - The "Hecke multiplicativity" sentence is vague: multiplicativity relates λ_F(a) and λ_F(b) for one F, not products
    over different F_i.

## §4 heuristics

- **"(U_1) holds unconditionally."** The ball-count sketch is incomplete.
  - Ball counts at radius c·2^{−τ/2} are indeed 2^{o(τ)}: the affine 3-wedge σ1σ2 ≤ Cr³R′⁶, then E5.
  - But F_2 ≲ E·Σ_{dyadic r ≥ φ}(φ/r)²B_τ(r) needs B_τ(r) up to r ≈ 2^{−τ/3}. The trivial covering loses 2^{τ/6} there.
    NOTES.md records clusters of size 2^{τ/9} at the mean spacing: harmless here, but not covered by the sketch.
  - Replace the sketch by "(U_1) holds under (R): Plancherel, 2^{η} = O(τ²)".
- **"Rate > 5/2 needs (U_5)."** True only for this route: write "via Theorem M".
- **"Integrality of v′×w′ forces exact parallelism only below R′⁻⁴."**
  - The quotients have numerators of size R′², so v′×w′ gives |σ1σ2| ≤ CR′⁸·angle. That forces parallelism only below
    ≍ R′⁻⁸.
  - The sharper thresholds come from the numerators themselves (see (b)): 3-wedge gives an exact arithmetic coset below
    R′⁻⁶; 4-determinant gives a common K-hyperplane below R′⁻⁴.
  - State these thresholds precisely.
- **The remaining §4 statements** are reasonable heuristics: the j = 3 description via W1·T_V, CM tori, and exact
  commuting giving harmless linear clusters.

## Numerics

**Code.** Sound.
- SU(2)→SO(3) map, H, S and T are checked.
- The words are T^e(HT|SHT)^d C, built by left multiplication, and the count equals N_τ. My hash check finds all
  72·2^τ−48 words distinct in PU(2) for τ ≤ 16.
- Chord 2sin(φ/2) with s(φ) = E/N_τ. The 3×3×3 grid with cell size = the largest chord is correct. Float storage is fine
  to τ = 22.
- The control is binomial; its ratios are 1.00 ± 0.02.

Minor points:
- The x1 seeds collide (mode·17 + id for nt = 24), so WORDS and CONTROL share some x1 streams. They are correlated but
  not biased.
- `momtable.py`'s "radius discretisation" comment is wrong. There is none: for every fixed x2 the x1-mean is exactly E.
- No per-x2 output is kept, so there are no error bars.

**The author's reading is contradicted by the data.** I reproduced the x2 sequence (same RNG) and attributed Σ N¹⁰ per
x2 at E = 1.

- **τ = 10.** One x2 carries **76%** of the order-10 moment (per-x2 M10 = 342× Poisson).
  - Its maximal cap holds 13 words with T-counts {4,5,5,6,6,…,10,10}.
  - All 12 quotients W0⁻¹Wi share one axis to 1.5·10⁻⁸ rad; the angles are irrational multiples of π; x2 lies 0.23φ
    from the axis.
  - This is a linear cluster W·h^j, j ∈ J, of a hyperbolic h: sizes 2(τ−m)+1, which is why the histogram excess sits at
    the odd values 9, 11 and 13.
  - The next x2's show the same pattern: sizes 11, 9 and 9, and one {4,4,6,6,8,8,10,10} with translation length 2.
- **τ = 15.** One x2 carries **91%**: a 21-word cluster {5,6,6,…,15,15}. My independent x1 draw gives M10/Poisson = 16.4,
  against 10.78 in `mom_tau15.txt`.
- **Clifford axes.** No sampled x2 comes within 2.4φ (τ = 10) or 48φ (τ = 15) of a Clifford face axis.
  - Their 8-groups contribute ≈ (2/3)E²8^{2k−1}/N_τ. At τ = 10, 2k = 10 that is about 1% of Bell(10), two orders of
    magnitude below the observed excess of about 3.
  - So the "x2 near Clifford axes, groups of 8" explanation is wrong. The deviations are CM-torus / P-unit linear
    clusters, i.e. exactly the "harmless" exact-commuting structures of §4.
- **Rigorous lower bound for one axis.** For the single axis of h = HT, I computed Σ_W |W⟨h⟩∩Λ_τ|^{2k−1} exactly
  (maximal cluster 2τ+1). The cap of radius φ/2 around ±axis alone contributes ≥ 0.23, 0.15, 0.067, 0.025 and 0.008 of
  Bell(10) at τ = 8, 10, 12, 14 and 16. Yet the chance that any of the 100 x2 even enters that cap is 6.8·10⁻⁴ at
  τ = 10 and 1.1·10⁻⁵ at τ = 16. There are about 2^s such axes at each tree level s.
- **Sampling limitation.**
  - An x2-event of measure below 1/100, or an (x1,x2)-event below 1/(2·10⁷), is typically invisible.
  - Yet a count n = 2τ+1 = 41 at τ = 20 doubles the order-10 moment already at probability ≈ 10⁻¹¹.
  - The flat 1.00 ratios at τ ≥ 18 are therefore what one would see whether or not the excess is there. The sequence is
    also non-monotone (τ = 15 at E = 1, τ = 12 at E = 4).
  - The exact single-axis computation suggests the *expected* linear-cluster excess does decay roughly like 2^{−τ}·poly(τ).
    It is bounded by poly(τ) anyway, so it is harmless for (U_k). But the data neither show convergence nor probe the
    open part (near-commuting).
- **Coverage gaps.**
  - E ∈ {4, 16, 64} is redundant given Prop. 3, and raw moments at E = 64 are insensitive to clusters.
  - The regime E ≪ 1 that Theorem M actually needs is not sampled.
  - Numerics at τ ≤ 22 cannot test a 2^{o(τ)} hypothesis except for gross failures.
- **Recommendations.**
  1. Output per-x2 power sums and bootstrap them.
  2. Report factorial moments F_j/E (the quantities in (U_k)) at E ≤ 1, including E ≪ 1.
  3. Importance-sample x2 near the axes of short elliptic and hyperbolic elements and their short conjugates, with
     exact weights (`moments_ax.cpp` is a start).
  4. Compute the linear-cluster sums exactly, as above.
  5. Rewrite the §5 reading accordingly.

## Safe as a paper remark

- **Safe.**
  - The exact Hopf form, Lemma 1 and Lemma 2 (with fixes).
  - Theorem M under (U_k), together with Prop. 3 and its corollary (with the τ2 fix), stated as conditional.
  - The implication "tube-H ⇒ (U_k)".
  - (U_1) under (R) via Plancherel.
- **Optional.** Theorem S. It is correct after the fixes, but it adds little.
- **Not safe.**
  - "5/2 already unconditional".
  - The §3 and §4 cancellation and "needs" claims, and the R′⁻⁴ parallelism claim.
  - Any numerical claim of convergence to Poisson or about Clifford-axis 8-groups.
