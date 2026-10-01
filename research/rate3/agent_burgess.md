# Burgess / bilinear forms on the congruence reformulation (agent report, 2026-10-01)

**Verdict: fails.** Neither Burgess's method nor a Type I/II (dispersion) treatment of a_χ gives
|S_χ| ≤ T₀N(H)^{−1/2−δ}, even in a restricted range of z₀ that is not already trivial. There are three separate
obstructions, each located precisely:
- **(B1) Range.** Drop the weight a_χ(Hn − N(g)) and keep only the pure character sum over the box. At the critical
  scale this model sits *exactly* at the Burgess endpoint. Its Poisson dual has length q^{1/4}, which equals the volume
  E. Every Burgess exponent r gives a bound ≥ √q, with equality only at r = 1 (Pólya–Vinogradov).
- **(B2) Weight.** The summation set is the lattice points of a definite quadric, which has no translation symmetry.
  Burgess's collection step needs the summand at g + ab to depend only on (āg mod H, b). The weight depends on
  N(g + ab) = N(g) + Tr(ḡab) + N(a)N(b).
- **(T) Bilinear.** For χ² ≠ 1, a_χ is a dihedral *cusp* form coefficient, so it has no divisor structure over O_K.
  The kernel 1[N(m) + N(g) = h] is non-negative and all the oscillation sits in completely multiplicative
  coefficients. Every Type II bound is therefore trivial. Dispersion needs an average over the shift h = Hn or over
  the modulus H, and neither exists.

Moreover, even the baseline θ = 1/2 for individual S_χ is **not known**. The determinant method bounds N(κ′) for
every κ′, and that says nothing about a single S_χ (§2).

## 1. The reformulation: checked

Write ⟨x, y⟩ = x₁ȳ₁ + x₂ȳ₂ on Z[ω]², with ω = ζ₈. Put f = (−p̄, q̄) and f^⊥ = (q, p). Then ⟨f, f^⊥⟩ = 0 and
|f|² = |f^⊥|² = H = |p|² + |q|². For s = (u, t):
- **Coordinates.** g = ⟨s, f⟩ = qt − pu and m = ⟨s, f^⊥⟩ = uq̄ + tp̄. So Hs = gf + mf^⊥ and H|s|² = |m|² + |g|²
  exactly.
- **Integrality.** s ∈ Z[ω]² ⟺ −gp̄ + mq ≡ 0 and gq̄ + mp ≡ 0 (mod H).
  - (⇒) Multiply by a and −b and add. With qa − pb = 1 this gives m ≡ κg, where κ = q̄b + p̄a.
  - (⇐) κq = p̄ + Hb ≡ p̄ and κp = Ha − q̄ ≡ −q̄.
- **Consequence.** s ↦ (λ, g) with m = Hλ + κg is a bijection Z[ω]² → Z[ω]², and
  |Hλ + κg|² + |g|² = H·|s|². ✔
- **Bonus identity.** κκ̄ + 1 = H(|a|² + |b|²). So κ lies on the conic xx̄ ≡ −1 (mod H), as the character expansion
  needs.
- **Numerical check (exact arithmetic in Z[ζ₈]).** 168 random (z₀, X, s) instances, including the converse map, all
  pass. So do all 3·2 sets of genuine tube points (all 2.1·10⁶ vectors of norm 2⁸, three random z₀, two ε). Every
  tube point lands in the g-box with m ≡ κg.
- **Primitivity.** In 32 of 200 trials the Minkowski vector (p, q) was not Z[ω]-primitive. Divide out d = gcd. Then
  T₀ and N(H) both drop by N_{Q(ω)/Q}(d), so E is unchanged.

**Box.** Since g = q(t − z₀u) + (qz₀ − p)u, we get |σ₁g| ≤ Y₁ := XεR + cRX^{−3} and |σ₂g| ≤ Y₂ := √2·XR.
- The box is the product of discs {σ_i N(g) ≤ Y_i²}, so it is a condition on ν = N(g) alone.
- #Box ≍ Y₁²Y₂² when Y₁Y₂ ≫ 1. When Y₁Y₂ < 1, Box = {0}, since N_{Q(ω)/Q}(g) ≥ 1 for g ≠ 0.
- Balancing by a unit (1+√2)^l turns it into a box of side (Y₁Y₂)^{1/2}.
- In σ₂ the condition is vacuous: Y₂² ≍ σ₂(Hn).

**Character expansion (bookkeeping).**
- *The coset.* For (g, H) = 1, the ratio x = m/g lies in the coset C = {xx̄ = −1} of T = {xx̄ = 1} ⊂ (Z[ω]/H)^×.
  Also κ ∈ C, and κ is a unit mod H.
- *Expansion.* Extend each χ ∈ T̂ arbitrarily to χ̃ on (Z[ω]/H)^×. Then
  N^{cop} = |T|^{−1} Σ_χ χ̃̄(κ) S_χ. Each product χ̃̄(κ)S_χ is independent of the extension.
- *Vanishing.* Replacing m by ωm shows S_χ = 0 unless χ(ω) = 1.
- *Non-coprime g.* Dividing (m, g) by the common factor reduces this part to the same problem with a smaller modulus.
  I did not track it, and it plays no role below.
- **S_χ does not depend on κ.** It depends on z₀ only through (H, Box). The adversary choosing z₀ chooses κ (any
  point of C) and essentially chooses H (for example, H prime is allowed).

**Automorphic form of S_χ.** Since the box is a condition on ν = N(g),
S_χ = Σ_{ν ∈ Ω} a_χ̄(ν) a_χ(h − ν), with h = Hn and Ω = {σ₁ν ≤ Y₁², σ₂ν ≤ Y₂²}.
- Up to a unimodular factor depending only on ν, a_χ(ν) = 8·Σ_{N𝔄=(ν)} ψ(𝔄). Here ψ is a Grössencharacter of Q(ω)
  of conductor dividing H. The infinity type absorbs χ(1 + √2).
- So a_χ is the coefficient of the weight-one CM Hilbert form θ_ψ over K. Its level has norm ≍ 4N(H)² ≍ q.
- **S_χ is therefore a shifted convolution of two weight-one dihedral forms of level ≍ q, at the single shift h, with
  ν confined to a window of relative size ε² in σ₁ and the full range in σ₂.**

**Exponent table.** Write ε = R^{−a} and X⁴ = 1/ε. Exponents of R:

| quantity | general a | a = 8/5 | in terms of q = N(H)² |
|---|---|---|---|
| N(H) | a | 8/5 | q^{1/2} |
| T₀ = #Box | 4 − a | 12/5 | q^{3/4} |
| volume E = T₀/N(H) | 4 − 2a | 4/5 | q^{1/4} |
| cells = 1/ε = PV size √q | a | 8/5 | q^{1/2} |
| target for Q1 | a − κ | 8/5 − κ | q^{1/2−η} |
| √T₀ (square-root cancellation) | 2 − a/2 | 6/5 | q^{3/8} |
| Y₁ (σ₁ radius of the box) | 1 − 3a/4 | −1/5 | – |
| balanced side; period | 1 − a/4; a/2 | 3/5; 4/5 | q^{3/16}; q^{1/4} |
| dual side; dual length q/T₀ | 3a/4 − 1; 3a − 4 | 1/5; 4/5 | q^{1/16}; q^{1/4} |

Dictionary: an error E_κ = T₀N(H)^{−θ} gives rate 2 + θ. At a = 8/5, θ > 1/2 ⟺ error < √q.

## 2. The per-χ target is strictly stronger than what is known

- A bound |S_χ| ≤ T₀N(H)^{−θ} for all χ ≠ 1 implies |N(κ) − main| ≤ T₀N(H)^{−θ} for every κ.
- The converse fails badly. Parseval gives Σ_χ |S_χ|² = |T|·Σ_{κ′} |N(κ′) − main|². So a uniform bound
  |N(κ′) − main| ≤ √q, which is what the determinant method proves (one point per cell), only gives
  RMS_χ |S_χ| ≤ N(H)^{1/2}√q = q^{3/4} = T₀. That is the trivial bound.
- I found no argument giving |S_χ| ≤ √q·R^{o(1)} even at θ = 1/2. Completion in g mod H needs the weight
  F(g) = a_χ(h − N(g)) to have a concentrated Fourier transform. For a rough F, the completion bound is
  (q·Σ|F|²)^{1/2} ≈ (qT₀)^{1/2}, which exceeds T₀.
- So the per-χ program would first have to reprove 5/2 analytically, and then beat it.

## 3. Burgess

**3.1 The weightless model, direct.** Take M_χ = Σ_{g ∈ Box} χ̄(g), with χ primitive mod H and a balanced 4-dim box.
- Assume Burgess extends to Z[ω]/H with the usual exponents, which is standard for boxes over an integral basis.
- Then |M_χ| ≪ T₀^{1−1/r} q^{(r+1)/(4r²)+ε}.
- At T₀ = q^{3/4} the exponent is (3r² − 2r + 1)/(4r²): 1/2, 9/16, 11/18, 41/64, … for r = 1, 2, 3, 4.
- This is ≥ 1/2, with equality only at r = 1, i.e. PV.
- Burgess beats PV only for T₀ < q^{(2r+1)/(4r)} ≤ q^{5/8}, i.e. for a > 16/9. So length q^{3/4} is inside Burgess's
  *nontrivial* range but outside the range where it *improves on PV*, and the target is to beat PV.

**3.2 The model, dual: the critical scale is the Burgess endpoint.** Use a smooth box weight and primitive χ, so
|τ(χ)| = √q. Poisson mod H gives
M_χ = (τ(χ̄)T₀/q)·Σ_ξ χ(ξ) w(ξ) (up to constants from the different),
where w is a bump on the dual box of side N(H)^{1/2}/B = q^{1/16}.
- So |M_χ| = (T₀/√q)·|dual sum of length q/T₀|, and at a = 8/5 the dual length is q^{1/4} = R^{4/5} = E.
- **Hence M_χ ≤ q^{1/2−η} ⟺ a power saving for a character sum of length exactly q^{1/4}.**
- That is the Burgess endpoint. Write the dual length as q^{1/4+β}, with β = (5a − 8)/(4a).
- Burgess applied to the dual gives |M_χ| ≪ √q·q^{(1−4βr)/(4r²)}.
  - For a ≤ 8/5 (β ≤ 0) this is > √q for every r.
  - For a = 8/5 + μ the best r ≈ 1/(2β) saves q^{−β²}, i.e. R^{−(16/5)(25μ/32)²}.
  - Examples: μ = 0.1 gives R^{−0.018}, against R^{−0.136} from the determinant method's 4μ/11 + μ relative to cells.
- **So for the model, the exact set of scales where Burgess-type arguments beat PV is a > 8/5.** That is the range
  where the determinant method already wins, and by much more. At a = 8/5 and below, nothing.
- In L-function terms (Perron, with the O(log) Grössencharacter twists λ that encode the box shape): beating PV at
  length C^{3/4} *would follow from* |L(½ + it, ψλ)| ≪ C^{1/8−η}. Via interpolation with σ = 0, anything weaker than
  exponent 1/8 gives nothing. That is beyond Weyl (C^{1/6}, Petrow–Young over Q), beyond
  Milićević's sub-Weyl for prime-power moduli (≈ C^{0.1645}), and far beyond Burgess-type bounds for Hecke characters
  (Wu).
- Smooth or powerful H (q-van der Corput, Postnikov) do not change this at length q^{3/4}. One A-process gives
  (T₀q^{1/3})^{1/2} = q^{13/24}. A q-analogue of the exponent pair (11/30, 16/30) would beat PV, but only for
  specially factored H, which the adversary does not have to supply.

**3.3 The actual S_χ: where Burgess's argument breaks.**
1. **Shifting.** Shifting the variable is harmless: S_χ = |𝒜ℬ|^{−1} Σ_{a,b} Σ_g χ̄(g + ab) F(g + ab) + O(boundary),
   with F(g) = a_χ(h − N(g)).
2. **Factoring.** χ̄(g + ab) = χ̄(a)χ̄(āg + b) also works.
3. **Collecting (this is where it fails).** Burgess groups by z = āg mod H and applies Hölder to Σ_z ν(z)|Σ_b χ̄(z + b)c_b|.
   That needs the coefficient to be a function of (z, b) only. Here it is F(g + ab), a function of the integer
   N(g) + Tr(ḡab) + N(a)N(b). It is neither a function of z = āg mod H nor of bounded variation.
   - With coefficients depending on (a, g, b), the 2r-th moment cannot be completed, so no Weil bound applies.
   - Arbitrary bounded c_{a,g,b} (for example c = χ(āg + b)) destroy all cancellation.
   - The deeper reason: in (m, g)-space the summation set is the lattice points of the definite quadric
     N(m) + N(g) = h. Translation by v ≠ 0 with m fixed keeps a point on the quadric only if Tr(ḡv) = −N(v), a
     codimension-2 condition. So no approximate translation invariance exists, and the boundary error in step 1 is of
     the size of S_χ itself if one tries to keep the quadric.
4. **Linearising the weight (circle method) does not help.** Writing F(g) = ∫ θ_χ(z)e(−(h − N(g))z)dz turns S_χ into
   ∫ θ_χ(z)e(−hz)·M_χ(z)dz, with mixed sums M_χ(z) = Σ_{g ∈ Box} χ̄(g)e(Tr N(g)z).
   - Burgess for mixed sums (Heath-Brown–Pierce type) bounds sup_z|M_χ(z)| at best like the model, and M_χ(0) *is*
     the model.
   - The outer integral is a *binary* additive problem at one shift. With absolute values one gets at best
     ∫|θ_χ||M_χ| ≈ ‖θ_χ‖₂‖M_χ‖₂ ≈ T₀ (heuristic), i.e. nothing.
5. **The non-abelian replacement.** The only symmetries of the quadric are U(2) rotations. The integral ones are a
   finite group; the rational ones are Hecke correspondences. Amplifying with those is the spectral route, which gives
   the Ramanujan error R² ≫ R^{8/5} (NOTES §3c, §6).

## 4. Type I/II and dispersion

**(i) No divisor structure.** a_χ is the coefficient of AI_{Q(ω)/K}(ψ). Since χ̃^c = χ̃̄ on T, this is cuspidal unless
χ² = 1.
- So for χ² ≠ 1 there is no factorisation a_χ = α ∗ β over O_K with a smooth factor.
- The χ with χ² = 1 (genus characters, where a_χ *is* a divisor sum) number |T[2]| ≤ 2^{ω(H)+O(1)} = R^{o(1)}.
  Bounding them trivially, their total contribution to N(κ) is ≤ R^{o(1)}·T₀/N(H) = R^{o(1)}E. They are harmless,
  and irrelevant.

**(ii) Bilinear forms are trivial here (rigorous, elementary).**
- Group the variables in any way: x = m, y = g; or x = m₁, y = (m₂, g) after m = m₁m₂ in Z[ω]; and so on. Then
  S_χ = Σ_{x,y} α_x β_y K(x, y), where K is a *non-negative* incidence (1[N(m) + N(g) = h],
  1[N(m₁)N(m₂) + N(g) = h], …) and α, β are unimodular (products of values of χ̃, χ̃̄).
- A Type II estimate uses only |α|, |β| ≤ 1 and the kernel. But α = β = 1 is admissible and gives ΣK, the trivial
  bound. So no Type II estimate can go below the trivial bound. Concretely, ‖K‖·‖α‖·‖β‖ ≥ ΣK.
- Moving χ into the kernel changes nothing. K_χ = diag(χ̃)·K·diag(χ̃̄) is a unitary conjugate of K on the coprime
  part, so ‖K_χ‖ = ‖K‖. For example, K_χK_χ* is block-diagonal on {N(m) = const}, and each block is the rank-one
  matrix χ̃(m)χ̃̄(m′) times a count.
- **Bilinear methods need the oscillation in the kernel (e(mn/c), χ(mn + a), …). Here the additive structure is a
  non-negative incidence and the oscillation is a completely multiplicative coefficient, so the two never interact.**

**(iii) Dispersion needs an average that is not there.** BFI/BDH gains come from averaging over moduli or over the
shift.
- *The shift.* h = Hn is a single value, with n = 2^k. Averaging over h is averaging over the norm, which is exactly
  the averaging Marshall uses and that the single prime forbids (NOTES §6).
- *The modulus.* H is fixed by z₀. The Dirichlet body has volume ≍ 1, so generic z₀ have O(1) admissible (p, q) per
  dyadic X.
- *Non-optimal approximations.* Allowing them, |σ₁(qz₀ − p)| ≤ M^{1/2}X^{−3}, gives ≍ M moduli H_i. But
  modulus H_i alone only localises s to the tube of width M^{1/2}ε around the rational circle t = (p_i/q_i)u.
  - Indeed min_i N_i ≤ M^{−1}Σ_i N_i = M^{−1}Σ_s w(s), where w(s) = #{i : s ∈ fat tube_i} is a weight of width
    M^{1/2}ε.
  - So averaging over the M moduli is the same as *fattening the tube* to a′ = a − log M/(2 log R) < 8/5.
  - There the best known bound (determinant boxes, R^{(8−2a′)/3}) is worse than at the critical scale. Nothing is
    gained.

## 5. What does hold (restricted ranges), and what remains

All of the following are elementary. None reaches generic z₀.
- **(a) Low-conductor characters.** Take the characters with conductor H′ | H, N(H′) ≤ N(H)^{1/2−δ}. Their total
  contribution is (|T_{H′}|/|T_H|)·#{m ≡ κg mod H′} ≤ R^{o(1)}T₀N(H′)/N(H) ≤ T₀N(H)^{−1/2−δ}R^{o(1)}. In the model,
  PV already handles every χ with N(cond) ≤ N(H)^{1−δ}. So the whole difficulty is the ≍ N(H) characters of
  essentially full conductor.
- **(b) z₀ near rationals of small height (divisor bound).** Suppose p/q ∈ Q(ω) has |σ_i p|, |σ_i q| ≤ Y and
  |σ₁(qz₀ − p)| ≤ Yε. Then N ≤ Σ_{g ∈ Box} r(Hn − N(g)) ≤ R^{o(1)}(1 + Y⁴ε²R⁴).
  - At a = 8/5 this is ≤ R^{8/5−κ+o(1)} for Y ≤ R^{1/5−κ/4}.
  - The g = 0 term is the R^{o(1)} points on the rational circle.
  - This is the paper's Clifford-frame bound with its 2^h loss. The set of such z₀ has measure ≈ R^{−8/5}.
- **(c) The model above the critical scale.** For a = 8/5 + μ, Burgess on the dual gives a q^{−β²} saving with
  β = 5μ/(4a) (§3.2). It is strictly weaker than the determinant method there and vanishes at μ = 0.

**What remains, stated sharply.**
1. Even for the pure character sum over the box, κ > 0 at a = 8/5 needs a power saving for character sums of
   length q^{1/4}, in a 4-dim box of side p^{1/4} with p = N(H)^{1/2}. This is uniform over moduli H, including prime
   H. It is the Burgess endpoint (it would follow from C^{1/8−η} subconvexity, which is beyond Weyl). No power saving at
   length q^{1/4} is known for general, e.g. prime, moduli.
2. For the true S_χ one must in addition control a binary additive problem (a shifted convolution of two level-q
   dihedral forms) at a *single* shift, over a range of length q^{3/4} below the level. No method handles this
   without averaging over the shift or the modulus, and both averages are structurally unavailable.
3. Coincidences worth recording:
   - E = q^{1/4} is both the volume and the dual length.
   - The model's Burgess threshold (a > 8/5) coincides with the determinant method's threshold.
   - So the "geometric 5/2" and the "analytic 5/2" are the *same* endpoint, and on the analytic side it is the
     q^{1/4} Burgess barrier, not merely PV.
