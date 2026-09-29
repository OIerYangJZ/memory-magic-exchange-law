# Diophantine dichotomy: α ≥ 17/7 in the full CW model (research notes, 2026-09-29)

**Result.** Let a₀ = 28/17. For every frame (G₁, G₂) and every T-count t ≤ (4/a₀)L = (17/7)L (so that
ε = 2^{−L} ≤ CR'^{−a₀}),

    #X_t ≤ R'^{T+o(1)} = 2^{(T/4)t + o(t)},   T = 11183/6800 = 1.644559… < a₀ = 1.647059….

So the per-T-gate count exponent is T/4 = 0.41114 < 7/17, valid up to t = (17/7)L. The generalised Gibbs lemma gives the rate
min(1/κ, 4/a₀) with κ = T/4. For Theorem elem this was min(1/0.4513, 20/9) ≥ 11/5; here it is
min(2.4323, 17/7) = **17/7 ≈ 2.4286**. The same certificate at
a₀ = 5/3 gives α = 12/5 with a larger margin: T = 41/25, worst case 1.64 against a₀ = 1.6667.
Previously we had 11/5 as a theorem and 2.23 as a remark. The numerical limit of the lemma set below is
≈ 2.43, and the level-1 ceiling of the whole box method is still 5/2.

Everything below is proved except where marked. The exact rational certificate is
`cert_dioph.py`; its outputs are `cert_17_7_output.txt` and `cert_12_5_output.txt`. `fast.py` is the float
optimiser used to find the parameters. `cond_generic.py` is the conditional generic-frame variant
(§9).

Notation is as in App. E of the paper:
- R' ≍ 2^{t/4} is the σ₁-radius. Exponents are log base R'. ε = R'^{−a₀}.
- Boxes have transverse side ρR' = R'^{1−a} and along-length λR' = R'^{1−b}.
- A sphere section Y = S³_{R'} ∩ H has σ₁-radius r = R'^{1−g} and σ₂-radius r₂ ≤ 2R'.
- O is the maximal order and O* = {x : ⟨x, O⟩ ⊂ O_K}. Then 4O ⊂ O* ⊂ O_K⁴ ⊂ O.
- For a vector v, **h(v) := |σ₁v|·|σ₂v|**; it is invariant under units.
- "Minimum" means a successive minimum of a Z-lattice with respect to a symmetric convex body.
- "K-independent" means independent over K = Q(√2).
- C denotes absolute constants. All exponent inequalities hold with fixed margins, so constants are absorbed.

---

## D0. Height of a hyperplane
Let H be a K-rational affine hyperplane, H₀ its direction, n ∈ O* a primitive normal (O_K n = Kn ∩ O*), and
M := O ∩ H₀.
- (i) covol_{R⁶}(M) ≍ h(n).
- (ii) M* := {ℓ ∈ H₀⊗K : ⟨ℓ, M⟩ ⊂ O_K} has covol ≍ 1/h(n).

*Proof.*
- (i) Every functional O → O_K is ⟨x, ·⟩ for a unique x ∈ O*. If n is primitive then f_n(O) = O_K, because O_K is a PID and f_n/α ∈ O* for any α with f_n(O) = (α). So 0 → M → O → O_K → 0 is exact. The quotient lattice is O_K scaled by diag(1/|σ₁n|, 1/|σ₂n|) in the plane spanned by (σ₁n, σ₂n). Therefore covol(O) = covol(M)·covol(O_K)/h(n).
- (ii) Under the trace form Σ_σ⟨σx, σy⟩, the Z-dual of M is δ^{−1}M* (δ = 2√2 generates the different), and covol(M^∨) = 1/covol(M). ∎

## D1. Level 1, Dirichlet form
Assume a ≥ a₀, a ≥ b and a₀+a ≥ 2b. Then the points of X_t in a box lie on at most C(1 + R'^κ) parallel
rational hyperplanes, where κ = max(0, 2 − a/2 − 3b/4).

*Proof.* The box lies in a 4-box with sides (2ρR', 2ρR', λR', Cλ²R'). This is Lemma E1's geometry; the last side uses
a₀+a ≥ 2b. Apply Minkowski's first theorem in O* ⊂ R⁸ to {|⟨x₁, e_i⟩| ≤ T₁/d_i} × {|x₂| ≤ T₂}. This gives
ℓ ≠ 0 with T₁T₂ ≍ (Πd_i)^{1/4} = R'^{1−a/2−3b/4}. For w, w' in the box, ⟨ℓ, w−w'⟩ ∈ O_K with
|σ₁| ≤ 8T₁ and |σ₂| ≤ 4T₂R'. An O_K-box [X]×[Y] contains ≤ 1 + CXY points: rebalance by a unit and use
|N ξ| ≥ 1. ∎

For 2a+3b > 8 this is Lemma E1 (κ = 0). The box count is E₁ := 2(a−a₀) + b + κ.

## D2. Level 2 with height: per-cell fibration
Let Y = S³ ∩ H have class g and height h(n) = R'^η. Take cell exponents with X_S ≤ X_L < 1−g and 2X_L − 1 ≤ X_S.
The points of X_t ∩ Y in one cell lie on at most

    1 + C·R'^{P},   P = (X_L + X_S + min(X_S, 2X_L − (1−g)) − η)/3 + 1

circles, hence on ≤ R'^{o(1)}(1 + R'^P) points. With Lemma E4, the level-2 cost of Y ∩ B is at most
N(X_L, X_S) + max(0, P).

*Proof.* The proof of E3 shows three things:
- the cell Q lies within Cx_S of a segment of direction u and length x_L;
- Q ∩ Y lies within Cx_L²/r of the tangent plane at any of its points (normal N ∈ H₀);
- hence Q∩Y − Q∩Y ⊂ D_Q := {v ∈ H₀ : |⟨v,u⟩| ≤ Cx_L, |P_{u⊥}v| ≤ Cx_S, |⟨v,N⟩| ≤ Cx_L²/r}.

Volume of D_Q: 3-sections of a 4-box have volume ≤ 4·(product of the 3 largest sides), and 2-sections have
area ≤ 6·(product of the 2 largest sides). With Fubini across the slab this gives
vol₃(D_Q) ≤ C x_L x_S min(x_S, x_L²/r).

Minkowski in M* with the body T₁D_Q° × B(T₂): John's theorem gives vol(D°) ≥ c/vol(D), so there is ℓ ≠ 0 with
T₁T₂ ≍ (vol D_Q / h(n))^{1/3}. For w, w' ∈ X_t ∩ Y ∩ Q we have w−w' ∈ M, so ⟨ℓ, w−w'⟩ ∈ O_K with
|σ₁| ≤ T₁ and |σ₂| ≤ 2T₂r₂. That leaves ≤ 1 + CT₁T₂R' values. Each fibre is Y ∩ (a rational affine plane), a circle,
and Lemma E5 bounds its points. ∎

At η = 0 and P < 0 this is Lemma E3 exactly: one value, one circle. So D2 is a smooth continuation of E3, and the
coplanarity threshold moves by η.

## D3. Patch fibration
If 2b ≥ a, then the points of X_t ∩ Y ∩ B lie on ≤ C(1 + R'^{(S₀+2U₀−η)/3+1}) circles.

*Proof.* By Lemma E2 the core parameters of Y ∩ B lie in ≤ 2 intervals of length ≤ e_s/R'. Fix one of them.
- Its points lie in a 4-box in the frame (C(s_k), C'(s_k), Π^⊥) with sides (δ_C, e_s, 2e_u, 2e_u).
- Here δ_C ≤ e_s²/R' + Cεe_u ≤ Ce_u; this uses λ² ≤ ρ (i.e. 2b ≥ a) and εe_u ≤ e_s.
- If r ≤ R'/2, then |⟨ν, C(s_k)⟩| ≥ 0.8.
- A central section of a symmetric convex body satisfies vol_{n−1}(K∩ν^⊥) ≤ n·vol_n(K)/w_ν(K). Hence
  vol₃((B₄−B₄)∩H₀) ≤ C e_s e_u².
- If r > R'/2, the three largest sides already give this.

Then argue as in D2. ∎

## D4. Crossing depth
Fix a normal n with unit σ₁-direction ν, and put A := R'|P_Πν|. The sphere with offset c has
**Δ := A − c** (the maximum of h(s) = ⟨ν, R'C(s)⟩ − c). Let δ' := 2ε²R' + 2εr, of exponent
lde = max(1−2a₀, 1−g−a₀), and δ'_B := the same with ε → ρ, of exponent ldb.

- **Tangent class** (Δ ≤ C'δ'). Along the tube the core length satisfies e_s^T ≤ C√(R'δ'). In a box it is bounded by the
  E2 value. The number of offsets c per normal is ≤ 1 + C h(n) δ' R'.
- **Crossing class** (Δ ∈ [Δ₀, 2Δ₀], Δ₀ ≥ C'δ'). Then e_s^T ≤ Cδ'√(R'/Δ₀), and in a box
  e_s^B ≤ C min(λR', √(R'δ'_B), δ'_B√(R'/Δ₀), r). The number of offsets is ≤ 1 + C h(n) Δ₀ R'. Always Δ ≤ r²/R'.

*Proof.*
- The sublevel set {|h| ≤ δ'} consists of the two preimages of an interval of length 2δ' under
  u ↦ A(1−cos u), on u ∈ [0, π/3].
- The derivative there is ≥ 0.8·A·u, with u ≥ √(2(Δ−δ')/A), which gives the crossing bound. The tangent bound is E2.
- For a box, |Δ_B − Δ| ≤ R'φ₀²/2 + R'φ₀|P_{Π⊥}ν| ≤ Cδ', so Δ_B ≥ Δ/2 in the crossing class.
- Offsets: σ₁(⟨n,w⟩) = |σ₁n|·c ranges over an interval of length |σ₁n|·(Δ-range), and σ₂ over one of length 2|σ₂n|R'.
  Count in the O_K-box. ∎

**Per-sphere count** (sphere of class g and height ≥ R'^{η_c}). The number of points of X_t ∩ Y in the tube is
at most the minimum of two bounds:
- **tube bound:** R'^{max(0, (S_T + 2U_T − η_c)/3 + 1)}, by D3 applied to Y ∩ tube, with
  S_T = min((1+lde)/2, lde + (1−D)/2, 1−g) and U_T = min(1−a₀, 1−g);
- **box bound:** #boxes × (per-box cost), with
  #boxes ≤ R'^{max(0, S_T−(1−b)) + 2max(0, U_T−(1−a))}, and the per-box cost from D2/D3 with the D4 value of S₀.

## D5. Low-height normals near Π (successive minima with wedge constraints)
Let γ := min(g, a₀) and h := R'^η. The primitive normals (modulo units) of class-g spheres that meet the tube and have
h(n) < h number at most C·R'^{max(η, 2η, 3η−γ, 4η−2γ)}.

*Proof.*
1. The spherical centre of Y is within ψ_Y + 2ε ≤ Cψ of a point of the core, with ψ := R'^{−γ}. So a balanced
   representative of n lies in 𝒦 = {|P_Πσ₁x| ≤ C√h, |P_{Π⊥}σ₁x| ≤ Cψ√h, |σ₂x| ≤ C√h}.
2. Let μ̃_j be the "K-minima": μ̃_j is the least λ such that λ𝒦 contains j K-independent vectors of O*. Then
   μ̃_j ≤ λ_{2j−1} ≤ λ_{2j} ≤ 2.42μ̃_j, using multiplication by 1+√2.
3. Henk's bound gives #(𝒦 ∩ O*) ≤ CΠ_j max(1, μ̃_j^{−2}).
4. Wedges of vectors of O* have O_K coordinates, so a nonzero wedge has h ≥ 1. Write σ₁v_i = p_i + q_i with p_i ∈ Π and
   |q_i| ≤ Cψ|σ₁v_i|. A 3-wedge then has σ₁-norm ≤ Cψ·Π|σ₁v_i|, because p₁∧p₂∧p₃ = 0. A 4-determinant has
   σ₁-norm ≤ Cψ²·Π|σ₁v_i|.
5. Hence μ̃₁ ≥ ch^{−1/2}, μ̃₁μ̃₂ ≥ c/h, (μ̃₁μ̃₂μ̃₃)² ≥ c/(ψh³) and (μ̃₁…μ̃₄)² ≥ c/(ψ²h⁴).
6. If exactly k minima are < 1, the count is ≤ C(μ̃₁⋯μ̃_k)^{−2}. ∎

For η ≤ γ this is R'^{2η}. The worst case is a rational plane V near Π that carries h² lattice points.

## D6. Rank-2 Diophantine dichotomy
In the regime k = 2 (μ̃₃ > 1), all lattice points of 𝒦 lie in the K-plane V spanned by the first two minima.
Let H_V = h(primitive generator of Λ²(V∩O*)) ≥ 1, and let θ₁ ≤ θ₂ be the principal angles between V₁ = σ₁(V)⊗R
and Π. Then:

- **(a) Sector count.** #(𝒦 ∩ V) ≤ C max(1, h, β h²/H_V), with β = min(1, ψ/sin θ₂).
  A unit vector of V₁ within ψ of Π makes angle |sin t| ≤ ψ/sin θ₂ with the principal direction e₁. So
  σ₁(𝒦∩V) lies in a rectangle √h × β√h. Two K-independent vectors in μ𝒦 have
  h(v₁∧v₂) ≤ C μ₁²μ₂² β h² and h(v₁∧v₂) ≥ H_V.
- **(b) Global fibration.**
  #X_t ≤ R'^{o(1)}(1 + C H_V R'⁴ (sin θ₂ + ε)²).
  - Take a reduced basis ν₁, ν₂ of V^⊥ ∩ O*. Its covolume is ≍ H_V by lattice duality under the trace form, since
    4O ⊂ O* ⊂ O. By Minkowski II, h(ν₁)h(ν₂) ≤ CH_V.
  - Every unit vector of V₁^⊥ is within sin θ₂ of Π^⊥. So ⟨ν_i, w⟩ over the tube takes
    ≤ 1 + C h(ν_i) R'²(sin θ₂ + 2ε) values. The fibres are circles (E5), and εR'² ≥ 1.

**Dichotomy.** Fix F* = T.
- If (b) is ≤ R'^{F*} for some (g-interval, layer), then #X_t ≤ R'^{T+o(1)} and we are done.
- Otherwise X := log H_V and τ := log sin θ₂ satisfy X + 4 + 2max(τ, −a₀) > F*, with X ≥ 0 and τ ≤ 0. Maximising the
  exponent of (a), 2η − X − max(0, τ+γ), over this region gives 2η − G with **G = max(0, γ + F*/2 − 2)**.
  The maximum is at τ = (F*−4)/2, X = 0.

So in D5 the term 2η can be replaced by 2η − G, giving

    normals(η) = max(η, 2η − G, 3η − γ, 4η − 2γ).

## D7. Assembly (HIGH/LOW split with height layers)
Fix a g-interval [g_l, g_r] and layers 0 = η₀ < η₁ < … < η_m = η*.
- **HIGH.** For a box, take its ≤ R'^κ hyperplanes whose sphere is of class g and height ≥ R'^{η*}. They contribute
  ≤ R'^{E₁ + cost_{η*}(g)}, where the cost is the minimum of D2 and D3.
- **LOW.** Spheres of class g with height in [R'^{η_{j−1}}, R'^{η_j}) are counted globally, over all boxes. Their
  number is
  ≤ R'^{normals(η_j)} × (offsets, D4 with h ≤ R'^{η_j}) × (per-sphere count, D4 with height ≥ R'^{η_{j−1}}),
  maximised over the tangent class and the crossing classes. The crossing classes are a finite partition of
  D ∈ [lde(g_r), 1−2g_l]; offsets are taken at the right end and per-sphere counts at the left end.
- Every point of X_t lies in some box on one of its hyperplanes (D1), so it is counted by HIGH or by LOW.

Hence

    log #X_t ≤ max( zero case,  max over intervals of min_{η*} max(HIGH, LOW₁, …, LOW_m) ) + o(1),

or ≤ F* if some fibration (D6b) succeeds. Two regimes need separate treatment:
- g > 2.01: whole spheres are coplanar, so the cost is 0 and the exponent is E₁ ≤ T.
- The zero case (r > R'/2, g < 1/den): HIGH with S₀ = 1−b.

**Monotonicity (for the grid).** S₀, U₀, S_T, U_T, lde and ldb are non-increasing in g, and N is non-decreasing in them,
so these are evaluated at g_l. The sagitta term min(X_S, 2X_L−(1−g)) and the constraint X_L < 1−g are evaluated at g_r.
normals is non-increasing in γ, so γ = min(g_l, a₀). Offsets are increasing in D and per-sphere counts decreasing, so
each crossing class is treated at its end points as stated above.

## Certificate
Run `python cert_dioph.py 28/17 28/17 157/100 11183/6800 400 100 12`. The parameters are:
- a₀ = a = 28/17, b = 157/100, κ = 0, E₁ = 1.57;
- 804 g-intervals of width 1/400 up to 2.01;
- η layers of width 1/100 and 12 crossing classes.

All inequalities are exact over Q. Floats are used only to propose cell witnesses (X_L, X_S), which are then
rationalised and re-checked. Worst exponent 67097/40800 = 1.644534 ≤ T = 1.644559 < a₀ = 1.647059, at g ≈ 1.50 and
η* = 0.91; max η* = 0.92. The run takes seconds.

Regression check: with η = 0 the new level-2 cost at the paper's parameters (a₀ = 9/5, a = 91/50, b = 73/50) is
sup_g ≤ 0.300, against 0.30501 in `research/cert`. This is consistent, since D2/D3 only add options.

## What each ingredient buys (float optimiser, approximate)
| model | α |
|---|---|
| paper lemmas E1–E5 | 2.234 |
| + D2/D3 at η = 0 | 2.23–2.24 (the worst case moves to the tangent rod, g ≈ 1.15) |
| + HIGH/LOW split, crude normal count 2η + (η−γ/3)₊ + (η−γ/2)₊ | ≈ 2.30 |
| + height layers and crossing classes (D4) | ≈ 2.37 |
| + K-minima normal count (D5) | ≈ 2.39 |
| + rank-2 dichotomy (D6) | **≈ 2.43 (certified 17/7)** |

## Where it stops, and conditional variants
- The binding pair is HIGH = E₁ + cost_η at g ≈ 1.45–1.55 (small tangent spheres, patch fibration) against the
  rank-2 normal count 2η − G.
- Level 1 caps everything at 5/2 (E₁ ≥ (8−2a₀)/3).
- **Generic frames (conditional, `cond_generic.py`).** If the low-height normal count is at most the volume term
  R'^{max(0,4η−2γ)} — i.e. Π has no unusually good rational approximations at the relevant scales, which holds for
  Haar-almost-every frame by Borel–Cantelli (not written out) — the same lemmas give α ≈ 2.49, essentially the level-1
  ceiling. But frames in the CW model are chosen by the process (e.g. Clifford+T prefix/suffix), so this is not a
  worst-case statement.
- Not done:
  - the rank-3 dichotomy (a single rational functional near Π^⊥, then a band problem on a great sphere);
  - any improvement of D2/D3 beyond the Minkowski exponent 1/3. The heuristic count is e_s e_u r₂/(r h), far below it.

## Status
- **Fully proved** (modulo writing at paper standard): D0–D7 and the certificate. The only external inputs are
  Minkowski I/II, Henk's lattice-point bound, John's theorem, and the paper's Lemmas E2, E4 and E5.
- Nothing in the 17/7 bound is conditional.
- Worth an independent check before any paper use:
  - the sector bound D6(a);
  - the duality covol(V^⊥ ∩ O*) ≍ H_V with the O vs O* normalisations;
  - that LOW is summed over all boxes without the (1 + R'^κ) factor. It is: spheres are counted once each.

## Errata after independent review (verdict: VALID WITH FIXABLE GAPS; α ≥ 17/7 holds; scripts in research/verify4/)
- **D0.** Verified exactly on the explicit order: covol O = 4, covol O* = 1024 and 2O ⊂ O*. covol(O ∩ n^⊥) = √2·h(n) for primitive n, and covol(M*) = 512/covol(M).
- **D1.** At 17/7 we have κ = 0 (2a+3b = 8.004), so D1 reduces to Lemma E1. The remark "level 1 caps at 5/2 (E₁ ≥ (8−2a₀)/3)" is false once κ > 0, since E₁ can go down to 2 − a₀/2. It is a heuristic ceiling, not a bound.
- **D3.** Replace δ_C by max(δ_C, e_u). For r ≤ R'/2 the section-width bound vol₃ ≤ 4vol₄/w_ν already gives C e_s e_u², so 2b ≥ a is needed only for r > R'/2. D3 holds for any patch region, including Y ∩ tube (used in LOW).
- **D5, step 4.** The claim |q_i| ≤ Cψ|σ₁v_i| is false for short vectors. Use instead |p_i| ≤ Cμ_i√h and |q_i| ≤ Cμ_iψ√h, which is what step 5 needs.
- **D6(a).** The proof should be phrased in absolute terms: P_{Π⊥}f₁ ⊥ P_{Π⊥}f₂ gives |u| ≤ Cψ√h/sin θ₂. In D6(b), read (sin θ₂ + 2ε). The dual covolumes satisfy covol(V^⊥ ∩ O*) = covol(V ∩ O*) = 8H_V, checked exactly.
- **D7 / rate.** τ_* = ⌈(17/7)log₂K⌉ can exceed (17/7)L by O(1). There is no slack in a₀, so this is absorbed by allowing ε ≤ CR'^{−a₀} with bounded C, which all lemmas tolerate.
- **Certificate.** "worst" is the first η* that succeeds, not the optimum. An independent pointwise re-implementation gives a sup of 1.6347 (g ≈ 1.47, η* ≈ 0.92) against a₀ = 1.6471, and 1.6255 against 1.6667 at the 12/5 parameters. The float limit of the lemma set is ≈ 2.435.
- **Practical.** Purely asymptotic. The margins are tiny (≈ 0.0025), so L₀ is astronomically large.
