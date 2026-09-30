# Beyond 17/7: annular fibration and rational projections (research notes, 2026-09-29)

**Result (exact certificate; independent review in `review/REVIEW.md`: A1, P1, P2, P3, G3 and
`cert_proj.py` all VALID, and the write-up gaps it lists are fixed below).** In the full CW model, for every frame,
the tube bound of Theorem `thm:dioph` holds with `a0 = 400/249`, i.e. every bit forgone costs at
least **α ≥ 249/100 = 2.49** committed T gates asymptotically, unconditionally.
With the first new lemma alone (A1), α ≥ 27/11 ≈ 2.4545 (margin a0 − T = 1/300, larger than the
1/400 of the 17/7 certificate).

The lemma set now stops essentially at the level-1 ceiling 5/2 of Lemma E1 (one hyperplane per box),
see "Where it stops" below.

Notation as in App. E/F of `main.tex` and `../dioph/DIOPH_NOTES.md`: exponents are logs to base R',
ε = R'^{−a0}, δ := εR' is the tube radius, Π the core plane, h(v) = |σ1 v||σ2 v|, O* the dual order,
boxes have transverse side ρR' = R'^{1−a} and core length λR' = R'^{1−b}, e_s, e_u the patch
lengths of Lemma E2, γ = min(g, a0), ψ = R'^{−γ}. C is an absolute constant.

| file | what |
|---|---|
| `../dioph/cert_annulus.py` | the 17/7 certificate with G replaced by G3 (A1 only) |
| `cert_proj.py` | exact certificate with A1 + P1 + P2 + P3 |
| `cert_A1only_27_11_output.txt` | `cert_annulus.py 44/27 44/27 317/200 4391/2700 400 100 12` |
| `cert_proj_27_11_output.txt` | same parameters with `cert_proj.py` |
| `cert_249_100_output.txt` | `cert_proj.py 400/249 400/249 1192747/747000 319751/199200 400 100 12` |
| `indep_model_annulus.py` | pointwise float model (A1 only), from `../verify4/indep_model.py` |
| `model2.py` | pointwise float model with all four lemmas; `breakdown*.py` print the binding terms |

All certificates have κ = 0 (2a+3b > 8), so they use Lemma E1 exactly as proved in the paper, not
the relaxed D1.

---

## A1. Annular fibration (replaces F5(b))

**Lemma A1.** Let V ⊂ K^4 be a K-plane of height H_V, and θ1 ≤ θ2 the principal angles between
V1 := σ1(V)⊗R and Π. Then

    #X_t ≤ 2^{O(t/log t)} · C · (1 + (sin θ2/ε)^{1/2} + H_V R'^4 ε (sin θ2 + ε)).

F5(b) had H_V R'^4 (sin θ2 + 2ε)^2: one factor sin θ2 is replaced by ε.

*Proof.*
1. *Fibres.* The standard form is anisotropic over K, so K^4 = V ⊕ V⊥. Put v(w) := P_{V⊥} w. If
   v(w) = v(w') then w' ∈ w + V, an affine K-plane meeting the sphere in a circle; by Lemma E5 it
   carries ≤ 2^{O(t/log t)} points of X_t. (E5 needs a plane spanned by points of X_t: a fibre with ≥ 3
   points is spanned by them, since three distinct points of a circle are affinely independent; fibres
   with ≤ 2 points are trivial. F5(b) of the paper uses the same step silently.) So
   #X_t ≤ 2^{O(t/log t)} · #v(X_t).
2. *Where v(w) lies.* σ1 v(w) = P_{V1⊥} σ1 w. Write σ1 w = R' cos φ C(s) + R' sin φ n' with n' ∈ Π⊥ a
   unit vector and sin φ ≤ 2ε. The map P_{V1⊥}|_Π has singular values sin θ1 ≤ sin θ2 with
   orthonormal left singular vectors u1, u2 (principal vectors). So R' P_{V1⊥} C(s) lies on the
   ellipse E = {R'(sin θ1 cos s' u1 + sin θ2 sin s' u2)}, and σ1 v(w) lies in the δ1-neighbourhood
   E_{δ1} of E, δ1 := 3δ. Also |σ2 v(w)| ≤ |σ2 w| ≤ 2R'. **The tube is a thickened curve, so its
   projection is a thickened ellipse, not the filled box used in F5(b).**
3. *Covering.* Put a2 := R' sin θ2. If a2 ≤ δ1, E_{δ1} lies in a disc of radius 2δ1. If
   R' sin θ1 ≤ δ1, E_{δ1} lies in one rectangle (2a2 + 2δ1) × 4δ1. Otherwise split E into arcs of
   length ≤ ℓ0 and total turning ≤ Θ0 with ℓ0Θ0 = δ1, Θ0 = (2πδ1/perim)^{1/2} < π/2 (perim ≥ 4a2 > 4δ1);
   each arc is within ℓ0 sin Θ0 ≤ δ1 of its chord, and every rectangle has length ℓ'_k ≥ 4δ1. With
   ℓ0 = (δ1 · perim/2π)^{1/2} this gives N ≤ C(1 + (a2/δ1)^{1/2}) rectangles P_k of sides ℓ'_k × 4δ1,
   Σ ℓ'_k ≤ C(a2 + δ1).
4. *One rectangle.* Λ := P_{V⊥}(O) is an O_K-module, a lattice of rank 4. Let 𝔅_k = {σ1 x ∈ P_k,
   |σ2 x| ≤ 2R'} and D_k := 𝔅_k − 𝔅_k (symmetric, convex, a product in σ1 and σ2). Let
   λ1 ≤ … ≤ λ4 be the successive minima of Λ with respect to D_k. Multiplication by 1 + √2 gives
   λ2 ≤ 2.42 λ1 and λ4 ≤ 2.42 λ3. (For λ4: the Q-span of three Z-independent vectors attaining λ1..λ3 has
   odd dimension, so it is not stable under √2; hence some (1+√2)x_i, i ≤ 3, is independent of them.)
   - If λ3 ≤ 1: Henk and Minkowski II give #(Λ ∩ D_k) ≤ C vol(D_k)/covol(Λ).
   - If λ3 > 1: # ≤ C max(1, λ1^{−2}). A vector x1 attaining λ1 has |σ1 x1| ≤ 8λ1 ℓ'_k and
     |σ2 x1| ≤ 4λ1 R', and h(x1) ≥ 1/(C H_V) by step 5. So λ1^{−2} ≤ C H_V ℓ'_k R'.

   Hence #(Λ ∩ 𝔅_k) ≤ C(1 + H_V ℓ'_k R' + H_V ℓ'_k δ1 R'^2). The middle term is ≤ the last because
   δ1 R' = 3εR'^2 ≥ 1 (a0 < 2).
5. *Constants.*
   - covol(O) = covol(O ∩ V) covol(Λ). O_K^4 ∩ V ⊂ O ∩ V ⊂ ½(O_K^4 ∩ V), and the saturated module
     has covolume 8H_V, so covol(Λ) ≍ 1/H_V.
   - For 0 ≠ x = P_{V⊥} w ∈ Λ, take K-independent 𝔫1, 𝔫2 ∈ V⊥ ∩ O* with h(𝔫1)h(𝔫2) ≤ CH_V (as in
     F5(b)). Then ⟨𝔫_i, x⟩ = ⟨𝔫_i, w⟩ ∈ O_K, and one of them is nonzero, so
     1 ≤ h(𝔫_i) h(x) and h(x) ≥ 1/(CH_V).
6. Summing over k: #v(X_t) ≤ C(N + H_V R'^2 δ1 (a2 + δ1)) = C(1 + (a2/δ1)^{1/2} + 9 H_V R'^4 ε (sin θ2 + ε)). ∎

**New dichotomy constant.** Use A1 in F5(c) in place of (b). Put s2 := max(log sin θ2, −a0) and
X := log H_V. The fibration fails only if X + 4 − a0 + s2 > E*; note (sin θ2/ε)^{1/2} ≤ R'^{a0/2} < R'^{E*}.
On that region X + max(0, log sin θ2 + γ) ≥ G3 := max(0, γ + a0 + E* − 4):
- if log sin θ2 ≥ −a0, then X + log sin θ2 + γ > γ + a0 + E* − 4;
- if log sin θ2 < −a0, then X > E* + 2a0 − 4 ≥ G3, using γ ≤ a0.

So in Lemma F5, G = max(0, γ + E*/2 − 2) becomes G3. Compared with G: G3 − G = a0 + E*/2 − 2 > 0.
At a0 ≈ E* ≈ 1.63 and γ = 1.5, G3 = 0.76 against G = 0.32. This is the only change in
`../dioph/cert_annulus.py`.

## P1. Rank ≤ 1: O(1) normals

In Lemma F5 the set 𝒩 consists of *primitive* normals modulo units and sign. If at most one K-minimum
of O* with respect to 𝒦 is ≤ 1 (d ≤ 1), every vector of 𝒦 ∩ O* lies in one K-line Kv1. A primitive
n in that line generates Kv1 ∩ O* over O_K, so it is unique up to units, and #𝒩 ≤ 1.

The paper bounds #𝒩 by #(𝒦 ∩ O*) ≤ Ch, which gives the term `y` in nor(y). (The one-line case of the
sector count D6(a) never occurs when d = 2: 𝒦 ∩ V lies in the sector body, so μ′2 ≤ μ2 ≤ 1. The real gain of
P1 is y → 0 for d = 1.)

## P2, P3. Per-box projection along a rational subspace orthogonal to n_B

Let B be a box whose points span the hyperplane H_B = {⟨n_B, x⟩ = m_B}, with n_B ∈ O* primitive,
height h(n_B), sphere section Y_B, and patch lengths e_s, e_u (Lemma E2). By F3(b), the patch
Y_B ∩ B lies in ≤ 2 four-boxes with sides (Ce_u, e_s, 2e_u, 2e_u). This uses e_s^2/R' ≤ λ^2R' ≤ e_u,
which follows from 2b ≥ a and g ≤ 2b. Put ē := max(e_s, e_u).

**Lemma P3 (projection along a rational line).** If 𝔪 ∈ O_K^4 ∖ 0 is primitive and ⟨n_B, 𝔪⟩ = 0, then

    #(X_t ∩ B) ≤ 2^{O(t/log t)} · C · (1 + ē e_u R'^2 · h(𝔪)/h(n_B)).

*Proof.*
1. Let W := 𝔪⊥ and P_W the K-orthogonal projection. For w ∈ X_t ∩ B put u := P_W w ∈ Λ_W := P_W(O).
   Since n_B ∈ W, ⟨n_B, u⟩ = ⟨n_B, w⟩ = m_B. So u lies on the affine K-plane
   P_B := {u ∈ W : ⟨n_B, u⟩ = m_B} (dimension 2).
2. Given u, w = u + t𝔪 with nrd(w) = ⟨u,u⟩ + t^2⟨𝔪,𝔪⟩ = n_t, so there are ≤ 2 values of t, and
   #(X_t ∩ B) ≤ 2 #{u}.
3. σ1 u = P_{W1} σ1 w lies in the section, by the plane σ1(P_B), of the projection of the two
   four-boxes. That section is contained in (segment of length ≤ ē) + (disc of radius Ce_u), so it
   has area ≤ C ē e_u and diameter ≤ Cē. Also |σ2 u| ≤ 2R'.
4. The u's lie in u0 + L_B with L_B := Λ_W ∩ n_B⊥. This is a rank-2 O_K-module. To get its covolume:
   - covol(Λ_W) = covol(O)/covol(O ∩ K𝔪) ≍ 1/h(𝔪);
   - f := ⟨n_B, ·⟩ maps Λ_W onto f(O) = O_K (F1, since f∘P_W = f on O);
   - the quotient lattice in span(σ1 n_B, σ2 n_B) has covolume 2√2/h(n_B).

   So covol(L_B) ≍ h(n_B)/h(𝔪).
5. Count in D := {σ1 x ∈ A − A, |σ2 x| ≤ 4R'} with the K-minima μ1 ≤ μ2 of L_B:
   - If μ2 ≤ 1: all λ_i ≤ 2.42 and the count is ≤ C vol(D)/covol(L_B) ≤ C ē e_u R'^2 h(𝔪)/h(n_B).
   - Otherwise the differences of the u's lie in one K-line. Then the u's lie on an affine K-line ℓ,
     and w ∈ ℓ + K𝔪, an affine K-plane, i.e. a circle. Lemma E5 gives ≤ 2^{O(t/log t)} points (with the
     spanning remark of A1, step 1).

   Steps 3–5 are applied to each of the two four-boxes separately, which costs a factor 2. ∎

**Lemma P2 (projection along a rational plane).** If V ⊂ K^4 is a K-plane of height H_V and n_B ∈ V, then

    #(X_t ∩ B) ≤ 2^{O(t/log t)} · C · (1 + ē R' · H_V/h(n_B)).

*Proof.* Project along V⊥ onto V.
1. u := P_V w satisfies ⟨n_B, u⟩ = m_B, which is an affine K-line in V.
2. The fibre (u + V⊥) ∩ sphere is a circle, handled by E5.
3. The line lattice K_B := P_V(O) ∩ n_B⊥ = O_K κ has covolume ≍ h(n_B)/H_V, by the same exact
   sequence with covol(P_V O) ≍ 1/H_V. So h(κ) ≍ h(n_B)/H_V.
4. σ1 u lies in an interval of length ≤ Cē and |σ2 u| ≤ 2R'. By (H3) there are
   ≤ 1 + C ē R'/h(κ) values. As in P3, the fibres use E5 with the spanning remark, and the two four-boxes
   are treated separately. ∎

**Use in the LOW layers.** Fix a class interval [g_i, g_{i+1}] and a layer j, with heights in
[R'^{y_{j−1}}, R'^{y_j}). Let d be the number of K-minima ≤ 1 of O* with respect to 𝒦(h_j, ψ).
Every box of the layer has n_B ∈ 𝒦 ∩ O* (the body argument of F5).

| d | normals in | route N (normals + offsets + per-sphere, as in F7) | route P (all boxes, per-box bound) |
|---|---|---|---|
| 0, 1 | — | 0 + pm (P1) | — |
| 2 | V, height R'^X | max(0, 2y_j − max(X, G3′)) + pm, with G3′ = T − 4 + a0 + γ (A1 + D6(a) + P1) | E1 + max(0, ē + 1 + X − y_{j−1}) (P2) |
| 3 | W = 𝔪⊥, height R'^X | max(0, 3y_j − γ − X) + pm | E1 + max(0, ē + U0 + 2 + X − y_{j−1}) (P3) |
| 4 | — | max(0, 4y_j − 2γ) + pm | — |

In the table, ē is written as an exponent. pm := max over the tangent class and the crossing classes
of [offsets(y_j) + per-sphere(y_{j−1})], exactly as in `cert_dioph.py`.

Rank-3 route N uses (μ1μ2μ3)^2 ψ h^3 ≥ c H_W. The reason: v1∧v2∧v3 is a nonzero O_K-multiple of the
primitive 3-vector of W ∩ O_K^4, whose height is h(𝔪) = H_W (Hodge star). This is D5 with H_W kept.

Layer j is bounded by max over d of min(route N, route P), maximised over the unknown X (and θ2).
Both routes are monotone in X, so the max–min is attained at a breakpoint or at the crossing.
`cert_proj.py` evaluates it exactly. HIGH boxes and the rest of F7 are unchanged.

**Why this is not circular.** Route P counts each box once with its own hyperplane. The number of
boxes is R'^{E1} with κ = 0, and every point of X_t lies in some box on its hyperplane. The rank d is a
property of the lattice O* and the body 𝒦, not of the process. The adversary may choose it, and we
take the maximum over d.

## Numbers

| lemma set | α (float, pointwise) | certified |
|---|---|---|
| App. F (paper) | 2.435 | 17/7 = 2.4286 |
| + A1 | 2.458 | 27/11 (1/300), 86/35 = 2.4571 (1/1600, den 800) |
| + A1, P1, P2, P3 | ≈ 2.495 | **249/100** (1/800; also at den 800, η-den 200) |
| + all, κ = 0 is forced by E1 | ceiling 5/2 | — |

`breakdown2.py 2.49` shows the binding terms. HIGH and both route-P terms sit at E1 = b, and b is
forced up to (8 − 2a0)/3 by 2a + 3b > 8, so E1 < a0 requires a0 > 8/5.

## Second round (2026-09-30): the continuum limit is exactly 5/2

- **Certified α ≥ 499/200 = 2.495** (`cert_499_200_output.txt`:
  `cert_proj.py 800/499 800/499 9569497/5988000 3199501/1996000 800 400 12`, κ = 0, margin 1/4000).
  The same holds at η-den 500. α = 2.497 fails at g ≈ 0.786 with η-den 1000 and den 800.
- **The remaining loss is discretisation only.**
  - Near 5/2 the pointwise model (`model2f.py`, `probe25f.py`, `probe079.py`) has value(g) = E1 for every
    g in the continuum. At g ≈ 1.2–1.6 HIGH reaches E1 before route P/N of the LOW layers exceed it. At
    g ≈ 0.79 the crossing of HIGH (slope −1/3 in η) and LOW (slope ≈ 3) is also at E1.
  - The excess over T is ∝ the layer width: 0.0033 at Δη = 0.01 and 0.00066 at Δη = 0.002. It comes
    from route N using y_j and route P using y_{j−1}, plus the g-interval endpoints.
  - So every α < 5/2 has a certificate with fine enough grids. 5/2 is the exact limit of the lemma
    set, and E1 = b > (8 − 2a0)/3 is the only obstruction.
  - A single statement "α ≥ 5/2 − o(1)" would need an analytic, uniform-in-α version of the
    certificate. This is not done.
- **κ > 0 fails even with the extra structure of Minkowski's normal** (`model3.py`). The D1 normal
  ℓ_B of a box satisfies h(ℓ_B) ≥ 1 and |σ2 ℓ_B| ≤ T2, so |σ1 ℓ_B| ≥ 1/T2. Hence:
  - transverse tilt ≤ R'^{κ−2+a}, i.e. γ_eff = max(γ, 2 − a − κ);
  - height ≤ R'^{κ−2+max(2b,a)};
  - class g ≥ min(2 − a − κ, b).

  With these constraints the worst value is still 0.1–0.3 above a0 at α = 2.52 and 2.55. The binding
  pieces are medium and large tangent sheets (λR' × δ) on low-height sphere sections, at g ≈ 0.7–1.25.
  Neither the determinant nor the projections control them. This is the S^2-band problem again.
- **Degree-2 level 1** (ceiling 28/11). At α slightly above 5/2 the degree-2 boxes are nearly
  δ-cubes. Level 2 on S_B = Q ∩ S^3 is free at the scale of the box when the curvature of S_B is
  ≤ R'^{0.42} (at a0 = 1.58). A thin torus of radius ≤ δ around the core (curvature ≥ 1/δ) costs
  R'^{≈0.05} at level 2 alone, which exceeds the slack a0 − E1 ≈ 0.015.
  - Q comes from 13×13 minors, so its height is not small, and no Diophantine rarity of torus-like Q
    is available.
  - A "quadric dichotomy" would be needed: torus-like Q ⇒ Π is well approximated by a rational
    quadric ⇒ a global fibration by the values of Q, with near-Clifford tori as level sets. Not
    attempted.

## Where it stops, and the road to α = 3

1. **5/2: level 1 with one hyperplane per box.** The lemma set above now reaches it up to 0.01.
   Allowing κ > 0 (D1, R'^κ parallel hyperplanes per box) does *not* help in the model. HIGH pays R'^κ,
   and the best configuration at α = 2.5 fails by +0.015 (`search2.py`). The D1 hyperplanes of a box
   can be tilted into bands on great-ish spheres, which is exactly the hard case.
2. **28/11 ≈ 2.545: level 1 with auxiliary quadrics.** The degree-2 condition is 12a + 16b > 44.
   Recomputed here: the 14 functions mod the sphere have σ1-size λ^{16}ρ^{12}R'^{22} and σ2-size
   R'^{22}. Then E1 = b < a0 needs a0 > 11/7. Level 2 would then have to count points on a K-rational
   quartic del Pezzo section Q ∩ S^3 in a thin box. Plane sections are genus-≤1 quartic curves, where
   Lemma E5 no longer gives R'^{o(1)}. Thin tori around the core are admissible Q. Not attempted.
3. **8/3: any determinant method with free level 2** (d → ∞: a + b > 3).
4. **3: needs the sharp volume law at the critical scale.** α = 3 is equivalent to
   n(τ) ≤ 2^{τ/3+o(τ)} for τ ≤ 3L, uniformly in the frame. At τ = 3L this is the tube volume law
   R'^{4/3+o(1)} for tubes of radius R'^{−4/3}, whose expected content is N^{1/3} of the N ≍ R'^4
   points. For r = 2 it is equivalent (NOTES.md §7.1) to a uniform small-cap bound for ternary
   Z[√2]-spheres. That bound is the Z[√2] analogue of the Linnik/Bourgain–Rudnick small-cap problem,
   open over Q at these scales. There is no route to it by the elementary methods here.

   Conditionally, the D1 slicing at b → 0 gives E1 = 2 − a0/2 < a0 ⟺ a0 > 4/3. So *if every
   K-rational sphere section met the tube in R'^{o(1)}·max(1, expected) points*, α → 3. That
   hypothesis is a ternary (one dimension lower) form of Conjecture H, not a new route.
