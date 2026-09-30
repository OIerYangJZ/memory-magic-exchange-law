# Referee report: A1, P1, P2, P3 and the 249/100 certificate (research/annulus/NOTES.md)

**Overall verdict.** Everything I checked holds: the four lemmas, the new dichotomy constant G3, the LOW-layer assembly and `cert_proj.py`. With them, **α ≥ 249/100 follows from App. E/F** (which the earlier review accepted), and α ≥ 27/11 follows from A1 alone. I found no counterexample, no broken step and no case the argument misses. Each lemma is **VALID**. The remaining gaps are all in the write-up, and all can be fixed (listed below). The margins are very small: a0 − T = 1/800, and the first threshold the certificate accepts gives 1.605111 against T = 1.605176. Every check is exact, so this is sound, but the result is purely asymptotic.

| item | verdict |
|---|---|
| A1 (annular fibration) | VALID, wording gaps fixable |
| G3 = max(0, γ + a0 + E* − 4) (new F5(c)) | VALID |
| P1 (rank ≤ 1: #𝒩 ≤ 1) | VALID. The half about the one-line case of F5(a) changes nothing (see below). |
| P2 (projection along V⊥) | VALID, wording gaps fixable |
| P3 (projection along 𝔪) | VALID, wording gaps fixable |
| body claim "n_B ∈ 𝒦 ∩ O*, hence n_B ∈ V (d=2) or W (d=3)" | VALID |
| dichotomy / route N vs route P / max over d, X, θ2 | VALID, no circularity, no case missing |
| `cert_proj.py` (and `../dioph/cert_annulus.py`) | VALID. It implements the NOTES table exactly, with exact arithmetic and conservative endpoints. |

---

## A1 — VALID (fixable gaps)

I checked every step.
- **Fibres.** The form is anisotropic, so K⁴ = V ⊕ V⊥. P_{V⊥} is K-linear and commutes with σ_i, so σ1 P_{V⊥} = P_{V1⊥} σ1.
- **Singular values.** P_{V1⊥}|_Π has Gram matrix diag(sin²θ_i) in the principal basis. So its singular values are sin θ1 ≤ sin θ2 and its left singular vectors are orthonormal. The ellipse claim follows: the tube point is within R'(1 − cos φ) + R' sin φ ≤ 3δ of E.
- **Covering.** The three cases are correct. In the third case, N ≤ perim/ℓ0 + 2π/Θ0 + 1 ≍ (a2/δ1)^{1/2}, and Σℓ'_k ≤ perim + 4δ1N ≤ C(a2 + δ1).
- **Step 4.** The rest of step 4 checks out: (1+√2) scales σ1 by 2.414 and σ2 by 0.414, and the Henk + Minkowski II split at λ3 ≤ 1 versus λ3 > 1 is right.
- **Step 5.** O_K⁴∩V ⊂ O∩V ⊂ ½(O_K⁴∩V) holds. So does ⟨𝔫_i, P_{V⊥}w⟩ = ⟨𝔫_i, w⟩ ∈ O_K, and one of the two values is nonzero because the form is nondegenerate on V⊥. Hence h(x) ≥ 1/h(𝔫_i) ≥ 1/(CH_V).
- **Sum.** Σ_k C(1 + H_V ℓ'_k δ1 R'²) gives the stated bound.

Fixes:
1. **Applying E5 to a fibre.** E5 is stated for a plane *spanned by points of X_t*. A fibre w + V is a K-rational plane but is not a priori spanned by such points. Add: "a fibre with ≥ 3 points is spanned by them (three distinct points of a circle are affinely independent); fibres with ≤ 2 points are trivial". The same sentence is needed in P2 step 2 and P3 steps 2 and 5. Paper F5(b) uses the same step silently.
2. **Rectangle sizes and chord approximation.** Take every rectangle with ℓ'_k ≥ 4δ1 (true by construction: chord + 4δ1, or 4δ1 in the disc case). This is what |σ1 x1| ≤ 8λ1ℓ'_k needs. Also record Θ0 = (2πδ1/perim)^{1/2} < 1.26 < π/2, which follows from perim ≥ 4a2 > 4δ1, so that "within ℓ0 sin Θ0 ≤ δ1 of the chord" is justified.
3. **Why λ4 ≤ 2.42λ3.** The parenthetical reason given is incomplete. The correct reason: the Q-span of three Z-independent vectors has odd Q-dimension, so it is not √2-stable. Hence some (1+√2)x_i with i ≤ 3 is independent of x1, x2, x3. The claim itself is true (checked on 74 instances in `check_count.py`).
4. **Cosmetic.**
   - |σ2 w| ≤ R': the σ2-radius is ≤ the σ1-radius. NOTES writes 2R'; it is only a constant.
   - δ1R'²(a2 + δ1) = 3εR'⁴(sin θ2 + 3ε).
   - The docstring of `cert_proj.py` omits the (sin θ2/ε)^{1/2} term. This is harmless because it is ≤ R'^{a0/2} < R'^{E*}.

**New F5(c) (G3).** I re-derived the adversary's optimum. For fixed X, with the failure region X + 4 − a0 + max(ϑ, −a0) > E*:

    min over ϑ of X + max(0, ϑ+γ)  =  max(X, G3′),   G3′ = E* + a0 + γ − 4.

The minimum over X is G3 = max(0, G3′). This uses γ ≤ a0 in the branch ϑ < −a0. `test_maxmin.py` cross-checks it against a direct maximisation over θ: 0 violations in 2000 random cases.

## P1 — VALID

- **The set 𝒩.** 𝒩 consists of primitive normals modulo units (paper, before F1, and the statement of F5). For d ≤ 1, 𝒦 ∩ O* ⊂ Kv1, and Kv1 ∩ O* = O_K n0. So 𝒩 has at most one element.
- **The assembly is consistent with this.** The paper counts offsets per primitive normal (F4(iii), m_H ∈ O_K), and each hyperplane has one primitive normal up to units. So #𝒩 × offsets × per-sphere with #𝒩 ≤ 1 is legitimate. The paper's #𝒩 ≤ #(𝒦∩O*) was just wasteful: it also counted the O_K-multiples of n0.
- **Remark (harmless).** The "one-line case of D6(a)" can never occur when d = 2. We have 𝒦 ∩ V ⊂ 𝒦′ (the body of F5(a)), so μ′_2 ≤ μ_2 ≤ 1. The h-term of F5(a) was therefore vacuous already, and the real gain of P1 is the d = 1 term (y → 0).

## P2 — VALID (fixable gaps)

- **Steps 1–2.** ⟨n_B, P_V w⟩ = ⟨n_B, w⟩ because n_B ∈ V. The fibres u + V⊥ lie in H_B, since V⊥ ⊂ n_B⊥.
- **Step 3.**
  - f = ⟨n_B, ·⟩ maps P_V(O) onto f(O) = O_K (F1). I checked this exactly.
  - Exact identity: covol(K_B)² = covol(P_V O)² h(n_B)²/8.
  - covol(P_V O)·H_V lies in [1, 4] (bounds [½, 8]), and covol(K_B)·H_V/h(n_B) lies in [0.35, 1.41] (`check_covol.py`, two seeds, 50 random planes).
- **Step 4.** P_{V1}(Q_k) has diameter ≤ C(e_s + e_u), and (H3) applied to differences gives the count. Monte Carlo gives a section length ≤ 2.75ē (`check_section.py`).

Fixes: (i) treat the two four-boxes of F3(b) separately (factor 2); (ii) add the E5 spanning remark.

## P3 — VALID (fixable gaps)

- **Step 2.** nrd(u + t𝔪) = ⟨u,u⟩ + t²⟨𝔪,𝔪⟩, so there are ≤ 2 values of t.
- **Step 3 (area).** Pick f2 ∈ σ1(P_B) orthogonal to the projected segment. The section has f2-width ≤ 2Ce_u and diameter ≤ Cē, so its area is ≤ Cēe_u. Monte Carlo gives area ≤ 4.2 ēe_u and diameter ≤ 3.0ē.
- **Step 4 (covolumes).** All checked exactly, over 50 lines 𝔪:
  - covol(O ∩ K𝔪)·covol(Λ_W) = 4;
  - covol(Λ_W)·h(𝔪) ∈ [√2, 4√2] (observed [1.41, 5.66]);
  - f(Λ_W) = O_K;
  - covol(L_B)² = covol(Λ_W)² h(n_B)²/8, hence covol(L_B)·h(𝔪)/h(n_B) ∈ [½, 2].
- **Step 5.** The K-minima case split is right.
  - If μ2 ≤ 1, all λ_i ≤ 2.42. Observed #/(vol/covol) ≤ 2.25.
  - If μ2 > 1, all differences lie in one K-line. I verified this *exactly* on the 45 such instances: the points lie on an affine K-line, and ℓ + K𝔪 is a circle.

Fixes: the same as for P2.

**Rank-3 route N (H_W kept).** v1∧v2∧v3 = ξ·(s1∧s2∧s3) with ξ ∈ O_K∖0, and s1∧s2∧s3 = ±unit·⋆𝔪. So H_W = h(𝔪). I checked exactly that covol(W ∩ O_K⁴) = 16√2·h(𝔪), and that h(v1∧v2∧v3)²/h(𝔪)² is a nonzero square integer.

## Body claim: n_B ∈ V (d=2), n_B ∈ W (d=3) — VALID

- **n_B ∈ 𝒦.** A LOW box of layer j in [g_i, g_{i+1}] has class ≥ g_i, meets the tube and has h(n_B) < R'^{y_j}. Then |P_{Π⊥}ν| ≤ 2ε + 2r/R' ≤ Cψ with ψ = R'^{−min(g_i, a0)}. Its balanced representative has |σ_i n| ≤ 2√h. So n_B ∈ 𝒦(h_j, ψ) ∩ O*.
- **Consequence.** If d = 2, any vector outside V would force μ3 ≤ 1; if d = 3, any vector outside W would force μ4 ≤ 1. The V (or W, with its X) used in route P is therefore the same one as in route N, for every box of the layer.

## Dichotomy and assembly logic — VALID

- **A1 is used once per case.** It is applied to the plane V_{ij} for the finitely many pairs (i, j) with d = 2. If any one succeeds, #X_t ≤ R'^{T+o(1)} globally. Otherwise the failure constraint holds in every such pair, and only route N uses it. Route P needs no constraint.
- **Both routes bound the same set.** Each bounds the points of X_t in the LOW boxes of layer j:
  - Route N counts sphere sections globally.
  - Route P counts boxes, R'^{E1} with κ = 0, times the per-box bound.
  - So the minimum of the two is legitimate.
- **Rank and adversary.** d is a property of (O*, 𝒦) and is maximised over. X and θ2 are shared by the two routes within a layer, and the adversary maximises over them.
- **Upper bounds on X not used.** H_V ≤ Cβh² and H_W ≤ Cψh³ would cap X, but they are not used. This is conservative.
- **Cases covered.**
  - boxes spanning at most a plane (R'^{b+o(1)});
  - g > 2.01 (coplanar);
  - the zero case;
  - the no-split intervals;
  - HIGH (unchanged);
  - LOW layers j ≤ j* for d = 0…4, with the tangent and crossing classes inside pm;
  - heights ≥ 1, so the layers start at y = 0.
- **No circularity.** No quantity depends on the bound being proved.

## `cert_proj.py` audit — VALID

- **It implements the table exactly.**
  - v01 = pm; v4 = (4y − 2γ)₊ + pm.
  - v3 = maxmin((3y−γ−X)₊ + pm, E1 + (ē + U0 + 2 + X − y_{j−1})₊).
  - v2 = maxmin((2y − max(X, G3′))₊ + pm, E1 + (ē + 1 + X − y_{j−1})₊).
  - ē = max(S0, U0) is the E2 value at g_i. The docstring says "e_s", but the code is right.
- **Endpoints are conservative.**
  - ē, U0 and γ (hence G3′) at g_i;
  - route N at y_j, route P at y_{j−1};
  - pm, HIGH and the cells exactly as in `cert_dioph.py`.
  - Every quantity is monotone in the direction used.
- **The max–min over X is exact.**
  - f is nonincreasing and continuous; g is nondecreasing, continuous and → ∞.
  - The candidates are 0, the breakpoints (A, −B, G3′) and every pairwise intersection of the affine pieces.
  - The X → ∞ limit (value pm) is covered.
  - Verbatim copies of these lines match a brute-force X-grid in 3000/3000 random cases for each rank, with no under-estimate (`test_maxmin.py`).
- **Arithmetic is exact.**
  - All arithmetic is in Fractions. Floats only propose cell witnesses, which are re-checked exactly.
  - One subtlety: the `Cost` cache ignores `need`, so a later call may receive a patch-only value. That value is still an upper bound, so the check stays conservative.
- **Reproduced.** `cert_proj.py 400/249 400/249 1192747/747000 319751/199200 400 100 12` → CERTIFIED (15 s), with output identical to `cert_249_100_output.txt`. `cert_annulus.py 44/27 …` → CERTIFIED 27/11.
- **The certificate is not vacuous.** The same script **fails** at α = 2.497 (interval [301/400, 151/200]).
- **Cross-check against the model.** At 615 sampled (interval, layer) pairs, the certificate's LOW and HIGH values are never below my pointwise model (difference ≥ −9e−16), as they must be (`compare_cert_model.py`).
- **The "worst exponent" printed is not the optimum.** It is the first threshold y* that succeeds, as in `cert_dioph.py`. At (1.5975, y* = 1.16) the model's optimum there is E1 = 1.5967.

## Independent float model (`indep_low.py`, written from NOTES + App. F)

How it differs from the authors' model:
- the cell cost is minimised exactly, by vertex enumeration;
- everything is pointwise in g;
- the rank-2 route N is maximised directly over (X, θ);
- the max–min uses bisection.

Results at the certificate's parameters (a0 = a = 400/249, b = 1192747/747000, T = 319751/199200):
- dy = 0.002, dg = 0.0025, 24 D-points: **sup_g = 1.59836** at g ≈ 1.20, y* ≈ 0.80. Margins: T − sup = +0.0068, a0 − sup = +0.0081.
- Layers of width 0.01, as in the certificate: sup = 1.60119. Margin: a0 − sup = +0.0052.
- The authors' `model2.py` gives 1.60000.

The binding terms sit at E1 = b plus a little: route P3/P2 against route N at g ≈ 1.2–1.6.

**Float limit.** With a0 = a = 4/α, b = (8−2a0)/3 + 1e−4, T = a0 − 1e−4, dy = 0.002:
- 2.498 passes (+0.0010);
- 2.499 passes (+0.0002);
- 2.4995 fails (−0.0002).

So the pointwise limit is **≈ 2.499**, essentially the 5/2 ceiling set by E1 = b > (8−2a0)/3. NOTES' "≈ 2.495" is the coarser-grid value. dy = 0.005 gives 2.498.

## Exact / numerical checks (all in this directory)

| script | checks |
|---|---|
| `klat.py` | exact Q(√2)-lattice toolkit (Fractions): O, O*, intersections, projections, HNF, exact LLL, enumeration |
| `check_covol.py` → `check_covol_output.txt` | covol O = 4, O* = 2O; covol(S) = 8H_V (vs h(v1∧v2)/index); covol(O∩V)covol(P_{V⊥}O) = 4; height(V⊥) = height(V); P2/P3 exact sequences; Hodge star; f(·) = O_K; min h on P_{V⊥}O ≥ 0.5/H_V |
| `check_ellipse.py` → `…_output.txt` | singular values = sin θ_i (1e−12); dist to E ≤ 0.997δ (≤ 3δ claimed); covering N/(1+√(a2/δ1)) ≤ 9.2, Σℓ′/(a2+δ1) ≤ 16, arc–chord deviation ≤ 0.125δ1 |
| `check_count.py` → `…_output.txt` | λ2 ≤ 2.4143λ1 and λ4 ≤ 2.4143λ3 (74 instances); #(Λ∩D) ≤ 1.05(1 + H_Vℓ′R′ + vol/covol); P3 lattice case #/(vol/covol) ≤ 2.25; one-K-line case exact |
| `check_section.py` → `…_output.txt` | P3 section area ≤ 4.2ēe_u, diameter ≤ 3.0ē; P2 segment ≤ 2.75ē |
| `test_maxmin.py`, `compare_cert_model.py` | certificate max–min exactness; certificate ≥ model |
| `indep_low.py` → `indep_low_249_output.txt`, `indep_alpha_scan.txt` | independent pointwise model, float limit |

## Caveats (not errors)

- **Everything is asymptotic.** The margins are ~1e−3 in the exponent (a0 − T = 1/800), and the first-admissible slack is 6.5e−5. Every exponent inequality the certificate checks is exact, so the argument is sound. But the implied R'₀ is astronomically large.
- **Presentation at paper standard.**
  - F5 must be restated to carry X = log H_V (d = 2) and X = log H_W (d = 3) into the assembly.
  - The assembly paragraph must add route P, with its endpoint conventions (ē, U0, γ at g_i; route P at y_{j−1}).
