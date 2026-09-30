# Referee report: Theorem U (α ≥ 5/2 − o(1)), `research/uniform52/NOTES.md`

**Overall verdict.**
- **Theorem U holds.** Its proof has three parts: the continuum model `cont_model.py`, the exact z3 check `verify_z3.py`, and the reduction Lemma U1. The model is a faithful pointwise limit of the reviewed assembly (App. F + A1/P1/P2/P3 + `cert_proj.py`). Every approximation goes in the safe direction. The z3 encoding of ¬C is correct and not vacuous. The HIGH/LOW reformulation and Lemma U1 are correct, but the write-up has gaps (all fixable).
- **The sharp constant E1 + (10/33)μ is confirmed exactly.** It holds even with the exact minimum over *all* admissible cells.
- **One stated claim is wrong: the "analytic check" of NOTES §3.** The value is *not* E1 + (10/33)μ on all of 2/5 ≤ g ≤ 8/5 − 7μ/3. On the lower part of the range the binding branch is rank 2, not rank 3. The supremum is unaffected.
- **Proposition S is valid as a conditional.** The sentence drawing its moral ("must show … sparse by a power") overstates what it proves.

| item | verdict |
|---|---|
| 1. Faithfulness of `cont_model.py` (patch_box, cellval, cost_upper, high, piece_total, routes) | **VALID** |
| 2. z3 encoding of ¬C (`verify_z3.py`), non-vacuity, threshold 10/33 | **VALID** |
| 3a. Reformulation "∃η*: HIGH ≤ T, LOW ≤ T on (0,η*]" ⟺ "no y with HIGH > T and LOW > T" | **VALID** |
| 3b. Lemma U1 (continuum ⇒ R'^{T+o(1)}), margins ϖ, R'^{o(1)} cases, Gibbs | **VALID WITH FIXABLE GAPS** (list in §3) |
| 4. Proposition S | **VALID**; the interpretive sentence needs rewording (§4) |
| 5. NOTES §3 "Where 10/33 comes from" | **INCORRECT AS STATED** (counterexample in §5); fixable, and the sup is right |

Every script is in this directory. Heavy runs were done on school-server (z3 5.1.0, numpy 1.26, `nice -n 10`).

---

## 1. Faithfulness — VALID

**Line-by-line comparison** against `../../dioph/cert_dioph.py`, `../../annulus/cert_proj.py`, App. F of `main.tex` and the table in `../../annulus/NOTES.md`:

| cont_model | source | check |
|---|---|---|
| `patch_box` | E2 (r ≤ R'/2): S0 = min(1−b, ½(2−a−min(g,a)), 1−g), U0 = min(1−a, 1−g); identical to `cert_dioph.patch_box` | same |
| `cellval` | E3/E4/D2: N = max(0, S0−XL, U0−XS, S0+U0−XL−XS, 2(U0−XS)), P = 1 + (XL+XS+min(XS, 2XL−R0) − y)/3, cost N + (P)₊; gated by `ok` | same, except XL ≤ R0 is non-strict (harmless, see below) |
| `cost_upper` | eq. (cost): min(patch (D3/F3(b)), cells) | the prover uses only 4 candidate cells, which can only weaken it (it is in fact exact, see below) |
| `high` | App. F HIGH: E1 + cost(S0, U0, 1−g, η*) | gl = gr |
| `piece_total` | F4(iii) offsets (y+lde+1 tangent, y+D+1 crossing); tube bound F2+F3(b)+F4(i) with S_T, U_T; box bound (S_T−(1−b))₊ + 2(U_T−(1−a))₊ + cost(S0b) with F4(ii) | same formulas as `low_pieces`/`cert_proj` with D_l = D_r, yp = yj |
| `routes` | NOTES table: d ≤ 1: pm (P1); d = 2: (2y − max(X, G3′))₊ + pm vs E1 + (ē+1+X−y)₊ (F5(a)+A1+P2); d = 3: (3y−γ−X)₊ + pm vs E1 + (ē+U0+2+X−y)₊ (P3); d = 4: (4y−2γ)₊ + pm | same; G3′ = T−4+a0+γ; ē = Sp = max(S0, U0) |

**Numerical agreement.**
- **Against the certificate code** (`compare_certproj.py`, 1500 random points, 3000 comparisons): I copied the per-layer LOW (v01, v2, v3, v4 with the certificate's own exact `maxmin`) and HIGH of `cert_proj.py` verbatim, with gl = gr = g, yp = yj = y and one crossing class D. They agree with my independent model to within 4·10⁻¹⁰ in both directions, over μ ∈ [0.001, 0.1] and T varied (T enters G3′).
- **Coverage.** This comparison includes the routes; the authors' `crosscheck.py` only checks the pieces and HIGH.

**Direction of every approximation.**
- *Prover side* (only weaker is allowed):
  - the cells are four explicit admissible candidates, each gated by `ok` (inadmissible gives BIG);
  - η* is existential;
  - route N vs route P, tube vs box, and patch vs cells are all minima.
- *Adversary side* (must be complete):
  - g ∈ (0, 2.01];
  - y > 0;
  - X ≥ 0;
  - ϑ = log sin θ2 is eliminated exactly by G3′. I re-derived this: min over the admissible ϑ of X + max(0, ϑ+γ) equals max(X, G3′).
  - D ∈ [lde, 1−2g];
  - piece ∈ {tangent, crossing};
  - rank d ∈ {≤1, 2, 3, 4}.
- The upper bounds H_V ≤ Cβh² and H_W ≤ Cψh³ are not used, which is conservative.

**Candidate cells.**
- **Admissibility.** Since U0 ≤ 1−g = R0 always, the candidates XL = (U0+R0)/2 and XL = XS = U0 are always admissible (non-strictly). XL = S0 and XL = min(S0, (U0+R0)/2) are admissible iff U0 ≤ S0. At U0 = S0 both coincide with XL = XS = U0, so `cost_upper` is continuous and piecewise linear, with slope in S0 between 0 and 1.
- **Exactness.** The candidates lose nothing:
  - `check_cellmin.py`: in 300 random cases they equal the exact minimum over all admissible cells (vertex enumeration of the line arrangement, validated against a dense 1801² grid);
  - the exact z3 optimum is also unchanged when the prover uses the full cell minimum (§2).

**Non-strict XL ≤ R0 versus strict XL < 1 − g_{i+1}.** Shift a continuum witness to (XL − 2ϖ, XS − 2ϖ):
- the shift preserves XS ≤ XL and 2XL − 1 ≤ XS;
- it gives XL ≤ 1 − g_{i+1} − ϖ;
- it raises N by at most 4ϖ and does not raise P.

So the non-strict inequality costs O(ϖ).

**Zero case, tiny spheres, crossing range.**
- *Zero case.* It is not part of C, but it is trivial. `zero_case.py` checks in z3 that E1 + cost(1−b, 1−a, 1−w, 0) = E1 for all μ ∈ (0, 1/10] and w ∈ [0, 1/100], via the whole-box cell with P = (w − 2/5 + μ)/3 < 0. It also checks that HIGH(g, 0) = E1 for g ≤ 2/5 − μ.
- *Tiny spheres* g > 2.01: the domain g ≤ 201/100 matches.
- *Crossing range.* The crossing range and the "no crossing for g ≥ a0" case match App. F.
- *Structural remark (not needed).* piece_total(D) is nonincreasing while y + D + 1 < 0 (offsets flat) and nondecreasing afterwards. So the maximum over D is at D = lde (= the tangent class) or at D = 1 − 2g. My first guess that "the deepest class always wins" was wrong; keeping D free, as the authors do, is correct.

## 2. z3 encoding — VALID

**Reproduction.** All five lines of `z3_outputs.txt` reproduce on school-server.

**Code audit of `verify_z3.py` + `cont_model.py`.**
- Every min, max, ite and And goes through the patched `Num`. A Python `min`/`max`/`if` on a z3 term raises an exception (tested in `mutate.py`), so it could not be silently evaluated.
- All divisions are Real divisions of Real terms.
- `c` and `mu0` are converted from `Fraction` to `z3.Q`.
- The `FR` override is dead code but harmless: `FR` is only used in `params()`, which `verify_z3` never calls.
- BIG = 10⁶ is only ever the value of an inadmissible candidate, and the patch term is always finite.
- The If-chains for max/min are correct.
- The existential over piece and rank is right, because route N is pm + (term) and route P does not depend on pm. So max over pieces of min(N, P) = min(N(max pm), P).
- **Truncations y ≤ 3 and X ≤ 6 are harmless.**
  - HIGH > T forces patch > 0, i.e. y < 3 + S0 + 2U0 ≤ 6/5 − 4μ/3.
  - For X ≥ 3y both route N's equal pm, which the disjunct `tot > T` already covers.
- The claim over μ ∈ (0, μ0] is linear in μ, so the check is quantifier-free linear real arithmetic.

**Independent encoding** (`indep.py` + `z3_indep.py`).
- **How it differs.** It is written from the paper and NOTES, not from `cont_model`. It has no bounds on y or X, and it can optionally use the exact full cell minimum. Each vertex of the arrangement is affine in (S0, U0, R0, y), because the line normals are constant, so the formula stays in QF-LRA.
- **Results** (`z3_indep_output.txt`). It gives the same sat/unsat pattern:
  - c = 1/2 unsat on (0, 1/10];
  - c = 10/33 unsat on (0, 1/300];
  - 303/1000 sat;
  - 0 sat.

  It also gives the same results in `full` mode.

**Exact threshold** (z3 Optimize of sup over (g, y, X, D, branch) of min(HIGH, LOW), with T self-consistent in G3′).
- The supremum is **exactly E1 + (10/33)μ** at μ = 1/1000, 1/300, 1/100, 1/30 and 1/10, in both `cand` and `full` mode. For example, at μ = 1/10 the sup is 86/55.
- Bisection agrees: at fixed μ = 1/10 and μ = 1/1000, c = 0.302999 is sat and c = 10/33 is unsat.
- **Stronger than claimed.** c = 10/33 is unsat on the whole range (0, 1/10]. So T = 8/5 − (4/11)μ holds for all μ ≤ 1/10, not only T = 8/5 − μ/6.

**Float brute force without z3** (`float_brute.py`, exact X-bisection, exact cells).

| μ | c found | note |
|---|---|---|
| 0.1 | 0.30303 | |
| 0.03 | 0.30303 | same in `cand` mode |
| 0.01 | 0.30300 | slightly low from the grid only |
| 0.003 | 0.30292 | slightly low from the grid only |

The float value is never above 10/33.

**Mutation tests** (`mutate.py`, at the tight c = 10/33 with μ ≤ 1/300).
- Adding +10⁻⁵ to any of HIGH, pm, the cell cost, N3, P3 or N2 turns unsat into sat. So does changing the route-P3 constant from 2 to 2.01, or lowering T by 1/10 inside G3′.
- P2 + 10⁻⁵ stays unsat, so P2 is not binding.
- Lowering N3 or P3 alone by 10⁻³ does not remove the counterexample at c = 0.303. The branch-restricted optimum (`branch_opt.py`) explains why:

| branch | exponent − E1, in units of μ |
|---|---|
| pm | 0 |
| v4 | 5/24 |
| rank 2 | 10/33 |
| rank 3 | 10/33 |

  **Both rank 2 and rank 3 bind.** The encoding reacts to every binding formula.

## 3. Reformulation and Lemma U1

**3a. Reformulation — VALID.**
- HIGH_c(g, ·) is continuous, piecewise linear and nonincreasing: every candidate term decreases in η, and admissibility does not depend on η.
- HIGH → E1 < T, so η* = min{η : HIGH ≤ T} exists and is attained.
- LOW_c(g, ·) is continuous: it is a supremum over the compact set of D and over the effectively compact X ∈ [0, 3y] of jointly continuous piecewise-linear functions. (Lower semicontinuity would already suffice.) So LOW ≤ T on (0, η*) implies LOW(η*) ≤ T.
- The converse direction is immediate.

**3b. Lemma U1 — VALID WITH FIXABLE GAPS.**

I checked that nothing in App. E/F or the annulus notes needs a *fixed* margin. Every use of a margin only needs R'^{margin} to beat an absolute constant:
- E1 with 2a0 + 3b_ϖ = 8 + 3ϖ;
- E3's X_L < 1 − g, i.e. 2x_L ≤ r/C;
- the zero interval [0, ϖ] containing all r > R'/2 (needs ϖ log R' ≥ 1);
- tiny sections (fixed 0.01, unchanged);
- A1's (sin θ2/ε)^{1/2} ≤ R'^{a0/2} < R'^{E*} (fixed);
- the tangent/crossing threshold C′δ′ (constants);
- the o(1) terms in F5(c).

The counting is also fine:
- The number of cases is O(ϖ⁻³) = (log log R')³ = R'^{o(1)}. Each case carries the same absolute constants and one E5 factor.
- F5 is applied to O(ϖ⁻²) pairs (g_i, y_j), each with its own V or W. Its first alternative is a global bound, so a union over the pairs is fine.
- **Gibbs step.** With λ = 1/α, the difference λ − T/4 = (a0 − T)/4 > 0 is fixed and η(τ) = o(τ) is uniform in the frame. So Z_λ = O_α(1): A is even O_α(1), and certainly o(L).

**Empirical convergence.** I ran the reviewed finite certificate `cert_proj.py` itself at μ = 3/100 (a0 = 163/100, b = E1 + 1/1000); output in `cert_convergence_output.txt`:
- T = E1 + μ/2 with grids 1/400 × 1/200: **CERTIFIED**;
- T = E1 + 0.4μ with grids 1/800 × 1/400: **CERTIFIED**;
- T = E1 + 0.29μ (below 10/33): **FAILS** at g ≈ 0.76, as it must for any grid.

This is exactly the behaviour Lemma U1 predicts.

**Gaps to fix in the write-up.**
1. **Missing hypotheses.** State the side conditions in U1: a0 + a ≥ 2b, a ≥ b, 2b ≥ a (for D3 and P2/P3; here 2b − a = 8/5 − 7μ/3 > 0), and T > E1 (needed for η* to exist and for the boxes spanning at most a plane). All hold for 0 < μ ≤ 1/10, but U1 as stated is for general T < a0.
2. **Step 3 ("each quantity exceeds its continuum value by ≤ Kϖ").** The continuum value must be taken at g_i. Spell out the three places where the certificate uses g_{i+1}:
   - the sagitta term and the constraint XL < 1 − g: use the witness shift of §1 (cost ≤ 4ϖ);
   - the lower end lde(g_{i+1}) < lde(g_i) of the crossing range: crossing classes with D ≤ lde(g_i) have exactly the tangent class's S_T and S0b and fewer offsets, so they are dominated by the tangent piece at g_i;
   - the no-crossing case 1 − 2g_i < lde(g_{i+1}): Δ ≥ C′δ′ ≥ R'^{lde(g_{i+1})} > R'^{1−2g_i} ≥ Δ is impossible.

   Also say explicitly that the Lipschitz statement needs `cost_upper` to be continuous across the flips of `ok`, which holds by the coincidence noted in §1.
3. **Step 4.** The last layer [y_{j*−1}, y_{j*}] contains η*. It closes because the discretized value is ≤ LOW_c(y_{j*−1}) + Kϖ, with y_{j*−1} < η*. State this.
4. **Zero case.** C does not contain the zero case. Add the one-line check of §1: the whole-box cell has P = (w − 2/5 + μ)/3 < 0.
5. **App. E's global sentence.** "Every inequality between exponents below is strict with a fixed margin" must be replaced, for U1, by "margins ≥ ϖ with R'^{ϖ} → ∞". Give the list of uses above.

## 4. Proposition S — VALID

Checked against Def. `def:model`, Def. `def:family`, the CW output model and Lemma `lem:loc`.
- **Legal process.** It reads each round once. The round-1 words W_u are emitted while reading round 1. The snapshot at the single cut is the vector ρ ∈ [0, γ)^m. Everything after the cut is a function of (σ, a2): V depends on ρ_j, a2_j and the public constant G, which belongs to the program, not the snapshot. The process is deterministic, so there is no tape to store, and it is zero-error on every valid stream, not only under 𝒟.
- **Snapshot bound.** S ≤ ⌈m log2 γ⌉ + O(log m) + O(1) = m(log2 γ + O(1)). The O(log(mQr)) of Prop. `prop:block` is not needed here: the cut position is fixed and r = 2.
- **Correctness.** d_proj is left- and right-invariant and satisfies the triangle inequality. R_z(ΔQ) = −I is a phase. So d(W_u V, R_z(Δx)) ≤ ε/2 + ε/2.
- **Rate.** T_pre = T_1 ≤ mτ_ε, and the bits forgone are m(log2 K − log2 γ − O(1)) = m(1 − o(1))L. So any bound of Theorem U's form with α > 5/2 fails for large L along the sequence. Combined with Theorem U, the rate is exactly 5/2 under the hypothesis.

**Fix the moral drawn from it.** "Any proof of a rate above 5/2 must show that tubes at T-count 5L/2 are sparse by a power" does not follow from Prop S.
- Prop S only forces that no tube at T-count (5/2 + o(1))L has maximal gap ≤ Q·2^{−(1−o(1))L}. A tube with 2^{(1−o(1))L} words concentrated on part of the circle does not satisfy S's hypothesis.
- Sparsity by a power is *sufficient* to block this construction, not *necessary*.
- A randomised variant makes the true condition stronger than "no dense tube". Store one shared random rotation offset s in the snapshot (log2 Q bits, amortised over m coordinates) and code ρ_j with an Elias-δ code. The maximal gap is then replaced by the length-weighted mean Σ_i (γ_i/Q) log2 γ_i + O(log L) per coordinate. By Jensen, this mean is ≥ log2(Q/|U|).
- So what a proof of a rate above 5/2 must exclude is tubes with small *weighted mean log-gap*. That condition lies between "no dense tube" and "sparse by a power".
- The heuristic paragraph on quadrics is labelled as a heuristic and was not reviewed.

## 5. NOTES §3 (analytic check) — INCORRECT AS STATED

**The claim.** §3 says that for every g ∈ [2/5, 8/5 − 7μ/3] the value is E1 + (10/33)μ, and that it comes from the rank-3 routes.

**The rank-3 computation needs its flat point to be admissible.** It assumes route P3 is flat up to X0 = y − 4/5 + μ/3. This requires X0 ≥ 0, i.e. y ≥ 4/5 − μ/3, i.e. g ≳ 6/5.

**Counterexample** (`fixed_g_opt_output.txt`, exact z3 at fixed g). The table gives exponent − E1 in units of μ.

| μ | g | rank 2 | rank 3 |
|---|---|---|---|
| 1/300 | 1/2 | 0 | 0 |
| 1/300 | 3/5 | 0 | 0 |
| 1/300 | 4/5 | 10/33 | 0 |
| 1/300 | 1 | 10/33 | 0 |
| 1/300 | 6/5 to 3/2 | 10/33 | 10/33 |
| 1/10 | 1/2 | 1/15 | 0 |
| 1/10 | 3/5 | 1/6 | 0 |

The float profile gives the same picture: at μ = 0.1 the value equals E1 + (10/33)μ only on g ∈ [0.74, 1.37], and at μ = 0.01 only on [0.79, 1.58].

**Rank 2 gives 10/33 through a different mechanism.** Put y = g − 2/5 + δ with 0 < δ < μ(1+c) and T = E1 + cμ, so G3′ = g − 4/5 + μ/3 + cμ lies above the start of the rise of P2. Then:
- route N2 is flat at E1 + 3δ + μ/3 − cμ;
- the balance with HIGH − E1 = (μ − δ)/3 gives δ = 3cμ/10;
- self-consistency, value = cμ, gives c(1 + 1/10) = 1/3, i.e. **c = 10/33 again**.

So 10/33 arises from two different balances: rank 3 for g ≳ 6/5, and rank 2 (through the T in G3′) on the range below.

**Where the value is attained.** The rank-2 mechanism needs G3′ ≥ 0, i.e. g ≥ 4/5 − μ/3 − (10/33)μ = 4/5 − 7μ/11. The upper end is the whole-box patch condition g ≤ 8/5 − 7μ/3. So the value is E1 + (10/33)μ exactly for

  g ∈ [4/5 − 7μ/11, 8/5 − 7μ/3],

and below that elsewhere. The float scan agrees: [0.738, 1.366] against the predicted [0.736, 1.367] at μ = 0.1, and [0.794, 1.576] against [0.7936, 1.5767] at μ = 0.01.

§3 should be corrected to say this. The supremum, the sharpness claim and Theorem U are unaffected.

**Minor points.**
- **Sharp-constant claim.** "The continuum exponent of the lemma set is exactly E1 + (10/33)μ" holds on all of μ ≤ 1/10 (tested at five values), and also with the exact cell minimum. It is the exponent of this assembly (a = a0, b = E1, κ = 0), not of the determinant method at large.
- **Stronger tube bound.** `z3_outputs.txt` could record that c = 10/33 already holds on (0, 1/10]. That gives the stronger tube bound R'^{8/5 − 4μ/11}.

## Files (this directory)

| script | output | what |
|---|---|---|
| `indep.py` | — | referee's continuum model (backend-agnostic; `cand` or exact `full` cells) |
| `z3_indep.py` | `z3_indep_output.txt` | independent z3 ¬C (no y/X truncation), fixed-μ checks, exact Optimize |
| `branch_opt.py`, `fixed_g_opt.py` | `branch_opt_output.txt`, `fixed_g_opt_output.txt` | exponent per branch; value at fixed g |
| `float_brute.py` | `float_brute_output.txt` | numpy brute force (no z3) |
| `check_cellmin.py` | `check_cellmin_output.txt` | vertex enumeration = dense grid; candidate = exact min |
| `compare_certproj.py` | `compare_certproj_output.txt` | per-layer agreement with the verbatim `cert_proj.py` code |
| `mutate.py` | `mutate_output.txt` | mutation tests of the authors' encoding |
| `zero_case.py` | `zero_case_output.txt` | zero case and no-split region |
| (runs of `../../annulus/cert_proj.py`) | `cert_convergence_output.txt` | finite certificates converge (c = 1/2, 0.4 pass; 0.29 fails) |
