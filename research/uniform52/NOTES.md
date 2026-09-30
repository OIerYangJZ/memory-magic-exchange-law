# α ≥ 5/2 − o(1), and why 5/2 is where this method stops (research notes, 2026-09-30)

**Review:** `review/REVIEW.md` (independent). Verdicts:
- faithfulness of `cont_model.py`, the z3 encoding, the reformulation and Proposition S: VALID;
- Lemma U1: VALID WITH FIXABLE GAPS;
- the old §3 analytic check: incorrect as stated, with the supremum right.

All fixes are applied below.

Builds on `../annulus/NOTES.md`: the lemmas A1, P1, P2 and P3, and their certificate `cert_proj.py`
(certified 249/100 and 499/200, independently reviewed).
Here the finite certificates are replaced by one exact statement that covers every α < 5/2.

## Results

**Theorem U (rate five halves, asymptotic).** For every α < 5/2 there is L0(α) such that, for every
process and every round t < r,
T̄_t ≥ α[(1−2δ)m log2 K − 2m h2(δ) − S − m A],  with A = o(L).
Equivalently, every bit forgone costs at least 5/2 − o(1) committed T gates. This holds
unconditionally, in the full CW model, for all frames.

**Tube bound behind it.** Put a0 = 8/5 + μ with 0 < μ ≤ 1/10. For every frame,

    #X_t ≤ R'^{8/5 − (4/11)μ + o(1)}   (tube radius ε ≤ C R'^{−a0}),

which is below R'^{a0} by R'^{(15/11)μ}. This is c = 10/33, which is unsat on all of (0, 1/10]. For
a0 > 17/10, use monotonicity in ε. The Gibbs term is A = O_α(1): the fixed gap λ − T/4 absorbs the
o(τ) terms.

**Sharp constant of the lemma set near the limit.** The continuum exponent of this assembly
(a = a0, b = E1, κ = 0, candidate cells equal to the exact cell minimum) is exactly
E1 + (10/33)μ = 8/5 − (4/11)(a0 − 8/5). This holds on all of 0 < μ ≤ 1/10: z3 Optimize at five values
of μ, and float brute force. It is the exponent of this assembly, not of the determinant method at
large. It is < a0 iff a0 > 8/5, so 5/2 is precisely the limit of the
method: approached, not attained.

**Proposition S (dense tubes would make 5/2 sharp).** See the last section. Suppose that for
arbitrarily small ε some tube at T-count (5/2 + o(1))L contains words at grid angles with maximal gap
≤ Q·2^{−(1−o(1))L}. Then an explicit r = 2 process attains exchange rate 5/2, so the true rate would
be exactly 5/2.

This only rules out *dense* tubes. A randomised variant replaces the maximal gap by a length-weighted
mean log-gap; see §4. What a proof beyond 5/2 must exclude lies between "no dense tube" and "sparse
by a power".

## 1. Continuum reduction (Lemma U1)

The certificates of App. F / `cert_proj.py` check finitely many cases: g-intervals of width 1/den,
height layers of width 1/eden, crossing classes of width ~1/nD, and κ = 0 with 2a+3b > 8 by a fixed
margin. Every exponent in them is piecewise linear, with slopes bounded by an absolute constant, in
(a0, a, b, T, g, y, D, X) and in the cell exponents (X_L, X_S).

**Lemma U1.** Suppose:
- E1 := (8 − 2a0)/3 < T < a0, with a = a0 and b = E1;
- the side conditions of the assembly hold: a ≥ b, a0 + a ≥ 2b, and 2b ≥ a. The last one is needed for
  D3 and for P2/P3; here 2b − a = 8/5 − 7μ/3 > 0.

All of these hold for 0 < μ ≤ 1/10. If moreover the continuum claim C(a0, T) of §2 holds, then
#X_t ≤ R'^{T+o(1)}.

*Proof.*
1. Take ϖ = ϖ(R') → 0 with ϖ log R' → ∞, e.g. ϖ = 1/log log R'.
2. Run the assembly of App. F (with A1, P1–P3) with:
   - b_ϖ = E1 + ϖ, so 2a0 + 3b_ϖ = 8 + 3ϖ, and Lemma E1 applies as soon as R'^{3ϖ} exceeds its absolute constant;
   - g-intervals, height layers and crossing classes all of width ϖ;
   - every strict exponent inequality that the lemmas need "with a fixed margin" imposed with margin ϖ.
     Each such inequality only has to beat an absolute constant, and R'^{ϖ} → ∞ does that.
   App. E's sentence "every inequality between exponents is strict with a fixed margin" is replaced by
   "margins ≥ ϖ with R'^{ϖ} → ∞". The uses are:
   - E1 with 2a0 + 3b_ϖ = 8 + 3ϖ;
   - E3's X_L < 1 − g, i.e. 2x_L ≤ r/C;
   - the zero interval [0, ϖ] containing all r > R'/2 (needs ϖ log R' ≥ 1);
   - the tiny sections (fixed 0.01);
   - A1's (sin θ2/ε)^{1/2} ≤ R'^{a0/2};
   - the tangent/crossing threshold C′δ′;
   - the o(1) terms of F5(c).
3. Estimate the discretisation error. Each checked quantity exceeds its continuum value **at g_i** by at
   most Kϖ, with K absolute (Lipschitz constants of the pieces). `cost_upper` is continuous across the
   admissibility flips of its candidates (`review/check_cellmin.py`). The three places where the
   certificate uses g_{i+1} are handled as follows:
   - the sagitta term and the constraint X_L < 1 − g: shift both cell exponents by −2ϖ (cost ≤ 4ϖ);
   - the lower end lde(g_{i+1}) of the crossing range: classes with D ≤ lde(g_i) have the tangent
     class's S_T and S0b and fewer offsets, so the tangent piece at g_i dominates them;
   - the no-crossing case 1 − 2g_i < lde(g_{i+1}): it is impossible, since Δ ≥ C′δ′ > R'^{1−2g_i} ≥ Δ.
4. Close each interval. By C(a0, T), let η* be the least height with HIGH_c ≤ T. It exists because
   HIGH_c is continuous, nonincreasing and tends to E1 < T. Every layer strictly below η* has
   LOW_c ≤ T. The last layer [y_{j*−1}, y_{j*}] contains η*, and its discretised value is
   ≤ LOW_c(y_{j*−1}) + Kϖ with y_{j*−1} < η*. So each interval closes at T + Kϖ.
   - *Zero case* (r > R'/2), not part of C: the whole-box cell has P = (w − 2/5 + μ)/3 < 0 on [0, ϖ], so
     its cost is 0 (`review/zero_case.py`).
5. Count. There are O(ϖ^{−3}) cases, Lemma F5 is applied O(ϖ^{−2}) times, and each case contributes
   C·R'^{T+Kϖ} or a global bound R'^{T+o(1)} (the first alternative of F5). The total is R'^{T+o(1)}. ∎

The Gibbs lemma with λ = 1/α and n(τ) ≤ 2^{(T/4)τ + o(τ)} (T/4 < 1/α) then gives the rate α with an
additive term A = 1 + log2 Z_λ = O_α(1): the gap λ − T/4 = (a0 − T)/4 is fixed, and the o(τ) terms are uniform in the frame.

## 2. The continuum claim and its exact verification

`cont_model.py` is the continuum version of `cert_proj.py`, written over a generic number type:
- pointwise in g;
- offsets and per-sphere counts at the same height y;
- a continuous crossing depth D ∈ [lde(g), 1 − 2g];
- the adversary's X ≥ 0 and angle handled exactly as in `cert_proj.py` (G3′ = T − 4 + a0 + γ);
- prover cell costs = min(patch fibration, four explicit candidate cells with X_S = U0 and
  X_L ∈ {S0, (U0+R0)/2, U0, min(S0, (U0+R0)/2)}), which is an upper bound for the true cell cost.

**C(a0, T):** there is no (g, y) with 0 < g ≤ 201/100 and y > 0 such that
HIGH_c(g, y) > T and LOW_c(g, y) > T.

*Why this is the right claim.*
- HIGH_c(g, η) is nonincreasing in η, so "∃η*: HIGH ≤ T and LOW ≤ T for all y ≤ η*" fails iff some
  y has HIGH(y) > T and LOW(y) > T.
- HIGH → E1 as η → ∞, so η* exists.
- LOW_c is a maximum over the rank d ∈ {≤1, 2, 3, 4}, the piece (tangent, or crossing D) and X, of a
  minimum of route N and route P. The adversary's choices are existential. The prover's minima
  (route N vs P, the candidate cells, tube vs box bound) become conjunctions.
- So ¬C(a0, T) is a quantifier-free linear-real-arithmetic formula.

**Verification.** `verify_z3.py μ0 c` checks C for all μ ∈ (0, μ0] at once: μ is a variable, with
a0 = 8/5 + μ, b = E1 = (8 − 2a0)/3 and T = E1 + cμ. Results (`z3_outputs.txt`, z3 5.1.0 on school-server):

| μ0 | c | result |
|---|---|---|
| 1/10 | 1/2 | **unsat** → C holds for all a0 ∈ (8/5, 17/10] with T = 8/5 − μ/6 < a0 |
| 1/10 | 10/33 | **unsat** → the stronger T = 8/5 − (4/11)μ on all of (8/5, 17/10] |
| 1/10 | 3/10 | sat |
| 1/300 | 10/33 | **unsat** (tight) |
| 1/300 | 38/125 = 0.304 | unsat |
| 1/300 | 303/1000 | sat (counterexample at g ≈ 1.592, y ≈ 1.193) |
| 1/300 | 0 | sat |

The encoding is not vacuous: the sat/unsat threshold is at c* = 10/33, which is the value predicted
analytically in §3.

`crosscheck.py`: at 3000 random points the continuum pieces and HIGH are never stronger than the
reviewed certificate code (`cert_dioph` pieces with gl = gr, yp = yj, D_l = D_r). The one flagged point
is a 3·10⁻⁹ rounding in the certificate's rationalised witness cell; the continuum witness is exact and
admissible.

## 3. Where 10/33 comes from (analytic check; corrected after review)

Write a0 = a = 8/5 + μ and b = 8/5 − (2/3)μ. Box length exponent: 1 − b = −3/5 + (2/3)μ. Transverse
exponent: 1 − a0 = −3/5 − μ.

For 2/5 ≤ g ≤ 8/5 − (7/3)μ:
- **Patch.** The patch of every section is the whole box: S0 = 1 − b, U0 = 1 − a0.
- **HIGH.** The whole-box cell has sagitta term 2X_L − (1 − g), which gives
  HIGH − E1 = max(0, (g − 2/5 + μ − η)/3). So η0 = g − 2/5 + μ.
- **Offsets and per-sphere counts.** The tangent class has no offsets below η0. The deepest crossing
  class (D = 1 − 2g, patch δ × δ, box bound 0) gives pm = y + 2 − 2g once y ≥ η0 − 5μ.
- **Rank-3 routes.** Route P is flat (= E1) for X ≤ y − 4/5 + μ/3. Route N at that point is
  2y − g + 4/5 − μ/3 + pm. Their crossing value is E1 + (3y − 3g + 6/5 + μ/3)/2.
- **Balance.** HIGH (slope −1/3 in η) against this crossing (slope +3/2) meets at η* = g − 2/5 + μ/11,
  with value E1 + (10/33)μ.
- **This mechanism needs its flat point to be admissible:** X ≥ 0, i.e. y ≥ 4/5 − μ/3, i.e. g ≳ 6/5.
- **Below that, rank 2 binds, through the T inside G3′.** Put y = g − 2/5 + δ and T = E1 + cμ. Then N2 is
  flat at E1 + 3δ + μ/3 − cμ. Balancing it with HIGH − E1 = (μ − δ)/3 gives δ = 3cμ/10. Requiring the
  value to be cμ gives c(1 + 1/10) = 1/3, i.e. the same c = 10/33.
- **Where the value is attained.** The value is exactly E1 + (10/33)μ for
  g ∈ [4/5 − 7μ/11, 8/5 − 7μ/3], and smaller elsewhere. The lower end is G3′ ≥ 0. Float scan: [0.738, 1.366]
  against the predicted [0.736, 1.367] at μ = 0.1 (`review/fixed_g_opt_output.txt`,
  `review/float_brute_output.txt`).

At μ = 0 all of this holds with equality at y = η0 = g − 2/5 (`limit25b.py` in ../annulus). This
explains why the finite certificates lose exactly O(layer width) near 5/2.

## 4. What 5/2 means, and Proposition S

- **Target = one word per δ-cube.** δ = εR' is the tube radius, and R'/δ = R'^{a0} for every α. So the
  rate-α target count 2^{t/α} = R'^{a0} at t = αL is exactly the number of δ-cubes along the core.
- **Why hyperplanes stop at 5/2.** The 5-point determinant makes a δ-cube's points coplanar iff
  δ^5 R'^3 < 1, i.e. a0 > 8/5, i.e. α < 5/2. So α = 5/2 is where the level-1 boxes are exactly δ-cubes.
- **What exceeding it requires.** Beyond 5/2 the bound must come from most δ-cubes being empty at
  T-count 5L/2, where the volume law predicts 2^{L/2} words against 2^L cubes.

**Proposition S.** Suppose that for ε in a sequence tending to 0 there are:
- a Clifford+T word G;
- words W_u (u ∈ U ⊂ Z_{Q_ε}) of T-count ≤ τ_ε = (5/2 + o(1))L with d_proj(W_u, R_z(Δu)G) ≤ ε/2;
- U has maximal gap γ ≤ Q_ε 2^{−(1−o(1))L}.

Then for r = 2 and every m there is a deterministic zero-error CW process with
S ≤ m(log2 γ + O(1)) and T_pre ≤ m τ_ε. Hence no bound T̄pre ≥ α(m log2 K − S − o(mL)) with α > 5/2
holds along that sequence.

*Proof.*
1. **Round 1.** For each coordinate, take u ∈ U with ρ := a1 − u ∈ [0, γ). Emit W_u, and keep only ρ in
   the snapshot. G is a public constant.
2. **Round 2.** Emit a word V with d_proj(V, G^{−1} R_z(Δ(ρ + a2))) ≤ ε/2. This costs nothing in T_pre.
3. **Correctness.** W_u V is within ε of R_z(Δu) G G^{−1} R_z(Δ(ρ + a2)) = R_z(Δ x), by unitary
   invariance and the triangle inequality.
4. **Rate.** The process forgoes m(log2 K − log2 γ − O(1)) = m(1 − o(1))L bits at committed cost
   m(5/2 + o(1))L. Compare with Theorem U's form of the bound. ∎

Theorem U and Proposition S together say that the true exchange rate is exactly 5/2 if such dense tubes
exist. The converse is weaker than "sparse by a power is necessary":
- Proposition S only excludes tubes with maximal gap ≤ Q·2^{−(1−o(1))L}. A tube with 2^{(1−o(1))L} words
  concentrated on part of the circle does not satisfy its hypothesis.
- A randomised variant stores one shared random rotation offset in the snapshot (log2 Q bits in total)
  and codes ρ_j by Elias-δ. It replaces the maximal gap by the length-weighted mean
  Σ_i (γ_i/Q) log2 γ_i, which is ≥ log2(Q/|U|) by Jensen (`review/REVIEW.md` §4).
- So what a proof of a rate above 5/2 must exclude is tubes at T-count 5L/2 with small weighted mean
  log-gap. That condition sits between "no dense tube" and "sparse by a power".

**Heuristic (not a theorem): higher-degree auxiliary surfaces do not help.** A K-rational quadric
through the points of a longer box can be a thin tube of radius ≤ δ around the core. On such a surface
the level-2 cells on which 5 points are cohyperplanar have length s^{−2/3}R'^{−1} (s = δ). The total
count is then R'^{2 + (2/3)(1−a0)} = R'^{(8−2a0)/3}, exactly the degree-1 box count. The degree-2 gain
(level-1 ceiling 28/11) is lost, and a helix model with degree-2 curves at level 2 and Bézout at level 3
also breaks even at 5/2. The obstruction is the configuration "points on a thin rational quadric tube
around the core", i.e. one point per δ-length. This sharpens the 8/3 ceiling of `../push`, which
assumed a free level 2.

## Files
- `cont_model.py` — continuum lemma set (float/z3 backends); `python3 cont_model.py μ c` scans g, η.
- `verify_z3.py` — the exact check; `z3_outputs.txt` — outputs; `z3_version.txt`.
- `crosscheck.py`, `crosscheck_output.txt` — agreement with the reviewed certificate code.
- `../annulus/limit25.py`, `limit25b.py` — the structure at μ = 0.
