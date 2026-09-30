# Referee report on `research/conditional3/NOTES.md` (2026-09-30)

Scope: the rigorous covering principle (§1), the BKS transplant (§2), the sphere-section cluster heuristic (§3), and
the numerics (§4). Checked against main.tex (Lemma gibbs, the proofs of Thms rate2/rate3seg/main2, Lemma E1, E5,
Remark barrier) and BKS (arXiv 1609.06097, §§2–5). New checks were run on school-server in `~/review_c3`. Scripts and
outputs are in this directory. Nothing outside `review/` was modified.

## Verdicts

| Claim | Verdict |
|---|---|
| (1) Lemma 1, Thm 1 (B_α ⇒ rate α), Thm 1′ (★_θ ⇒ rate 1/θ), minimal spacing d ≥ 2^{−τ}/4 | **VALID WITH FIXES** (minor) |
| (1′) "Conjecture H ⇒ (★_{1/3})", "strictly weaker", the "unconditional ladder" | **VALID WITH FIXES** (H is grid-restricted; ranges) |
| (2) twisted Linnik over Q(√2) gives only N^{1/2} for tubes and ρN^{1/2} for balls | **VALID as a heuristic**; the no-go wording is **OVERSTATED** |
| (3) clusters of size 2^{τ(2/3−5/(3α))}; (B_α) false for α > 5/2; threshold = R'^{−8/5}; compatible with H | exponents and mechanism **VALID as a heuristic** (now confirmed arithmetically at τ ≤ 26); "(B_α) is false" and "no frame-free hypothesis" are **OVERSTATED**; the constants are off by ≈ 3–70× (median 14) |
| (4) code (`balls.cpp`, `ballpts.cpp`, `clusan.py`) | **VALID** (no bugs found; counts independently reproduced) |
| (4) readings of the tables and clusters | **OVERSTATED / partly WRONG** (the statistic, the control and the τ = 26 identification) |

## 1. Rigorous part

All steps check out against main.tex:
- **Haar measure.** 2(θ − sin θ cos θ)/π with θ = 2 arcsin(ρ/2) is right: two antipodal caps of S³. Its value is ≈ 0.424ρ³ for small ρ and equals 1 at ρ = √2, so c₁ ≈ 0.354 and c₂ ≈ 0.424.
- **Lemma 1.** Bi-invariance holds, and d(R_z θ, R_z θ′) = 2|sin((θ−θ′)/4)| ≤ |θ−θ′|/2. Using θ ∈ [0, 2π) is enough in PU(2). The packing/covering step is right.
- **Theorem 1.** Lemma gibbs is applied with λ = 1/α and τ* = ⌈α log₂K⌉ ≥ λ^{−1}log₂K. From ε < sin(π/2Q) < π/2Q and K < Q we get log₂K < L + log₂(π/2). Hence for τ ≤ τ*, r = 2^{−τ/α} ≥ (2/π)2^{−1/α}ε ≥ 0.31ε, and the O(1) extra balls are justified. The final bound α[(1−2δ)m log₂K − 2m h₂(δ) − S − mA] is exactly the rate-2 bookkeeping with ½ → 1/α.
- **Theorem 1′.** For θ ≥ 1/3, ε²2^{(1−θ)τ*} = O(1).
- **Minimal spacing.** Every step is correct:
  - y = ⟨w₁,w₂⟩ = ½trd(w₁w̄₂) ∈ ½O_K and x = n₁n₂ − y² ∈ ¼O_K.
  - x is totally ≥ 0 (Cauchy–Schwarz at both definite places), and x = 0 iff W₁ = W₂ projectively.
  - σ₁(x)σ₂(x) ≥ 1/16.
  - d² = 2 − 2c ≥ 1 − c² = σ₁(x)/σ₁(n₁n₂), since (1−c)² ≥ 0.
  - nrd = (2+√2)^t·unit = (√2)^t·unit has |N| = 2^t, so N(n₁n₂) ≤ 4^τ.

**Fixes.**
1. **Conjecture H counts grid-restricted words** (∃d ∈ Z_Q, Q ≤ Q_ε). It does not bound the tube count n(τ) of Lemma 1. So "H ⇒ (★_{1/3})" needs a covering step. Four grids with Q = Q_{2ε}, shifted by 2π/(4Q), cover T_ε by grid-neighbourhoods of accuracy 2ε, which gives n_tube(τ, ε) ≤ 4(c₀ + c₁τ + 4c 2^τ ε²). Alternatively, state (★) for the grid-restricted alphabet, which is all that the Gibbs step uses (as in the proof of Thm main2), and carry the hypothesis Q ≤ Q_ε.
2. **"Strictly weaker"** → "weaker (implied by H, not known to imply it)".
3. **The ladder.** Theorems elem and dioph, and uniform52, give θ = 0.4513, 0.4111 and → 2/5 only for τ ≤ (20/9)L, τ ≤ (17/7)(L+2) and τ ≤ (5/2 − μ)L respectively. So state (★_θ) for τ ≤ τ* only, which is all that Lemma gibbs needs.
4. **Cosmetic.** Use ε + r/4; say explicitly that for r ≳ 0.8 the covering is trivial (the volume bound is used only for ρ ≤ √2).
5. **Worth adding: (B_α) for α < 2 is unconditional.**
   - After rescaling by powers of √2, the words of one T-count parity lie on the sphere nrd = n_τ.
   - The 3-wedge of four points in a ball of radius ρ has σ₁·σ₂ size ≤ Cρ³R'^6, and its coordinates are in (1/8)O_K. So it vanishes for ρ < c·2^{−τ/2}.
   - The points therefore lie on a circle, and Lemma E5 bounds them by 2^{O(τ/log τ)}.
   - So Theorem 1 already gives every rate α < 2 without (R). This sharpens the §2 remark that twisted Linnik "gives (B_2)": that is no better than the elementary bound.

## 2. The BKS transplant

**What checks.** The scaling bookkeeping matches BKS.
- The mnemonic E ≈ M·#C/N(Q) reproduces their dominant error: M = ε³N, #C ≈ 1/ε (|ĉ₁,₂,₃| ≲ N^δ, |ĉ₄| ≲ N^δ/ε, Lemma 4.1 and #𝒞 = O(ε^{−1}N^{4δ})), and Q = ε√N give εN^{1/2}. This is exactly the large-q contribution E_i(w,R) ≪ ε⁴N²/(Q²R) (the step after (5.5)).
- **Ball over K.** N(Q) = ρR² after unit rebalancing. Q₁ = ρR < 1 is harmless because only the norm matters. This gives E ≈ ρN^{1/2}.
- **Tube.** The torus weight has ⊥-width εR and in-plane radial width ε²R, so Q₁ = εR. The 2-D annulus Fourier decay is (R|c_Π|/q)^{−1/2}. The weighted dual count is Σ_{|c_Π|≤1/ε}(ε/|c_Π|)^{1/2} ≈ 1/ε, from ≈ 1/ε² vectors of mean weight ε, not "1/ε vectors". So E ≈ N^{1/2}, independent of ε. This agrees with the spectral bound (k+1)2^{k/2} of Thm tube.
- The Bessel phase e(±R|c_Π|/q) has twist parameter |c_Π|/√F(c) ≤ 1, inside the range of Conj. 1.1.
- The claim that thinner σ₂-shells worsen #C/Q₂ by 1/Δ checks at the same heuristic level, once the |ζ|^{−3/2} decay of the sphere is included.

**What is overstated.**
- The Verdict (items 1–2) states the N^{1/2} error as a fact. It is a heuristic about one natural implementation: Linnik-type cancellation in q for each dual vector c, then |·| over c.
- "Beating N^{1/2} needs cancellation between different c … not a statement about Kloosterman sums" is fair as a description of this framework. But the obvious candidate, bilinear sums Σ_c Σ_q a_c S(n, F*(c); q) via Kuznetsov / the spectral large sieve (Deshouillers–Iwaniec), *is* a statement about Kloosterman sums. It leads back to sums of Hecke eigenvalues along the theta lift, i.e. the spectral side of Remark barrier. Say "we see no input that is not equivalent to the count", not "cannot".
- The same applies to §5's "cannot be converted into rate > 2 by the circle or spectral method".
- **Ternary bookkeeping.** v = w k′ w̄ has nrd n²·nrd(g₂), of norm 2^{2τ+h}. So there are about 2^{τ+h/2} points (not 2^{τ+h}), with square root 2^{τ/2+h/4}. The conclusion "≥ 2^{τ/2}" is unchanged. The identity also needs an arithmetic frame G₂ ∈ Γ, so it does not cover general frames.

## 3. The cluster heuristic, checked by exact arithmetic

`check_sections.py` enumerates, exactly in O_K-coordinates:
- Z = {z ∈ O : Re z = s, nrd z = nl};
- all g ∈ O with nrd g = l;
- |Y_{g,s}| = #{z : z g/l ∈ O}.

`ballpts` counts the words directly.

**Confirmed.**
- **Event probability.** It is lattice counting and is well calibrated: at μ = 16, τ = 20–26, H ≤ 40, there are 36 events against 39.4 expected. At μ = 4 (τ = 16–26, H ≤ 60): 52 against 61. At μ = 1: 16 against 24.
- **Factorization.** Each P-primitive z gives exactly 48 pairs (w,g): mean |Y_g| = 48#Z/#g for all 83 events with N(l) odd. When P | l the P-parts of z and g interact, and the count is off by a factor 2–4 (21 events). The ≈ H centres share the words uniformly for prime l, by the left O¹ action.
- **Local conditions.** D is unramified at every finite place and O has class number 1. So there are no local or spinor obstructions: every event had #Z > 0.
- **The factor ½ in s and the balancing of odd τ** do not affect the box area √N(nl).
- **The optimisation and the threshold** are correct arithmetic. ρ^{5/3}2^{2τ/3} ≥ 1 ⟺ ρ ≥ 2^{−2τ/5} = R'^{−8/5}, and the 5-point ball version of Lemma E1 (curvature ρ²R') does give coplanarity exactly below R'^{−8/5}. That the two thresholds coincide is correct arithmetic, not an explanation.

**Problems.**
- **(a) Constants.** Over the 104 census events, #Z/N(m)^{1/2} ranges over 18–42 (median 29). The largest per-centre size at T-count τ is 3–74× the formula ρ2^{τ/2}/√H (median 14). For example, τ = 25, N(l) = 7: 126 words against the formula's 5.5 (the author's in-progress `arith_tau30–60.txt` on the server show similar ratios). The exponents are unaffected, but no numerical "prediction" in §4 is calibrated.
- **(b) Double counting.** l and 4l (g ↦ 2g) give the same centres and sections, and P-power multiples reuse the same z at other τ. For example, the τ = 21 section S = (73,32) reappears at (τ, l) = (17, 4) and (19, 2). Count primitive g only. This is a constant factor.
- **(c) Rigour.** The lower bound r(s,m) ≫ N(m)^{1/2−ε} is Siegel, hence ineffective. The existence of near-square events at H ≈ H_min is a random-model assumption (equidistribution of s² ≈ nl in σ₁ as l varies). So "(B_α) is false for every α > 5/2" is heuristic and must be labelled so.
- **(d) What is actually shown.** The defensible conclusion is narrower than "no frame-free hypothesis".
  - Heuristically, M(ρ) := max_ξ #(Λ_τ ∩ B(ξ,ρ)) ≳ max(1, ρ³2^τ, ρ^{5/3}2^{2τ/3}).
  - So every ball-covering bound n ≤ M(ρ)/ρ (over any ρ ≥ ε) is ≥ 2^{2τ/5}.
  - That is: Lemma-1-type arguments cannot beat 5/2. Frame-free structural hypotheses (few sections per ball, height-dependent section bounds, or averaged/L² ball counts) are not excluded.
- **(e) Compatibility with H.** It is checked only for tubes with ε ≥ r. For ε < r one needs equidistribution on Y (Duke-type). Then a tube of radius ε meets ≲ s(ε/r)^{3/2} cluster words, which is compatible with H. Add this.
- **(f) δ-cubes.** "The one point per δ-cube barrier describes the problem" is overstated. Rate > 5/2 needs most cubes *along one tube* to be empty. A few heavy cubes do not contradict that, as §5 itself says.

## 4. Numerics

**Code.** No bugs found.
- The Matsumoto–Amano enumeration is right: left-multiplied syllables HT and SHT, an optional leading T, and the Clifford on the right. It gives 3·2^τ − 2 prefixes (confirmed in `balls_err.txt`).
- The SU(2) forms of −iH, e^{−iπ/4}S and e^{−iπ/8}T are right, and so are `mul`, `inv` and `ip`.
- **Right-Clifford reduction.** x ∈ F iff v[h] = vmax. A point y within ρ of x ∈ F has |⟨y,1⟩| ≥ vmax(y) − 2ρ, so MARGIN = 2ρ_max suffices.
- The sign normalisation is safe, since every kept point has x₀ > 0.8. The 3-D grid of cell size ρ_max with ±1 neighbours is complete.
- Independent confirmation: `ballpts` at ξ_g (τ = 25, ρ = 0.0024989) returns exactly the 126 words of the arithmetic section, and at radius 0.0026 (τ = 26) the 126 + 219.

**Readings that do not hold.**
1. **Wrong statistic.**
   - The code computes the maximum over *word-centred* balls, not max_ξ. Section clusters sit around Hecke points ξ_g = g/|g|, which are not Clifford+T words when N(l) is odd, so a word-centred ball captures only a cap.
   - At τ = 25, μ = 16 the table says 49, but B(ξ_g, ρ) with g = ((3+2√2)+i+j−k)/2 contains **126** words (Haar mean 16).
   - The census finds B(ξ_g, ρ) ≥ 64, 168, 41, 48, 48, 126, 48 for τ = 20…26 (H ≤ 40).
2. **The τ = 26 row is misidentified.**
   - The listed "centre" is a word at distance ≈ 0.0023 from ξ_g.
   - The two sections lie at 1.13ρ and 1.28ρ from ξ_g. Their per-centre sizes are 126 (T-count 25) and 219 (T-count 26), and B(ξ_g, ρ) itself is **empty**.
   - They are not τ = 26 events at radius ρ. They are the same z-sets as the events at (τ, l) = (25, 3−√2) and (24, 2(3−√2)), since n_τ l is the same: 2^12(4+√2) and 2^13(3−√2). This also explains why two T-count classes share one g.
   - The actual τ = 26 events at radius ρ with H ≤ 40 are at N(l) = 17 and 31 (2 found, 3.9 expected).
   - "H ≈ 7.4 matched" and "0.09 per l … one was found" are therefore not tests. Given the constants of 3(a) and the sparse set of norms {1, 2, 4, 7, 8, 9, 14, …}, any small H would "match".
3. **τ = 21 is the l = 1 exceptional case** (N(m) = 42, 168 words, exact), not the generic mechanism. §3 calls these rare.
4. **Under-symmetrised control.**
   - Λ_τ is invariant under left and right Cliffords and under inversion (order 1152). The `--random` control has only the right action, i.e. 48× more independent centres, which biases every "words vs random" comparison.
   - With a matched control (`balls_sym.cpp --symrandom`, two seeds, `symrandom.txt`), the maxima (μ = 1/4/16) are e.g. τ = 23: 13–18 / 26–31 / 48–51; τ = 25: 12–15 / 24–28 / 48–49; τ = 26: 12–16 / 21–24 / 52–53.
   - So the word maxima exceed the control only at τ = 21 (μ = 16) and τ = 26 (μ = 4, 16). "Exceed at every μ from τ = 21 on, and grow" is not supported.
   - "Repulsion at μ = 1, τ ≤ 20" is confounded. The matched control's maxima come from orbit collapse near fixed-point sets (e.g. near the y-axis torus), and the words show voids around low-height points: B(ξ_g, ρ) is empty at τ = 26.
   - Table-level statistics are therefore inconclusive. The section census is the real evidence.
5. **Small-τ maxima are the same mechanism.** The "Clifford/identity structures" at τ ≤ 15 (and e.g. τ = 19, μ = 4 at a low-T-count word) are the same section mechanism with l = 1 or a P-power, not a separate effect.

## 5. Required fixes (summary)

- **§1:** fixes 1–5 above.
- **Verdict and §2, §5:** label the N^{1/2} statements as heuristic, and soften "cannot / no Kloosterman input". Fix the ternary exponent and note that it needs an arithmetic G₂.
- **§3:**
  - state the constant (#Z ≈ 30·N(m)^{1/2}; per centre ≈ 48#Z/r_O(l) = #Z/σ(l) for N(l) odd);
  - count primitive g only;
  - label "(B_α) false" as heuristic;
  - replace "no frame-free hypothesis" with the ball-covering statement of 3(d);
  - add the ε < r case to the compatibility check with H.
- **§4:**
  - say that the statistic is word-centred;
  - replace the τ = 26 row and the "prediction matched" sentence with the exact section data (`known_output.txt`, `census_*.txt`);
  - use a symmetry-matched control or drop the "repulsion / exceed / grow" readings;
  - call τ = 21 the l = 1 exception.

**Safe for a paper remark:**
- Lemma 1 / Thm 1 and Thm 1′ with the fixes, and the minimal spacing lemma.
- The elementary (B_{2−}), and the 5-point coplanarity of balls of radius < c·2^{−2τ/5}.
- A labelled heuristic: "near rational points of low height the words lie on sphere sections; heuristically max ball counts at radius 2^{−τ/α} grow like 2^{τ(2/3−5/(3α))} for α > 5/2, so bounds that cover a tube by balls cannot give more than 5/2; this is why Conjecture H is stated for tubes."
- One exact example, e.g. "126 words of T-count 25 in the ball of radius 0.0025 (Haar mean 16) around g/|g|, g = ((3+2√2)+i+j−k)/2, all on ⟨w,g⟩ = (357+256√2)/2".

**Not safe as written:** "No frame-free hypothesis can give more than 5/2", "(B_α) is false for α > 5/2", "twisted Linnik cannot give more than rate 2", and the §4 readings.

Files: `check_sections.py`, `known_output.txt`, `census_mu{1,4}_H60.txt`, `census_mu16_H40.txt`, `balls_sym.cpp`, `symrandom.txt`, `bp_25_0.0024989.txt`, `bp_26_0.0026.txt` (server copies are in `~/review_c3`).
