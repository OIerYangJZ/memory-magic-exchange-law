# Linnik basic lemma / entropy / non-concentration: S-arithmetic route (agent_linnik, 2026-10-01)

**Verdict: negative on all three items.** Every mechanism in this family gives a tube bound with exponent
s ≤ 3/2 in μ_t(T_ε) ≤ ε^s, and most give much less. Q1 needs s = 3/2 + κ. The determinant method already gets
s = 3/2, and volume is s = 2. The failure points are exact and listed below. Throughout: ε = 2^{-2t/5}, Q ≍ 1/ε
cells, R' = 2^{t/4}, and N = N_t(T, ε) is the number of points in the tube. d_i is the PU(2) distance in σ_i.
"Prefix c" means two words share exactly c leading T-syllables (Gromov product c in the tree). m := t − c.

## 0. Three elementary lemmas (proved; the first is also checked numerically)

**L1 (separation).** If W₁ ≠ W₂ have T-count t (more generally, the same parity and T-count ≤ t) and prefix c, then
d₁(W₁,W₂) ≥ c₀ 2^{-(t-c)/2}. With consistent lifts (same sign), also d₁·d₂ ≥ c₀ 2^{-(t-c)/2}.
*Proof.*
- Write W_i = P U_i, so d_i(W₁,W₂) = d_i(U₁,U₂).
- Take numerators u_i with nrd u_i = n_m. If the parities agree but the T-counts differ, rescale by a power of √2;
  the ratio n_s/n_{s'} is a power of 2 when s ≡ s' mod 2.
- Then x = u₁ ∓ u₂ ∈ O∖0, so |σ₁x|·|σ₂x| ≥ 1, while |σ_i x| ≲ 2^{m/4}. ∎

*Numerics.* /tmp script, exhaustive over all 36·2^m words, floats. min d₁ · 2^{m/2} = 0.83, 1.17, 1.17, 0.90, 1.27,
1.27, 0.97, 0.74, 1.05, 1.05 for m = 2..11. So the exponent is sharp.

*Corollary.* Every non-identity E ∈ Γ of **even** T-count k has d₁(E, I) ≥ c₀ 2^{-k/4}. Factor E = U₁⁻¹U₂ at the midpoint of
its geodesic. This is twice the exponent of the naive bound 2^{-k/2} that comes from E's own numerator. For odd k the two halves
have different parity, and only 2^{-k/2} is proved. Below, only L1 for same-parity words is used.

**L2 (joint S-arithmetic basic lemma).**
- *Bound.* The number of pairs (W₁ ≠ W₂) at level t with prefix exactly c, d₁ ≤ δ₁ and d₂ ≤ δ₂ is
  ≤ 2^{c+o(t)} (1 + Δ2^m)(1 + Δ²2^m), where Δ = δ₁δ₂.
- *Vanishing.* The count is 0 if Δ < c₀2^{-m/2} (same-sign lifts). For mixed signs it is 0 unless each δ_i ≥ c₀2^{-m/2}.

*Proof.*
- (W₁,W₂) ↦ (P, D = U₁⁻¹U₂) is 24-to-1, and D has T-count exactly 2m.
- Write 2d = A + Bi + Cj + Dk over O_K, with nrd d = n_{2m}.
- The B-coordinate lies in a box of area ≍ δ₁δ₂2^m, and A lies in a box of area ≍ δ₁²δ₂²2^m.
- Given (A, B), C² + D² = r is a norm equation in Z[ζ₈] (class number 1, finitely many relative-norm-one units).
  So the divisor bound gives 2^{o(m)} solutions.
- Vanishing follows from L1. ∎

*Summed over c (σ₁ only).*
- **P_t(δ) := #{W₁≠W₂ : d₁ ≤ δ} ≤ 2^{o(t)} 2^{2t}δ³** for δ ≥ c₀2^{-t/2}, and P_t(δ) = 0 below that.
- This is the Haar (volume) law at every scale: μ_t has optimal L²-dimension and is 2^{-t/2}-separated.
- Numerics: at δ = 3·2^{-m/2} the pair counts are within a factor 1.3 of Haar (m = 9, 10, 11).

**L3 (one-point tree lemma).** Any set of σ₁-diameter r contains ≤ C(1 + r²2^t) words of T-count t.
*Proof.* If it contains M points, pigeonhole them into the 3·2^{c−1} prefix classes. Two of them share a prefix
c ≈ log₂(M/C), and L1 then gives r ≥ c₀2^{-(t-c)/2}. ∎

## 1. What the basic lemma gives on tubes (item 1)

| method | bound on N | exponent s in μ_t(T_ε) ≤ ε^s | at critical ε |
|---|---|---|---|
| L3, tube = Q cells | N ≤ C 2^t ε | 1 | 2^{3t/5} |
| pairs + L2 (any 2-adic/σ₂ refinement) | N ≲ 2^tε + t2^{t/2} | 1 | 2^{3t/5} |
| Ramanujan (R) | 2^tε² + t2^{t/2} | 5/4 | 2^{t/2} |
| determinant (D) | ε^{-1} | 3/2 | 2^{2t/5} |
| **needed (Q1)** | ε^{-1+κ} | 3/2 + κ | 2^{(2/5−κ)t} |
| volume | 2^tε² | 2 | 2^{t/5} |

**Pairs.**
- All pairs in T have W₁⁻¹W₂ in the 2ε-tube of the one-parameter group H₂ = G₂⁻¹R_zG₂, at T-count 2(t−c).
- So N² ≤ Σ_c 2^c · 24 · N^{(D)}_{2(t−c)}(H₂, 2ε).
- Insert Ramanujan for the D-tubes, or even the conjectural volume law 2^{2m}ε². Either way N² ≲ 2^{2t}ε² + t²2^t.
- This ceiling is intrinsic: a tube with N points needs ≈ N² pairs, while the D-tube has Haar mass 2^{2t}ε².

**2-adic restriction.**
- Cauchy–Schwarz over 2^c prefix classes gives N² ≤ 2^c · #{pairs sharing ≥ c}.
- But #{pairs sharing ≥ c} ≲ Σ_{c'≥c} 2^{c'}2^{2(t−c')}ε² = 2^{2t−c}ε².
- Result: N² ≲ 2^{2t}ε² **for every c**. No gain.

**σ₂ restriction.**
- Pigeonholing S³ into η-caps gives N² ≤ η^{-3}(N + #{pairs in T with d₂ ≤ η}).
- The main term η^{-3}·2^{2t}ε²η³ is unchanged.
- The only possible gain is that the off-diagonal pairs vanish. By L2 that needs Jεη < c₀2^{-t/2} for pairs within
  J cells, i.e. Jη < 2^{-t/10}. This is incompatible with η ≥ M^{-1/3} ≥ J^{-1/3}.

**L² over all tubes.**
- Σ_{T ∈ ε-net} N(T)² ≲ 2^{o}(2^tε^{-2} + 2^{2t}). This is computed from P_t(r) and the fact that two points at
  distance r lie in ≍ r^{-2} common tubes.
- So all but a fraction ε^{1−2κ} of the ≍ ε^{-4} tubes satisfy N ≤ ε^{-1+κ}.
- This is an almost-all statement. By NOTES §2/§7 the adversary chooses the frame, so it is worthless for Q1.

## 2. Gap principle across cells (item 2)

**Setup.** In a window of J consecutive cells with M ≤ J points (full occupancy means M ≈ J), split the points into
K classes. If M > (k−1)K, some k points share a class. A contradiction needs the shared class to force a separation or
determinant incompatible with lying in J cells. Write the gain as a function of K.

| class type (K classes) | what sharing gives | needed for a contradiction | deficit (rows 5–6 in determinant units) |
|---|---|---|---|
| prefix c (K = 2^c), pairs | separation 2^{-t/10}K^{1/2} cells (L1) | > J ≥ K | K^{1/2}2^{t/10} |
| σ₂ cap η (K = η^{-3}), pairs | separation 2^{-t/10}K^{1/3} cells | > J | K^{2/3}2^{t/10} |
| prefix × suffix (K = 2^{c₁+c₂}) | identical to prefix c₁ + c₂: d = p(u₁−u₂)s | — | same as prefix |
| prefix × σ₂ cap | ≤ 2^{-t/10}K^{1/2} cells | > J | needs η > J·2^{t/5} > 1 |
| 5-point det + prefix | det ∈ 2^{2c}O_K/16: gain K² | det ≍ J³ must be beaten, M ≤ J | factor J |
| 5-point det + σ₂ cap | gain η^{-5} = K^{5/3} | K³ | factor J^{4/3} |
| heap identity (BG-type, §3b) + prefix | separation 2^{-3(t−c)/2}, size Jε² | c > 7t/15 + (2/3)log J | needs log J > 7t/5 |

*Why the 5-point determinant rows fail.*
- k = 5 points sharing prefix c make the determinant nonvanishing only if J³ ≳ 2^{2c} ≈ (M/12)².
- So cohyperplanarity is never forced across more than O(1) cells.
- This is the 2-adic analogue of NOTES §5(ii) "J > J^{9/5}"; route (e) of NOTES §3 is the same loss seen from the other side.

**Why nothing survives.**
- In a 3-regular tree, splitting into K classes buys at most K^{1/2} in the archimedean product |σ₁d|·|σ₂d|, and at most
  K² in a 4×4 determinant.
- A gap principle along a 1-dimensional tube needs gain K¹ for pairs, and K³ for determinants.
- There is also a fixed offset: at c = 0 the words are separated only at 2^{-t/2} = ε·2^{-t/10}, which is below one
  cell. Removing that offset by sharing costs K = 2^{t/5} classes, and that is exactly the volume E = 2^tε².

**Global pigeonhole.**
- Prefix classes of length c = t/5 + 2k have their points ≥ 2^k cells apart. With N ≈ Q, a class has ≈ 2^{t/5−2k}
  points spread over Q·2^{-k} slots: occupancy 2^{-t/5−k}, so there is no tension.
- The class decomposition is just the Hecke identity N_t(T) = Σ_{P∈Λ_c} N_{t−c}(P⁻¹T).
- The tubes P⁻¹T are left translates, so they are Clifford-parallel and impose no mutual constraint.

## 3. BG flattening, product theorems, projections, Furstenberg (item 3)

**(a) Flattening is already optimal and blind to tubes.**
- L2 gives ‖μ_t^{(ε)}‖₂² ≍ ε^{-3}2^{-t} and ‖μ_t∗μ_t‖₂ ≍ 1, so μ_t flattens in one step.
- *No-go (rigorous).* Add ε^{-1} points along any tube, one per cell, to a set obeying μ(B(x,r)) ≤ Cr^s for r ≥ ε:
  - It still obeys the same bound for every s ≤ log₂(2^t)/log₂(1/ε) = 5/2, since r/ε ≤ 2^t r^s for r ≥ ε.
  - The pair energy at scale r changes by ≲ ε^{-2}r + ε^{-1}·2^{t+o(t)}r³ ≪ 2^{2t}r³.
- So no theorem whose hypotheses are ball non-concentration or L²-dimension (flattening, Bourgain/He projection
  theorems, discretised sum-product) can decide Q1.

**(b) BG escape from subgroups (Diophantine commutators), computed exactly.**
- *Heap identity.* For W₁, W₂, W₃ in T, put V = W₂W₁⁻¹W₃ and V' = W₃W₁⁻¹W₂. Then d₁(V, V') ≤ 6ε, since the torus is
  abelian.
- If the three points lie within J cells, the first-order terms cancel and d₁(V, V') ≤ CJε². The residual is the
  commutators [sn, f] and [f, f].
- V and V' have T-count ≤ 3t and the same parity, so V ≠ V' ⇒ d₁ ≥ c₀2^{-3t/2} (L1).
- Forcing V = V' for all triples puts the quotients W₁⁻¹W_i in one centraliser K(x)^×. That is S-units of a CM
  quadratic extension modulo K, of rank ≤ 1, so the set has O(t) points.
- *Result: N_t(T, ε) = O(t) for ε ≤ c·2^{-3t/2} = cR'^{-6}.*
- *Window version.* A J-cell window lies in one abelian coset if ε² J < c2^{-3t/2}. This is a genuine gap principle, but
  only for ε ≲ R'^{-3}.
- *Prefix refinement.* Summed over windows it gives N ≲ t·ε^{2/3}2^t, which is worse than L3.
- *Comparison.* The determinant method gives cohyperplanar boxes of J = R'^{5a/3−8/3} cells, against
  R'^{2a−6} for commutators, at ε = R'^{-a}. So BG-type escape is dominated for every a ≤ 10 and is useless near
  a = 8/5.
- *Degree count.* The heap relation has σ₁-size Jε²R'^6, against J³ε⁵R'^8 = J³ for the 5-point determinant. At
  critical it is off by R'^{14/5}.

**(c) Product theorems** (BG, de Saxcé, Breuillard–Green–Tao).
- They take non-concentration on cosets of proper closed subgroups as an *input* and output growth of A·A·A.
- A rich tube is by definition A ⊂ ε-neighbourhood of a torus coset, i.e. the exceptional case. A⁻¹A sits in the
  ε-tube of H₂ and does not grow; the theorems say nothing about |A|.

**(d) Translates of a rich tube.**
- *Left or right translates are parallel.* γT = (γG₁G₂)H₂ is a left coset of the same H₂: a Hopf fibre.
  T·γ is a right coset of H₁. Disjoint parallel tubes carry no incidence geometry.
- *Amplification is exactly neutral.* Apply the Hopf-form Ramanujan bound to the union of the 2^j caps:
  - c2^jN_t ≤ C2^{t+2j}ε² + C(t+j)2^{(t+j)/2}·ε^{-1}(2^jε²)^{1/2}.
  - Hence N_t ≤ C2^{t+j}ε² + C(t+j)2^{t/2}.
  - Signal and spectral error both scale by 2^j. This is the same square-root barrier as NOTES §3(c).
- *Two-sided translates γTγ'* with γ ∈ Λ_j, γ' ∈ Λ_{j'}:
  - They give up to 2^{j+j'} rich tubes in general position at level t' = t + j + j'. Great circles ↔ lines of RP³,
    so this is a Furstenberg/Kakeya configuration in R³.
  - But the total incidence count is 2^{j+j'}ε^{-1} = |Λ_{t'}|·2^{-3t/5}.
  - Any Furstenberg or Kakeya lower bound for the union is ≤ (number of tubes)·ε^{-1}. So it can never exceed
    |Λ_{t'}|.
  - Even a hypothetical Szemerédi–Trotter bound #(r-rich frames) ≲ n²/r³ + n/r, for the 2^t surfaces
    {(R(W)x₂, x₂)} in S²×S², fails. The linear term n/r = 2^{t'}ε exceeds 2^{j+j'} for every j, j'.
  - Without that term it would only give N ≤ 2^{2t/3} = ε^{-5/3}, i.e. rate 3/2.
- *Ren–Wang* is planar. In Hopf form the tube becomes a single cap, and the question becomes multiplicity, not
  incidence.

## 4. Bottom line

The S-arithmetic Linnik lemma (L1–L3) is perfectly sharp for balls: volume law at every scale, plus exact
2^{-(t−c)/2} separation. But it controls tubes only in L², or through covering by cells. That yields N ≤ C2^tε
(s = 1, rate 2) for every tube, and volume only for almost every tube.

Tree classes are too cheap to force a gap principle: the gain is K^{1/2} per K classes, against the K¹ needed. They
also carry the fixed offset 2^{t/10} = E^{1/2}.

Sparse-regime tools either require ball/subgroup non-concentration as input, which cannot see an ε^{-1}-point tube, or
work by Diophantine commutators, which are effective only for ε ≤ R'^{-3}…R'^{-6}.

So nothing here goes past the determinant's s = 3/2. The missing ingredient is the same one NOTES §§5–7 identify:
single-norm arithmetic information about a generic torus. No information about non-concentration of μ_t, joint or
2-adic, can supply it, because every such statement is also true for Λ_t with a full tube added.
