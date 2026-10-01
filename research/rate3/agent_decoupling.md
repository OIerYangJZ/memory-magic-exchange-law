# Decoupling, restriction, polynomial and Huxley-type methods for Q1 (agent report, 2026-10-01)

**Verdict.** None of (a) to (c) gives a power saving over one point per cell, in the toy model or in the real
problem. Each method fails at a definite, computable point (§§2 and 3). §1 gives one structural reason that covers
all of them. These methods are *Archimedean*: they use the norm only through a window of size ≍ w². At the critical
width, the union over that window genuinely has about one point per cell. In the literature (§4) I found nothing for
this configuration beyond Jarník-type and slicing bounds. Those are weaker than the S³ determinant bound n^{2/5}. The
closest new reference is Zhang–Zhu (arXiv:2606.08650, June 2026), which also gives an exact dictionary between the
toy count and eigenfunction restriction.

## 0. Normalisation

- **Toy model.** R = n^{1/2} and w = n^{1/10} = R^{1/5}, so ε = w/R = n^{−2/5}. The plane Π is a 2-plane through 0,
  and P_⊥ is the orthogonal projection onto Π^⊥.
- **Tube count.** N_Π(n) = #{x ∈ Z⁴ : |x|² = n, |P_⊥x| ≤ w}.
- **Cells.** Q ≍ R/w = n^{2/5}.
- **Volume.** E_Π N_Π(n) = r₄(n)·w²/n ≍ n^{1/5}. This is exact: for Y uniform on S³, |P_⊥Y|² is uniform on [0,1].
- **Real problem.** R' = 2^{t/4} and ε = R'^{−8/5}. The σ₁-tube has width w' = εR' = R'^{−3/5}, which is below 1.
  Cells number ε^{−1} = R'^{8/5}, and the volume is R'^{4/5}.

## 1. The norm-window barrier (main structural finding)

**Proposition 1 (toy model).** Put c₁ = 1/250 and U_Π := {x ∈ Z⁴ : |P_⊥x| ≤ w, | |x|² − n | ≤ w²}.

**(i) The determinant input holds for U_Π.** Any five points of U_Π whose Π-angles lie in an arc of length c₁w/R are
cohyperplanar, i.e. det(1, x_i) = 0.

*Proof.* Use frame coordinates (e_u, e_v, Π^⊥) at the midpoint of the arc. The four differences x_i − x₁ lie in a box
of sides 1.01c₁w (tangential), 2w and 2w (in Π^⊥), and 2w²/R (radial). The radial side uses
|P_Πx|² ∈ [n − 2w², n + w²] and cos φ ≥ 1 − c₁²w²/(8R²). Hence |det| ≤ 4!·(1.01c₁w)(2w)²(2w²/R) ≤ 200c₁·w⁵/R < 1,
since w⁵ = R. The determinant is an integer, so it vanishes. ∎

**(ii) The union genuinely fills the cells.** Average over SO(4) and use Archimedes:
E_g #U_{gΠ} = Σ_{|m−n|≤w²} r₄(m)·w²/m. For odd m, r₄(m) = 8σ(m) ≥ 8m, so the sum is at least 8w²(w² − 1) ≍ w⁴ ≍ Q.
(The true average is ≈ 2π²w⁴ ≈ Q/80 with this c₁.) Since #U_Π is bounded, #U_Π ≥ c·Q/2 on a set of planes of positive measure.

**(iii) Heuristic.** Cells of U_Π have volume ≍ c₁, so a standard Poisson/second-moment heuristic gives a positive
proportion of occupied cells. I did not complete the second moment rigorously; it needs equidistribution of E_m in
thin slabs.

**Real problem (sketch).** Put
U' = {x ∈ O : σ₁x ∈ w'-tube, |σ₁(nrd x) − σ₁(n)| ≤ w'², |σ₂x| ≤ CR'}.
- *Determinant input.* The 5-point K-determinant bounds (σ₁: ℓ³ε²R'⁴, σ₂: R'⁴) hold verbatim. The σ₁ radial spread is
  (ℓ²R'² + 2w'²)/R', and the σ₂ bound only uses |σ₂x| ≲ R'.
- *Average count.* Using r_O(m) ≍ N(m)^{1+o(1)}, the average of #U' is ≍ w'⁴R'⁴ = ε⁴R'⁸ = ε^{−1} = #cells.
- *Which norms.* U' is the union over the ≍ w'²R'² = R'^{4/5} = ε^{−1/2} totally positive m ∈ O_K with σ₁(m) ≈ σ₁(n)
  and σ₂(m) ≲ R'². Most of these m are not √2-powers.

**Congruence refinement.** Add the condition |x|² ≡ n (mod M) for some M ≤ w². Then (i) is unchanged, and the
average count becomes ≍ Q/M. In the real problem M = (√2)^j, i.e. tree depth j, and the count becomes ≍ ε^{−1}2^{−j}.

**Consequence.**
- Suppose an argument uses the tube points only through integrality, Archimedean size bounds of auxiliary
  determinants, and membership in a set insensitive to the norm within ±w² (and within congruences mod M). Then it
  cannot prove N ≤ Q^{1−κ} unless M ≥ Q^{κ}.
- Equivalently, the determinant bound "one point per cell" is the *volume* of the window-union.
- So Q1 is a statement about the tube count *not concentrating on the single norm n* among the ≍ n^{1/5} (real
  problem: ε^{−1/2}) neighbouring norms. A proof must see the exact norm arithmetically: Hecke or tree structure at
  depth ≳ κ·log(1/ε), divisor bounds in a fibration, or L-functions.
- Within-cell steps (Apps. E–H) do use the exact sphere. The barrier concerns only the *global* (cross-cell) saving.
- This is the geometric twin of NOTES §6 ("one prime"): geometric methods implicitly average over ≍ ε^{−1/2} norms,
  and the arithmetic methods that would see a single √2-power norm stop at the square-root barrier.

## 2. (a) Discrete restriction and decoupling: exponent bookkeeping

**Level-set lemma.** Let A = tube points and f_A = Σ_{ξ∈A} e(ξ·y) on T⁴. Then |f_A| ≥ |A|/2 on
{|P_Πy| ≤ c/R, |P_⊥y| ≤ c/w}, a set of measure ≍ (Rw)^{−2}. If ‖f_a‖_p ≤ K_p‖a‖₂ for all a, then

  N ≤ 4K_p² (Rw)^{4/p},  with Rw = n^{3/5}.

This is exactly Zhang–Zhu's Prop. 6.3 combined with their Thm. 2.5, rescaled to discs of radius 1/w.

Bourgain–Demeter conjecture K_p ≈ R^{max(0, 1−4/p)+ε} on S³ (IMRN 2015, Conj. 1.1, with radius N = R). Both lower
bounds are forced: a single Dirac mass gives K_p ≥ 1, and the constant coefficient vector on E_n gives
K_p ≥ cR^{1−4/p}. The resulting bounds:

| input | status | bound on N (toy) |
|---|---|---|
| p = 10/3, K = n^ε (ℓ² decoupling, BD Annals 2015) | proven | n^{18/25} = n^{0.72} |
| p > 44/7, K = R^{1−4/p+ε} (BD IMRN 2015, Thm 1.2) | proven | n^{41/55} ≈ n^{0.745} |
| p = 4 endpoint, K = n^ε | conjectural | n^{3/5} |
| any p, with K_p at its forced lower bound | best possible | ≥ n^{3/5} |
| determinant method / volume | — | n^{2/5} / n^{1/5} |

So even the sharp discrete restriction conjecture does not reach the trivial cell count.

*Cap form (S²(n), radius r = n^{3/5}).* The bound is N ≤ K_p²(r⁴/n)^{2/p}:
- p = 4 in d = 3 is proven by decoupling and gives n^{7/10}. This is exactly the Bourgain–Rudnick Jarník-type cap
  lemma.
- The p = 6 endpoint (conjectural) gives n^{7/15}.
- Both are worse than the Hopf-lifted n^{2/5}.

**Why no refinement helps.**
1. *ℓ² decoupling is blind to occupancy.* Project to Π. For fixed η = P_⊥x, the points lie on circles of radius
   √(n − |η|²), which together form an annulus of width w²/R. The canonical decoupling caps of that annulus have
   length √(R·w²/R) = w, i.e. *exactly the cells*. Decoupling is an inequality valid for every configuration, and it is
   saturated by one wave packet per cap. So it cannot certify empty caps.
2. *Restriction estimates do not see the arrangement.* The only robust lower bound for ‖f_A‖_p for A in a tube is the
   single peak above. A spread configuration (one point per cell) should have additive energy ≈ 2|A|², because sums
   ξ₁ + ξ₂ = ξ₃ + ξ₄ on S³ inside a tube are essentially trivial (heuristic, not proved). If so, it is consistent with
   every conjectured K_p.
3. *Coarser scales hit the barrier.* Decoupling at any scale coarser than the single sphere sees only the window-union
   U_Π of §1.
4. *Other decoupling families do not apply.* Small-cap decoupling (Demeter–Guth–Wang) concerns full lattice
   paraboloids, not single-norm sets. Vinogradov-type decoupling (BDG) needs a non-degenerate curve, but the core is a
   *planar* circle in R⁴.

**Exact dictionary with eigenfunction restriction (toy).**
- *Statement.* Adapt Zhang–Zhu Thm. 2.5 (Schur test plus a matching lower bound) to flat 2-discs D ⊂ T⁴ of radius
  ≍ 1/w in direction Π^⊥:
  sup_{e_λ} ‖e_λ‖²_{L²(D)}/‖e_λ‖² ≍ w^{−2}·max_{m∈E_n, Π} #{x ∈ E_n : dist(x, m+Π) ≤ cw}, with λ = R.
- *Our tubes are in the family.* Great-circle tubes containing a lattice point are of this form after doubling w.
- *Translation of the bounds.* With D of radius λ^{−1/5}:
  - the volume law corresponds to ‖e‖²_{L²(D)} ≲ λ^ε;
  - Q1 corresponds to ≲ λ^{2/5−2κ};
  - the determinant bound corresponds to λ^{2/5}.
- *What the restriction literature gives.* The best bound there, for unit-size totally geodesic 2-dim Σ ⊂ T⁴, is
  λ^{3/4+ε} (Zhang–Zhu Thm. 1.3, α(2,4) = 3/8). The general bound is λ log λ (Burq–Gérard–Tzvetkov).
- *Upshot.* In this problem harmonic analysis takes its input *from* lattice-point counts, and is currently behind
  the determinant bound.

**Real problem.** The barrier of §1 applies verbatim. The level-set bookkeeping is also worse here:
- The σ₁ alignment in the P_⊥ directions needs |P_⊥y₁| ≤ c/w' = cR'^{3/5}, which is wider than the torus, so it
  constrains nothing. The aligned set therefore has measure ≍ R'^{−2}·R'^{−4} = R'^{−6}.
- The forced lower bound for this 8-dimensional E is K_p ≥ max(1, cR'^{2−8/p}).
- Together these give N ≤ 4K_p²R'^{12/p}, which is at least R'³ = ε^{−15/8} for every p. The cell count is
  ε^{−1} = R'^{8/5}.

The free second embedding destroys the Knapp localisation altogether.

## 3. (b) and (c) Polynomial method, polynomial partitioning, and Huxley-type arguments

**Degree-D auxiliary forms, leading order (toy).**
- *Setup.* Take boxes of core length L = R^λ and width w = R^ω. Use local coordinates u (size L), s and s'
  (size w), and h (radial, size H = R^{2max(λ,ω)−1}).
- *Monomials.* On the exact quadric, monomials are u^a s^b s'^c h^f with f ≤ 1, and there are
  C(D+3,3) + C(D+2,3) of them.
- *Log-size of the determinant.* It is (C(D+3,4) + C(D+2,4))(λ + 2ω) + C(D+2,3)(2max(λ,ω) − 1), up to constants.
- *Vanishing criterion for λ ≤ ω.* The determinant is forced to vanish iff λ < 2(1 − 2ω)/(D+1) − 2ω.
- *At the critical width ω = 1/5* this reads λ < λ_D = (4 − 2D)/(5(D+1)):

| D | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| λ_D | 1/5 (one cell) | 0 | −1/10 | −4/25 |

- *No quadric reduction.* Using all C(D+4,4) monomials, which also covers U_Π, the criterion is 3λ + 2ω < 1 for
  every D. That is the same box as D = 1.

So D = 1 is optimal, and every D ≥ 2 forces vanishing in no box of length ≥ 1.

**Real problem at the critical scale.**
- *One-cell boxes.* Take λ = log_{R'}L = −3/5. The log-norm of the degree-D determinant is
  (6/5)·C(D+2,3)·(D − 1)/2, which is 0 for D = 1 (borderline: exactly one cell) and 12/5, 12, 36 for D = 2, 3, 4.
- *Sub-cell boxes.* Higher D forces vanishing only for L < R'^{12/(5(D+1)) − 9/5}, i.e. R'^{−1} or R'^{−6/5} for
  D = 2, 3. These boxes are strictly smaller than a cell (R'^{−3/5}).
- *Comparison.* This confirms NOTES §3(d) with an exact formula.

**Reason.** The only small coordinate is the normal one (radial), and modulo the quadric it enters at most linearly.
The free directions have size at least 1 in norm: the tube cross-section in the toy model, and the σ₂ sphere in the
real problem. So every higher monomial costs more than it gains.

**Polynomial partitioning (Guth–Katz).** There is no incidence structure to exploit.
- Partitioning a near-1-dimensional set by Z(P) yields cells along the core, which are just boxes. The algebraic part
  (points on Z(P)) is the determinant method above.
- Incidence bounds for points against Hecke-translated tubes reduce to the multiplicative energy, i.e. NOTES §3(b),
  which gives rate 2.
- Bourgain–Demeter's incidence-theory L⁴ argument bounds energies, not single tubes.

**Huxley, Swinnerton-Dyer and Bombieri–Pila.**
- *The planar model case.* For points within δ of a curve of length N with |f''| ≍ 1/N, the 3-point determinant gives
  N^{2/3} + Nδ. At δ_c = N^{−1/3} this is the volume of the δ_c-neighbourhood, so it is sharp there: the same barrier
  as §1.
- *Higher derivatives do not help for circles.* Huxley–Sargos (k ≥ 3) gives
  R ≪ Nλ_k^{2/(k(k+1))} + Nδ^{2/(k(k−1))} + (δ/λ_k)^{1/k} + 1. This improves only the δ-term; for circle-like curves
  Nλ_k^{2/(k(k+1))} ≥ N^{2/3}.
- *Where the genuine gains come from.* Bombieri–Pila (N^{1/2+ε}) and Swinnerton-Dyer (N^{3/5+ε}) need two things: the
  points lie a power *closer* than δ_c, and the set is one-dimensional at the lattice scale, so that higher-degree
  interpolation is cheap.
- *What our problem satisfies.* The exact norm does put the points a power below the critical radial window
  (0 versus w²/R). But at the critical width the set is *not* one-dimensional: the cross-section w equals the cell
  length. So degree raising loses (table above).
- *A further obstruction.* The quadric is exactly quadratic, so there is no third-derivative (Swinnerton-Dyer) step.

**Space curves and rational points near manifolds (Huang, Beresnevich–Vaughan–Velani, Beresnevich–Yang).**
- *Degenerate core.* These results need a non-degenerate curve or manifold (non-vanishing torsion, curvature). The
  core is a planar circle in R⁴, a geodesic of S³, and a *line* in P³. Near a flat core, counts are governed by the
  Diophantine type of Π (Schmidt), not by curvature.
- *Wrong averaging.* These results also sum over all denominators q ≤ Q. Our points are x/√n with a single fixed
  "denominator", which is exactly the single-norm problem.
- *Conclusion.* No Huxley-type improvement exists over the first-order determinant at this width, and §1 explains
  why none can exist within Archimedean methods.

## 4. (d) Literature: what is known for these tubes and caps

Each item lists the source and what it gives at our scale.

1. **Bourgain–Rudnick (GAFA 22 (2012)), Lemmas 2.1 and 2.2.** For a cap of radius λ on the sphere of radius R in Z³:
   R^ε(1 + λ²/R^{1/2}) (Jarník) and R^ε(1 + λ) (slicing). For S²(n) at r = n^{3/5} these give n^{7/10} and
   **n^{3/5}**.
   - *Correction to an earlier draft and to NOTES.* The general BR bound at this scale is n^{3/5+o(1)}, via Lemma 2.2,
     not n^{7/10}. It is still worse than n^{2/5}.
   - *On NOTES §8(d).* Monotonicity in r gives the Hopf-lift bound n^{2/5+o(1)} for *all* r ≤ n^{3/5}, which is
     stronger than n/r there. It beats min(r²n^{−1/2}, r) exactly for r ∈ (n^{9/20}, n^{3/5}].
   - *Prior art.* I found no prior statement of this lift for square radii.
2. **Slicing by Diophantine vectors.**
   - *Sources.* Maffucci, JFA 2017 (arXiv:1611.00571, Props. 6.2–6.3), and Ortiz Ramírez, Monatsh. 2020
     (arXiv:1912.10462). Both treat codimension-1 bands: ψ ≪ κ_d(R)(1 + Rθ^{1/d}).
   - *Adapted to our tube.* For generic Π (all successive minima ≍ 1), take two Minkowski vectors b₁, b₂ with |b_i| ≤ (R/w)^{1/2} and |P_Πb_i| ≲ (w/R)^{1/2}. The
     points then lie on ≲ Rw circles, each with ≤ n^{o(1)} points, giving N ≤ n^{3/5+o(1)}. This is the same exponent
     as the restriction endpoint.
3. **Discrete restriction.** Bourgain–Demeter, Annals 182 (2015) (decoupling; p ≤ 10/3 on S³), and IMRN 2015
   (arXiv:1310.5244; p > 44/7, endpoint p = 4 open). See §2: at best n^{3/5}.
4. **Huang–Zhang (Anal. PDE 14 (2021)).** Sharp restriction to *rational* totally geodesic submanifolds. This is the
   analogue of the Clifford-frame volume law (C).
5. **Zhang–Zhu (arXiv:2606.08650, 2026).** This is the closest reference.
   - *Their quantity.* A_{2,4,λ} is the maximal number of points of λS³ in the 1-neighbourhood of a 2-plane, i.e. a
     unit-width tube around any circle. They prove A_{2,4,λ} ≲ λ^{3/4+ε} (Prop. 3.9) and conjecture λ^ε
     (Conj. 2.4).
   - *Their method.* Successive minima and "slicing and packing", plus the Magyar–Stein–Wainger multiplier.
   - *At our width.* Slicing the critical tube into w² unit tubes gives N ≲ w²λ^{3/4} = n^{23/40}, which is worse than
     n^{2/5}.
   - *Their conjecture would settle the toy questions.* It gives N_Π(w) ≲ w²n^ε for all w ≥ 1 and all planes, i.e.
     the toy volume law and hence toy Q1 and Q2. So the toy questions are special cases of an open conjecture in this
     literature, whose best known bound is weaker than ours.
6. **Average and conditional results** (already in NOTES): Humphries–Radziwiłł (arXiv:1910.01360), Lutsko
   (arXiv:2402.12822), Shubin under GRH (arXiv:2108.00726), and Bourgain–Rudnick–Sarnak (Bull. Iranian Math. Soc.
   2017, arXiv:1606.05880; random caps only). None of these bounds individual caps.
7. **Space curves and rational points near manifolds.** Huang (Math. Ann. 374 (2019), arXiv:1809.07796; needs
   non-vanishing torsion) and Hickman–Srivastava (IMRN 2025). Inapplicable, as explained in §3.

**Bottom line for (d).** I found no individual bound for lattice points of S³(√n) near a great circle, or for
S²(n) caps (square radius) at r = n^{3/5}, beyond Jarník and slicing (n^{3/5}). The S³ five-point determinant
(n^{2/5}) is already ahead of the published bounds, Zhang–Zhu 2026 included.

## 5. What a proof of Q1 would have to use

By §1, it must distinguish the exact norm n among ≳ Q^{κ} of the ≍ n^{1/5} (real problem: ε^{−1/2}) norms in its
window. It must do so by information other than Archimedean sizes and congruences to moduli below Q^{κ}.
- *Toy model.* The available inputs of that kind are: unique factorisation in the Hecke tree at depth ≳ κ log Q
  (NOTES §3(e): loses), divisor bounds in a Dirichlet fibration (NOTES §3(a): PV-type completion, gives 5/2), and
  L-functions (NOTES §6: no averaging over norms at one prime).
- *Real problem.* Decoupling, restriction and polynomial methods contribute nothing beyond what is already proved.
- *Possible spin-off (unverified).* The 5-point determinant at width w = 1 gives boxes of length R^{1/3}. If the
  within-box count were R^{o(1)}, this would improve Zhang–Zhu's A_{2,4,λ} ≲ λ^{3/4} to λ^{2/3}. A naive within-box
  estimate via codimension-1 band bounds (λ^{1/12}) returns exactly λ^{3/4}, so this is not claimed.
