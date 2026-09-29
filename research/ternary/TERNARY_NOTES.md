# Small-cap / thin-rod equidistribution on ternary Z[√2]-spheres (2026-09-29)

## What is needed
Worst case of level 2 (App. E): a small sphere section Y = S³ ∩ H, with σ₁-radius r ≈ R'^{−0.1} and
σ₂-radius r₂ ≈ R'. The core meets Y tangentially, and the rim patch is a thin rod.

Balance the two embeddings by a unit of Z[√2]. Then:
- Y carries P ≍ r r₂ ≈ R'^{0.9} points.
- The rod has relative length P^{−0.4} and relative width P^{−0.8}, i.e. area P^{−0.2} × (mean spacing)².

The needed statement, **Rod hypothesis (NC)**: the points of X' on a sphere section inside such a rod lie on
R'^{o(1)} circles.

**Consequence (conditional).** Under NC, level 2 costs nothing, and the level-1 box count alone gives E(a₀) = (8−2a₀)/3,
so α ≥ 5/2 − o(1). Check at α = 5/2: a₀ = a = b = 1.6, so 2a+3b = 8, a ≥ b, and a₀+a = 2b; all constraints are tight, so
any α < 5/2 works with strict margins.

## Why it is beyond present methods
- The only elementary tool is 4-point coplanarity (Jarník / Bourgain–Rudnick), our Lemma E3. At this scale it loses
  exactly P^{1/3}. That is the observed R'^{0.3}, since 0.3/0.9 = 1/3.
- **The same P^{1/3} loss occurs over Q at the same relative scale.** For x²+y²+z² = n (P ≍ R = √n), a rod of
  relative size P^{−0.4} × P^{−0.8} loses R^{1/3} under the plane method. So the second embedding is not the
  obstruction. The obstruction is the plane method itself.
- Over Q the state of the art is Bourgain–Rudnick's cap lemma F₃(R,r) ≪ R^ε(1 + r²/R^{1/2}), against the expected
  r²/R. This is open, and it is stated as a conjecture of Bourgain–Rudnick; Humphries–Radziwiłł reach optimal
  equidistribution only for caps of relative radius ≫ n^{−1/24}.
- Our rod version is stronger than the cap conjecture. A cap of the rod's length would still lose P^{0.2}. What is
  needed is the "thin-region" version (count ≪ P^ε(1 + P·area)).
- Automorphic methods (Waldspurger / subconvexity for the relevant Hilbert modular forms) give equidistribution
  only at macroscopic or slowly shrinking scales. The needed scale is below the mean spacing.
- Tried and failed:
  - the height of H;
  - the σ₂ radius;
  - 3-point (wedge) and 4-point determinants;
  - quaternion products ξ'ξ of nearby points, which move the problem to caps around ±N but square the norm and
    make σ₂ worse;
  - the (a,b) ∈ Z³×Z³ form, which gives two quadrics in 6 variables, too few for a circle method.

## Numerical evidence for NC (`reps.py`, `cluster.py`, `rods.py`)
All representations of N as a sum of three squares in Z[√2], with σ₁(N) small and σ₂(N) large (4 instances, up to
1920 points):

- **No sub-spacing clusters:**
  - no ball of radius 0.5× the mean spacing contains 4 points (max 3);
  - at 0.3× the maximum is 2;
  - at ≤ 0.1× the maximum is 1.
- **Thin rods:** rods of length 2–8× the spacing and width 0.02–0.1× the spacing contain at most 2–4 points. This is
  Poisson-like, with mild repulsion.

## Status
- α ≥ 2.23 is rigorous.
- α ≥ 5/2 holds conditionally on NC, the thin-rod analogue for ternary Z[√2]-spheres of the Bourgain–Rudnick cap
  conjecture. NC is supported numerically.
- α = 3 needs Conjecture H, i.e. also removing the level-1 limitation.
