# Pushing beyond α ≈ 2.22 (working notes, 2026-09-28)

Baseline: Result E in ../NOTES.md §8, the plane-only chain hyperplane → sphere → circle, gives τ* ≈ 2.225L, i.e. α ≈ 2.22.

## Where the exponent is lost
The bits carried at cost τ equal (τ/4)·[E₁ + E₂]:
- E₁ = 2(a'−a₀) + b' is the level-1 box count. It is limited by the level-1 condition, which is 2a'+3b' > 8 at d = 1.
- E₂ is the worst-case level-2 overhead. It comes from strongly curved pieces: round-sphere sections of radius ≈ R'^{−0.3…−0.4} met near their rim, which cost about R'^{0.3}.

Ceilings (all exponents exact):
| level-1 degree | level-2 overhead | τ*/L | α |
|---|---|---|---|
| d = 1 | as proved (worst round sphere) | 2.225 | **2.22 (Result E)** |
| d = 1 | none (E₂ = 0) | 2.5 | 2.5 |
| d → ∞ | none (E₂ = 0) | 8/3 | 2.67 |
| d → ∞ | model: auxiliary surface ≈ worst round sphere | 2.31 | ≈ 2.31 (model, unproved) |

## What was tried
1. **High-degree determinants at level 2 on the rigid sphere, followed by an honest level 3** (`level3_bezout.py`). The points of a region lie on a curve Γ of degree ≤ 2d₂. Γ is then cut into pieces of diameter z. By Bézout, a plane meets Γ in ≤ 2 deg points unless the circle lies on Γ, in which case the divisor bound applies. Result: **no gain.** Level 3 needs isotropic pieces z ≈ R'^{1−(1.5+g/4)}, which costs as much as it saves. The earlier 2.24 had assumed level 3 was free.
2. **Level 1 of degree d**, exact box conditions (`level1_degree.py`):
   - d=1: 2a'+3b'>8
   - d=2: 12a'+16b'>44
   - d=3: 40a'+50b'>140
   - in general, → a'+b'>3.

   Even if a degree-d auxiliary surface behaved no worse than the worst round sphere at level 2, the rate would only reach ≈ 2.28 (d=2,4) and ≈ 2.31 (d=40). For d ≥ 2 this is a model, not a proof. The obstruction is surfaces with ridges or bumps along the core, which act like small spheres. There is a heuristic reason bounded degree limits this: a great-circle arc meets the surface in ≤ 2d points, so there are only O(1) bends per box. Turning it into a rigorous "sagitta lemma for bounded-degree surfaces in thin boxes" is open. The dangerous examples are torus-like quadric sections around the core. For rational core planes these are arithmetically harmless: points on them split into two binary problems, so the divisor bound applies. A rational quadric section can also only approximate an irrational torus.
3. **Localisation in σ₂.** Not useful. At level 1, one unit of σ₂-cap exponent gains 5 at a cost of 3, while the along direction gains 3 at a cost of 1. At level 2 it gains 4 at a cost of 2, while the long side gains 3 at a cost of 1.
4. **Level 3 by torsion/curvature of curves on spheres** (sketch only). A curve piece of length z has 4-point volume ≲ z⁶κ²τ_tors. With κ ≈ 1/r_Y and averaged torsion ≈ 1/x_L this saves ≈ 0.08 in the exponent, and it needs an unproved total-torsion lemma. Not pursued.

## Conclusion
The rigorous chain gives α ≈ 2.22. The next rigorous gain would need a genuinely new treatment of strongly curved level-2 pieces, i.e. small sphere sections met tangentially by the core. Without that, the ceiling of this method family is 2.5 (d = 1) or 8/3 (d → ∞), and realistic variants stop around 2.3.
