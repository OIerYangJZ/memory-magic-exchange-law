# Round 9 scratch: five frames in a common cap, different fibre arcs (brief item (i))

## Dictionary — T-count units ONLY (the round-5 slip came from mixing sde k and T-count t = 2k)
N = #states of T-count t ≍ 2^t.  n := |u|²+|t|² ≍ 2^{t/2} (numerator norm), R' = √n ≍ 2^{t/4} (App. E), Bloch radius R = 2^{t/2}.
Critical scale ε = 2^{-2t/5} = N^{-2/5}.  vol = ε²N = 2^{t/5}.  1/ε = 2^{2t/5} = vol².  B = εN = 2^{3t/5} = vol³.  √N = 2^{t/2} = vol^{5/2}.
Cells per band 1/ε; cells on S² 1/ε² = 2^{4t/5}.  Fibre ε-arcs over a cap: 1/ε = 2^{2t/5}.  Determinant level N_c ≤ 1/ε = vol².
Spike (model M_vol): M = 1/ε = 2^{2t/5} frames, one per arc.  M² = 2^{4t/5}.  Quotient tube at T-count ≤ 2t: 2^{o}(1+ε²4^t) = 2^{6t/5+o}.
⇒ §3(b) check: M² ≤ 2^{6t/5} holds with room; quotient counting gives only N_c ≤ 2^{3t/5} = B (rate 2). Consistent with the model.

## Lemma E1 (App. E) for five frames in one cap
Tube point x = R'(cos φ·C(s) + sin φ·n), φ ≤ 2ε (transverse offset ≤ 2εR'), core parameter s. Five frames in a common cap with
core parameters in an arc of length λ: coordinate variations e₁ ≤ Cλ²R' (+CεR'·ε), e₂ ≤ λR', Π^⊥: ≤ 4εR' each.
  |det(w_i − w_0)|_{σ₁} ≤ C λ³ ε² R'^4,   |det|_{σ₂} ≤ C R'^4,   det ∈ (1/16)Z[√2].
Vanishing iff C² λ³ ε² R'^8 < 2^{-8}, R'^8 ≍ 2^{2t}:  λ³ < c ε^{-2} 2^{-2t} = c 2^{4t/5 − 2t} = c 2^{-6t/5}  ⟺  λ < c' 2^{-2t/5} = c' ε.
So the K-determinant of five cap frames vanishes exactly when they lie in ONE ε-arc (up to the constant c'). For frames X arcs
apart (λ = Xε): |det|_{σ₁}|det|_{σ₂} ≤ C X³ ε⁵ 2^{t}·2^{t} = C X³ — nonzero allowed as soon as X ≥ c. No relation across arcs.
Full fibre (λ ≍ 1): vanishing iff ε < c 2^{-t}: the Liouville scale, useless.

## Other algebraic relations for frames in a common cap (all checked, none vanishes at ε = 2^{-2t/5})
- Bloch 4-point determinant det(1, ν_i), ν_i ∈ O_K³, |ν|² = 2^t: σ₁ ≤ Cε⁴2^{3t/2}, σ₂ ≤ C2^{3t/2}; vanishes iff ε < c2^{-3t/4}. No.
- 2×2 minors g = ut' − tu' ∈ Z[ζ₈]: |σ₁g| ≤ 2ε2^{t/2}, |σ₂g| ≤ 2^{t/2+1}; N_{Q(ζ₈)/Q}(g) ≤ 16ε²2^{2t} = 2^{6t/5+4} ≥ 1. Vanishes iff ε < c2^{-t}. No.
- Complex 3-point Δ over Z[ω] (§3(f)): vanishes only for ℓ < 1/(8εR'^4): weaker than E1. No.
- Gram determinant of the five numerators = (K-determinant)²: same condition as E1.
- Lagrange identity for quotients: |ūu'+t̄t'|² + |ut'−tu'|² = n²: an identity, no constraint beyond the minors.

## What this leaves for Σ_c N_c⁵ (5-point criterion ⇒ 13/5)
Target: Σ_c N_c⁵ ≤ 2^{o}·(1/ε²)·vol⁵ = 2^{o}2^{4t/5 + t} = 2^{9t/5+o}.
Write N_c = Σ_{arcs a} n_{c,a}, n_{c,a} ≤ 2^{o} (E1 + E5 within one arc). Same-arc 5-tuples: Σ_c Σ_a n_{c,a}⁵ ≤ 2^{o}·ε^{-2}·ε^{-1} = 2^{6t/5+o}: harmless.
Cross-arc 5-tuples (five frames in the same cap, five different arcs): no algebraic relation (above); they are exactly the 5-tuples
of the spike model (M⁵ = 2^{2t} > 2^{9t/5}). So the 5-point criterion is a statement purely about cross-arc tuples, on which the
determinant method is silent at every scale down to the Liouville scale 2^{-t}.
Mixed tuples (two arcs, etc.): the relation for "k frames in one arc + others elsewhere" is still E1 on the ≥5 frames of one arc only.

## Verdict (round 9)
Item (i) is negative with exact scales: the K-determinant (the whole engine of 5/2) constrains frames only within a single fibre
ε-arc (λ ≤ c'ε); across arcs the first scale where it bites is the Liouville scale. Σ_c N_c⁵ is governed by cross-arc tuples.
