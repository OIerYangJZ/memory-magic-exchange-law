# Round 5 scratch: P-adic Fourier form of the class count; depth-split determinant bound

Normalisation (§15(c)): states of T-count 2τ ↔ ν ∈ O_K³, |ν|² = 4^τ = P^{4τ}·unit... (P = (√2), N(P)=2, 4^τ = P^{4τ}).
Cap radius ε = 2^{-2τ/5}; N_cap = ε²·4^τ = 2^{6τ/5}; prefix classes P ∈ Λ_τ (≈ 2^τ); vol = N_cap/2^τ = 2^{τ/5}.
Class P ⟺ ν ∧ ν_P ≡ 0 (mod P^τ) ⟺ ν ∈ L_P := O_K ν_P + P^τ O_K³, a lattice of index 2^{2τ} in O_K³.

## 5.1 Additive-character expansion (exact)
1_{L_P}(ν) = 2^{-2τ} Σ_{η ∈ L_P^⊥} e_{P^τ}(⟨η,ν⟩),  L_P^⊥ = {η mod P^τ : ⟨η, ν_P⟩ ≡ 0 (P^τ)}, |L_P^⊥| = 2^{2τ}.
n(P) = 2^{-2τ} Σ_{η∈L_P^⊥} Ŵ(η),   Ŵ(η) := Σ_{ν ∈ X_{4^τ} ∩ cap} e_{P^τ}(⟨η,ν⟩).
η = 0 gives vol. So  n(P) − vol = 2^{-2τ} Σ_{η ≠ 0} Ŵ(η).
Sizes: trivial |Ŵ| ≤ N_cap ⇒ 2^{16τ/5}·2^{-2τ} = 2^{6τ/5} (= N_cap, trivial). Square-root per η with absolute values:
2^{-2τ}·2^{2τ}·2^{3τ/5} = 2^{3τ/5} = vol³ ≫ vol. So cancellation over η is also needed (random signs give 2^{-2τ/5}).
Grouping η by P-adic valuation j: Σ_{η ∈ P^j L_P^⊥} Ŵ(η) = 2^{2(τ−j)} n_{τ−j}(P|_{τ−j}): the Hecke recursion by depth.

## 5.2 Parseval over L_P^⊥
Σ_{η∈L_P^⊥}|Ŵ(η)|² = 2^{2τ}·#{(ν,ν') ∈ cap² : ν − ν' ∈ L_P}  ≥ 2^{2τ}(N_cap + n(P)²).
Diagonal N_cap = 2^{6τ/5} already exceeds the target vol² = 2^{2τ/5}. Dead (L² over a family loses the diagonal).
Random model check: N_cap²/2^{2τ} = 2^{2τ/5} = vol² off-diagonal, consistent with Poisson.

## 5.3 What Ŵ(η) is
Σ_ν e_{P^τ}(⟨η,ν⟩) 1_cap(ν) over X_{4^τ} = coefficient at n = 4^τ of a theta series with characteristic η/P^τ (a form of
level ≍ P^{2τ}) with harmonic weights for the cap. Coefficient index = level: depth aspect for the quaternion forms.
Spectrally this is again Σ_f c_f(η) λ_f(P^{2τ}), i.e. the same object as §14(b). No new input.

## 5.4 Depth-split determinant bound (rigorous)
Split V = PQ, P ∈ Λ_p, Q ∈ Λ_{2τ−p}. n_c(P) = #{Q : Q|0⟩ ∈ P^{-1}c}: a cap of radius ε for the depth-(2τ−p) orbit
(N'' = 2^{2τ−p} points). App. E determinant scale for N'' points: ε_det'' = N''^{-3/5}. ε ≤ ε_det'' ⟺ p ≥ 2τ/3
(critical scale ε = 2^{-2τ/5}: 2τ/5 ≥ 3(2τ−p)/5 ⟺ p ≥ 2τ/3... check: 2τ ≥ 6τ − 3p ⟺ p ≥ 4τ/3?? redo below)
