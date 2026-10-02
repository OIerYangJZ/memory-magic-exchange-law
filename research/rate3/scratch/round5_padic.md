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
Split V = PQ, P ∈ Λ_p, Q ∈ Λ_{2k−p} (§10 normalisation: sde k, word length 2k, ε = 2^{-4k/5}).
n_c(P) = #{Q : Q|0⟩ ∈ P^{-1}c}: a cap of radius ε for the depth-(2k−p) orbit, N'' = 2^{2k−p} points.
App. E determinant scale for N'' points: ε_det'' = N''^{-3/5}. ε ≤ ε_det'' ⟺ 4k/5 ≥ 3(2k−p)/5 ⟺ p ≥ 2k/3.
Sanity: p = 0 gives (ε/ε_det)² = 2^{-8k/5+12k/5} = 1/ε ✓ (the "one point per cell" bound N_c ≤ 1/ε).
So for p ≥ 2k/3 (prefix at least one third of the word) every prefix class contributes ≤ 2^{o(k)} points to the cell.
H_p := #{P ∈ Λ_p : n_c(P) ≥ 1}. Then N_c ≤ 2^{o(k)} H_p (p ≥ 2k/3); H_p ≤ H_{p+1} ≤ 2H_p; H_{2k} = N_c.
SA ⟺ H_{2k/3} ≤ vol·2^{o(k)}: the cell's points do not branch below depth 2k/3 in the class tree.
U := ε-neighbourhood of the fixed orbit Λ_{4k/3}|0⟩, measure 2^{4k/3}ε² = 2^{-4k/15} = N^{-2/15}.
H_{2k/3} = #{P ∈ Λ_{2k/3} : P^{-1}z₀ ∈ U}. Expected 2^{2k/3}·2^{-4k/15} = 2^{2k/5} = vol ✓.
Ramanujan: error 2^{k/3}·‖1_U‖₂·(#harmonics ≤ 1/ε) ≈ 2^{k/3}·2^{-2k/15}·2^{4k/5} = 2^k = √N. Same wall.
This is §10(D) meet-in-the-middle with the asymmetric split 1:2 chosen so that the determinant method kills multiplicity.

## 5.5 Deviation fields and the Hecke recursion (exact)
F_j(z) := N_{cap(z,ε)}(j) − ε²·#states(j), fixed ε. 3-regular tree: A_2A_{2j} = A_{2j+2} + A_{2j} + 4A_{2j−2}, hence
T_2 F_k = F_{k+1} + F_k + 4F_{k−1},  T_2 := Σ_{g∈Λ_2} g^* (6 terms). Ramanujan ⇒ spec(T_2 | L²₀) ⊂ [−3, 5].
Per eigencomponent F_k ∝ r^k, r² = (λ−1)r − 4; |λ−1| ≤ 4 ⇒ |r| = 2 exactly: Poisson growth 2^k saturates Ramanujan.
Backward propagation of a spike X·vol at level k: some translate has ≥ (2X/3)vol_{k−1}; after s steps ≥ X vol_k 6^{-s};
determinant bound at level k−s, radius ε: 2^{o}·max(1, 2^{4k/5−12s/5}); contradiction needs X > vol·(6/2^{12/5})^s = vol·1.137^s.
Best s = 0: nothing beyond the determinant bound. Forward propagation: §14(c).
Fixed z, k ↦ F_k(z)/2^k = Σ_f c_f(z)e^{2ikθ_f}-type almost periodic sequence, ≈ ε^{-2} Satake angles; max/RMS up to 1/ε allowed
unless the angles are spread mod π/k, i.e. at scale 1/log(family): vertical Sato–Tate (discrepancy ≍ 1/log) is silent there.

## 5.6 Literature (WebSearch only; arXiv, TAU, Wiley, UZH blocked by the proxy)
Bourgain–Rudnick: F₃(R,λ) ≪ R^ε(1+λ²) for caps of radius λ on x²+y²+z² = R² (plane/determinant method); nothing beyond it.
Humphries–Radziwiłł (CPA 2022): variance over random caps/annuli (almost all caps), parts under GRH. Not individual caps.
Burrin–Gröbner 2502.17678: rational points of height ≤ T (average over n), Hecke operators at varying primes.
Zhang–Zhu tube conjecture: not found. Verdict: the rational analogue is also stuck at the determinant method for individual caps.
