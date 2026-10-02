# Round 10 scratch: true properties of the orbit tested against the spike model M_vol (T-count units, see round9 dictionary)
Spike: M = 1/ε = 2^{2t/5} frames over one cap c of radius ε = 2^{-2t/5}, one per fibre ε-arc, at the single norm P^τ (T-count t).

(a) Exact multiplicativity X_{P^τ l} = X_{P^τ}·X_l/48, (l,P)=1, N(l) = Λ.
    Induced frames w·v (v ∈ X_l) at norm P^τl lie in the ≈ Λ caps c·v, ≥ M each (48-to-1). Volume per cap at norm P^τl: ε²·2^t·Λ^{1+o} = vol·Λ.
    Induced spikes are above volume iff Λ < vol. Known statements at composite norms: (1) norm-averaged volume law (§15(a)):
    the spikes at norms P^τl, Λ ≈ Λ₀, contribute Λ₀²M to a sum with main 2^τΛ₀²vol; ratio vol/2^τ ≪ 1, below the Poisson fluctuation
    (2^τΛ₀²vol)^{1/2} of the other norms since vol³ ≪ 2^τ. Consistent. (2) Hecke-orbit average over v ∈ X_l of N_{P^τ}(c'v^{-1}) equals
    48·vol·Λ^{1+o}(1+o(1)) for Λ ≥ ε^{-2-δ} (spectral, §15(a)); one translate carrying vol² against a main term ≥ vol·ε^{-2} = vol⁵: invisible.
    ⇒ absorbed.
(b) Σ_W D(W) = whole-sphere ε-pair deviation, main term known, error ≤ k²N = 2^t k² (Ramanujan). Spike: M² = 2^{4t/5} < 2^t. Absorbed.
(c) Unique factorisation of quotients V^{-1}V' (round 9): M² = 2^{4t/5} ≤ 2^{6t/5+o} (Theorem clifford at T-count ≤ 2t). Absorbed.
(d) Galois symmetry: the conjugate configuration has σ₁-positions = σ₂-shadows of the spike, which are unconstrained below √N. Absorbed.
(e) Flat tube bound at every level + Hecke recursion: forward, 3·2^{s−1} translates with ≥ M frames at T-count t+s, bound 2^{2(t+s)/5} ≥ M ✓;
    pair deviation at level t+s ≈ 1.5·2^s M² = 1.5·2^{s+4t/5} ≤ k²2^{t+s} ✓. Backward (round 5): ≤ 1.137^s loss, no contradiction. Absorbed.
(f) Exact band identity B_W = Σ_{m∈W} r(m)r(n−m) ≤ 2^{o}·2^{3t/5} and ≤ 2^{o} frames per norm m: the spike uses ≥ M 2^{-o} distinct m ∈ W,
    |W| ≈ 2^{3t/5} ≫ M. Absorbed. (No asymptotic for B_W is known at this scale anyway, §16(e).)
(g) Clifford-coset cube count (words of T-count ≤ 2t within ε of a fixed R_z(ψ)): spike pairs with angle difference ≈ ψ number ≈ M,
    cube volume ε³4^t = 2^{4t/5} ≥ M; L3 allows 2^{6t/5}. Absorbed even if a cube volume law were proved.
(h) Window-union count (agent_decoupling): the union over the ≈ ε^{-1/2} = vol norms in |nrd − n| ≲ w² has ≍ 1/ε points per cap on average
    (and ≤ 2^{o}/ε by the determinant bound). The spike puts 1/ε at one norm; the other norms then carry a deficit in c, of size ≪ their
    total. No per-cap lower bound is known for the union. Absorbed.
(i) Prefix-class structure + L1: classes at depth p carry ≤ 2^{2(t−p)/5+o} each (flat bound) and 2^p·2^{2(t−p)/5} ≥ M for all p ✓;
    L1 bounds shared prefixes of frames at distance Xε by t/5 + 2log₂X: no constraint on the distribution. Absorbed.
(j) Known arithmetic excesses (rich circles at |+⟩, H±, §11): divisor-type, ≤ 2^{o}, at arithmetic points only. Irrelevant to a generic cap.
(k) Higher Hecke relations (T^{(s)}), spectral positivity, Weyl-sum positivity: all consequences of (R)+(e). Absorbed.

Verdict: every true property listed is satisfied by M_vol. No provable property of the orbit distinguishes the spike.
Structural reason: each candidate is either (α) a Weyl-sum statement with error ≥ √N ≫ M², (β) an exact identity for a divisor-type
total (band, union, Hecke average) in which the spike is a lower-order part, or (γ) a counting bound (tube, cube, quotient) with capacity
≥ the spike's demand. A killer must be a statement about a *single* norm, a *single* cap and *sub-√N* accuracy: exactly SA-type.
