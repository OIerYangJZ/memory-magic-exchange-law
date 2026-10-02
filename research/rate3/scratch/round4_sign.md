# Round 4 scratch: exploiting the sign in NOTES §15(d)

## Normalisations (point form, §10(A); k = sde, N = 4^k states, ε = 2^{-4k/5} critical)
- band B = ε·N = 2^{6k/5};  cells per band 1/ε = 2^{4k/5};  vol = εB = 2^{2k/5}.
- Useful identities: 1/ε = vol^2,  B = vol^3,  N = vol^5,  √N = 2^k = vol^{5/2}.
- Determinant sup bound: N_c ≤ 1/ε = vol^2.  Ramanujan sup bound: N_c ≤ vol + √N = vol^{5/2}.
- S_ℓ(W) = Σ_{x ∈ band W} e^{iℓα_x} (azimuth Weyl sum), |ℓ| ≤ 1/ε, S_0 = B.
- D(W) := P(R) − B·vol = Σ_c (N_c − vol)^2 = Σ_{ℓ≠0} ĉ_ℓ |S_ℓ|^2,  ĉ_ℓ ≈ ε,  Σ_ℓ ĉ_ℓ ≈ 1.
- Target (⇒ α > 5/2): D(W) ≤ B·vol·2^{-δk}  ⟺  Σ_{0<|ℓ|≤1/ε} |S_ℓ|^2 ≤ B^{2−δ}.
- Rich cell N_c = X·vol  ⇒  (Cauchy–Schwarz in ℓ)  D(W) ≥ X^2 vol^2.  Need to refute X ≥ vol^δ.

## Provable inputs about {S_ℓ}
(A) per ℓ: |S_ℓ| ≤ √N k^C (Ramanujan; e^{iℓα}1_W has ~1/ε harmonics, C–S over them).
(B) large sieve: Σ_{|ℓ|≤Λ}|S_ℓ|^2 ≤ (Λ + B)·B  (points spaced ≥ 1/B at best). Sharp for Λ ≥ B.
(C) whole sphere, isotropic kernel: Σ_W D(W) = (ε-pair count on S²) − main ≤ N^2ε^2·k^C/vol = N k^C.
(D) long-range exact: Σ_{ℓ mod Λ}|S_ℓ|^2 = Λ·#{pairs, equal azimuth mod 1/Λ} → Λ·B·2^{o(k)} as Λ→∞
    (equal azimuth ⟺ ūt / ū't' ∈ Q(√2)^+, divisor bound).
Combining (A) over |ℓ| ≤ Λ: Σ ≤ Λ N;  beats (B) only for Λ ≤ vol = B^{1/3}, where it gives exactly B^2.

## Positivity transfers: D(W_0) ≤ Σ_{W∈F} D(W) for any family F ∋ W_0 of bands with D ≥ 0.
Each transfer costs |F| (the main terms add) and gains whatever power saving the family statistic has.

T1 (window sum, same axis, 1/ε windows): Σ_j D(W_j) = ε-pair count on S² − main ≤ N k^C.
    ⇒ D(W_0) ≤ N k^C = B·vol · (vol·k^C).   Deficit factor vol^{1+δ}.  Max bound: X vol ≤ √N. [Ramanujan]
T2 (axis average, z in ε-ball): ε-close pair lies in band_z for z in a set of measure ε (a band of axes).
    D(W_{z0}) ≤ ε^{-2}∫_{|z−z0|≤ε} D(W_z) ≤ ε^{-1}·(N k^C) : deficit vol^{3}.  Worse than T1.
T3 (Hecke translates g ∈ Λ_s of the axis): Σ_g D(W_{gz}) = Σ_{x,y ε-close} #{g : gz ∈ band ⊥ x} − 2^s B vol.
    Orbit-in-band count = 2^s ε + O(2^{s/2}k^C) uniformly only for 2^s ≥ ε^{-2-δ} = vol^{4+δ}.
    Then Σ_g D ≤ 2^s N k^C/vol·ε... main-term error 2^s·B·k^C, needs 2^s ≤ vol^{1−δ}. Incompatible.
    Deficit ≥ vol^{4}/vol = vol^3.  Worse than T1.
T4 (anisotropic: band width ε, azimuth cell η ≥ ε): D_{ε,η}(W) ≈ η Σ_{|ℓ|≤1/η}|S_ℓ|^2, main ηB^2.
    Majorise by isotropic η-pairs on S² (loses η/ε) and sum over 1/ε bands: D_{ε,η}(W) ≤ N k^C.
    Poisson with power saving iff ηB^2 ≥ N^{1+δ}, i.e. η ≥ vol^{-1+δ}. Max bound (D)^{1/2} = √N again.
T5 (Gram matrix over windows): G(W,W') = Σ_ℓ ĉ_ℓ S_ℓ(W)conj(S_ℓ(W')) is PSD; Σ_{W,W'}G = Σ_ℓ ĉ_ℓ|S_ℓ(S²)|^2
    with S_ℓ(S²) the complete shifted convolution, ≤ √N k^C per ℓ ⇒ total ≤ N k^C.  Off-diagonal entries
    (same-azimuth pairs at different latitudes) have no sign, so no bound on the diagonal. Nothing.

## Alignment (what a rich cell forces on {S_ℓ})
Rich ε-cell N_c = X vol stays rich in the coarser (ε × η) cell for η ≤ εX/2: deviation ≥ X vol/2.
C–S over |ℓ| ≤ 1/η:  Σ_{|ℓ|≤1/η}|S_ℓ|^2 ≥ X^2 vol^2/(4η)  for all η ≤ εX/2.
At η = εX/2: RMS_{|ℓ|≤2vol^2/X}|S_ℓ| ≥ X vol/2.  Against (A): X ≤ 2 vol^{3/2}k^C ⇒ N_c ≤ vol^{5/2}: Ramanujan.
Against (B) at Λ = 2vol^2/X: (Λ+B)B ≥ X^2vol^2/(4η)·… gives X ≤ vol: determinant level. Nothing new.
Pure spike (all deviation in one cell) ⟺ S_ℓ ≈ X vol e^{iℓα_0} for all |ℓ| ≤ 1/ε: contradicts only (A).

## Verdict so far
Every positivity transfer with a provable (isotropic, Ramanujan-type) majorant returns the Ramanujan sup
bound √N = vol^{5/2}; the sign is already implicit in that bound. The sign cannot connect the short ℓ-range
Λ ∈ [vol, vol^2] (where the gap lies) to the long range Λ ≥ B where exact Parseval/divisor bounds hold.
