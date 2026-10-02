# Round 12 scratch: the minimal open problem below SA (ℓ = 0), stated exactly

T-count units: N ≍ 2^t states, ε = N^{-2/5}, J := 1/ε = N^{2/5}, B = εN = N^{3/5}. Latitude z(x) = (|u|²−|t|²)/n = 2σ₁(m)/n − 1, m = |u|².
Zonal (Clifford-axis) Weyl sums Z_j := Σ_x Y_{j0}(z(x)) = Σ_{m ≺ n} r(m)r(n−m) Y_{j0}(2σ₁m/n − 1).

## Chain of reductions (each step rigorous)
(1) Band count for ALL bands with power saving  ⟸  variance over bands:  max_W|B_W − main| ≤ (Σ_W(B_W − main)²)^{1/2}, and the
    Poisson value of the variance is ≍ N·k² ≪ (εN)² = N^{6/5}, so any bound Σ_W(B_W−main)² ≤ N^{6/5−δ} gives the band asymptotic.
(2) Σ_W (B_W − main)² ≍ ε Σ_{1≤j≤J} |Z_j|²  (band-step functions resolve degrees ≤ 1/ε; ĝ_{W,j} ≍ ε, Σ_W ĝ_{W,j}ĝ̄_{W,j'} ≈ εδ_{jj'}).
    So the band asymptotic follows from  Σ_{j≤J}|Z_j|² ≤ N·J^{3/2−δ}  (= N^{8/5−δ'}).  Truth (random, Test E): ≍ N·J.  Ramanujan: ≤ N k² J².
    Gap to close: J^{1/2+δ} = N^{1/5+δ} on average over the zonal degrees j ≤ N^{2/5}.
(3) Equivalent 1-dimensional form: pair correlation of latitudes at scale ε: #{(x,y): |z(x)−z(y)| ≤ ε} = N²ε·(1 + O(N^{-2/5−δ})).
    Equivalent divisor form: Σ_{h ∈ Z[√2], |σ₁h| ≤ εn, |σ₂h| ≤ n} [C(h) − main(h)] ≤ N^{6/5−δ},  C(h) := Σ_m r(m)r(n−m)r(m+h)r(n−m−h):
    a 4-fold divisor correlation at shift h, averaged over ≍ N^{3/5} shifts short in σ₁. Needed per shift on average: error ≤ N^{3/5−δ}
    (each C(h) has N terms, Poisson error √N k²: room N^{1/10}).
(4) Spectral form. Z_j = Σ_{f∈H_j^C} ⟨Y_{j0},f⟩ f(ẑ) λ_f^{(t)},  λ_f^{(t)} = 2^{t/2}U_t(cos θ_f) (depth-t Hecke eigenvalue, Satake angle θ_f).
    ‖Π_C Y_{j0}‖² = (1/24)Σ_{g∈C} P_j(ẑ·gẑ) = (1+(−1)^j)/6 + (2/3)P_j(0): exactly 0 for odd j (z ↦ −z symmetry, X gate), ≈ 1/3 for even j.
    So Z_j (j even) is a weighted sum over the ≍ j/12 Clifford-invariant forms of degree j of U_t(cos θ_f), with explicit weights
    c_{j,f} = ⟨Y_{j0},f⟩f(ẑ), Σ_f|c_{j,f}|² ≍ j/(12·4π)·(1/3). Ramanujan is Cauchy–Schwarz with |U_t| ≤ t+1. The needed saving J^{1/2}
    on average over j is square-root cancellation among the U_t(cos θ_f), f ∈ H_j^C: equidistribution of the ≍ j/12 Satake angles
    at the scale 1/t ≍ 1/log N, weighted by c_{j,f}. This is the vertical Sato–Tate statement of §14(b) at its unavailable scale
    (Serre-type results have discrepancy ≍ 1/log, i.e. exactly this scale, with no power saving).

## Literature (WebSearch snippets only)
Full range over Q: Σ_{m≤x}d(m)d(m+h) = main + O(x^{2/3+ε}) (Motohashi, ASENS 1994, via Kuznetsov); conjectured x^{1/2+ε}.
Short intervals: only mean-square results (Jutila; Ivić: E[Δ₂(x,H)²] ~ H·F₃ for X^ε < H < X^{1/2−ε}, for the divisor problem itself;
Ivić–Motohashi: Σ_{f≤F}|D(N,f) − main|² ≪ N^{4/3+ε}F^{1/3}, F ≤ N^{1/2}, i.e. per-shift error ≈ N^{4/3}F^{−2/3} ≥ N·F^{−2/3}·N^{1/3}).
Matomäki–Radziwiłł–Tao: asymptotics for almost all shifts h ≤ H, H ≥ (log X)^{C}. No individual short-interval result below x^{2/3}.
Analogue of (3) over Q: per-shift error needed N^{3/5−δ} against Ivić–Motohashi's N^{4/3}F^{−2/3}|_{F=N^{3/5}} = N^{14/15}: gap N^{1/3}.

## Statement of the minimal open problem (ℓ = 0)
  (M₀)  For the single Hecke orbit X_n ⊂ S² (states of T-count t, N ≍ 2^t) and the Clifford axis ẑ:
        Σ_{j ≤ N^{2/5}} |Σ_{x∈X_n} Y_{j0}(ẑ·x)|²  ≤  N^{8/5−δ}  for some δ > 0.
  (M₀) ⟹ band-count asymptotic at the critical scale (all bands). Ramanujan gives N^{9/5}k². Numerically the left side is ≍ N^{7/5}.
  (M₀) involves no azimuth, no cancellation over ℓ, no cusp forms: only the latitude distribution of one orbit. It is strictly
  easier than every mean-square route to SA (which needs the same for the CM-twisted sums S_ℓ, 0 < |ℓ| ≤ 1/ε, i.e. for the
  non-zonal harmonics Y_{jℓ}), and nothing proves it.
