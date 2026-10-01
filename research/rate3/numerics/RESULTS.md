# Numerics for NOTES §10(A): the shifted convolutions S_ℓ (2026-10-01)

Code: `sl.c` (Tests A and B), `caps.c` (caps around fixed centres), `brute.py` (brute-force check), `analyze.py`.
Ran on school-server in `~/mmx_claude/rate3num` (outputs `c*` at β = 4/5, `d*` at β = 2/3).

`sl.c` agrees with brute force for k = 5, 6, 7:
- S_ℓ agrees to 1e-12.
- The band lengths L agree exactly.
- Cell counts agree except at the fold φ ≡ 0 (mod π/4), where floating-point rounding decides the bin.

The number of orbit pairs is exactly 2^{2k−1} − 1, so there are 8(2^{2k−1} − 1) states.

## Setting

- States are pairs (u, t) ∈ Z[ω]² with |u|² + |t|² = 2^k. Each sde corresponds to T-count ≈ 2k.
- Bands: W_j = {m : σ₁(m)/2^k ∈ [j/nb, (j+1)/nb)}, with nb = 2^{βk}. In z the band width is 2/nb.
- Frequencies: ℓ = 8n ≤ nb.
- The 5/2 critical scale is β = 4/5. The binding scale for α = 3 is β = 2/3.
- Normalised statistic: Z = S_ℓ/(8√S_0), where S_0 is the number of (u, t) pairs in the band. Under random
  phases, E Z² = 1.
- Bulk means |z| ≤ 0.8. Mirror bands j ↔ nb−1−j carry identical S_ℓ.

## Test A: S_ℓ at the critical scale β = 4/5

| k | nb | mean L | E Z² | E Z⁴/(E Z²)² | \|Z\| 99.9% | \|Z\| 99.99% | max \|Z\| (Gauss max) |
|---|----|--------|------|--------------|-------------|--------------|-----------------------|
| 12 | 776 | 913 | 0.99 | 3.34 | 3.44 | 4.20 | 5.54 (4.69) |
| 13 | 1351 | 1917 | 1.01 | 3.48 | 3.69 | 4.50 | 5.30 (4.92) |
| 14 | 2353 | 4056 | 1.00 | 3.42 | 3.66 | 4.48 | 5.52 (5.14) |
| 15 | 4096 | 8641 | 0.99 | 3.45 | 3.66 | 4.47 | 5.60 (5.35) |
| 16 | 7132 | 18503 | 1.00 | 3.47 | 3.69 | 4.54 | 6.07 (5.56) |

For reference, the Gaussian quantiles are 3.29 (99.9%) and 3.89 (99.99%).

- The normalised distribution does not depend on k. Its tails are slightly heavier than Gaussian, which matches the
  divisor-type weights r(m)r(m').
- From k = 13 to 16, L grows 9.7× while the 99.99% quantile moves from 4.50 to 4.54.
- If θ were 1/3, Z would have grown by about L^{1/6}/√k ≈ ×1.3–1.4 over this range.
- The raw maximum grows like L^{0.58}, because S_0/L ∝ k (the density of pairs of norms is ~1/k).
- **Effective θ = 1/2 up to logarithms.** In the dictionary of §10(A) this would mean rate 3. It is far from the
  θ = 1/3 line that decides 5/2.
- **At β = 2/3** (k = 12, 14, 16):
  - E Z² = 0.95, 1.00, 1.02 and the kurtosis is 3.3–3.4.
  - max |Z| = 3.65, 4.47, 6.08 against Gaussian maxima 4.19, 4.61, 5.00. The k = 16 maximum is a single band with
    heavy-tailed terms.

## Test B: direct counts in cells of side 2/nb, critical scale

| k | λ (volume) | var/mean | max cell (Poisson max) | max over 2×2 windows (Poisson max) |
|---|------------|----------|------------------------|------------------------------------|
| 12 | 35.5 | 1.006 | 65 (66) | 196 (198) |
| 13 | 46.8 | 1.014 | **112** (83) | 261 (255) |
| 14 | 61.7 | 1.002 | 104 (105) | 350 (328) |
| 15 | 81.5 | 1.000 | 131 (132) | 419 (423) |
| 16 | 107.5 | 1.002 | 166 (168) | 551 (545) |

- The determinant bound allows about (π/2)·nb ≈ 1.1·10⁴ states per cap at k = 16. A "5/2 is sharp" scenario would
  need caps of that size; the observed maximum is 166.
- At β = 2/3 the counts are also Poisson: var/mean 0.96–1.04, and the maxima match the Poisson maxima for k = 14, 16.

## Special points (caps.c): rich circles

The k = 13 outlier sits at |+⟩ and at the H-eigenstates (z = ±1/√2, φ ≡ 0).

Measured excesses:

| k | centre | radius | count | volume | ratio |
|---|--------|--------|-------|--------|-------|
| 11 | H± | ρ₀/2 | 37 | 5.3 | 7.0 |
| 12 | \|+⟩ | ρ₀ | 153 | 27.9 | 5.5 |
| 13 | \|+⟩, H± | 2ρ₀ | 345 | 147 | 2.35 |

Here ρ₀ = 2^{−0.8k}. For k = 14–16 the ratio is ≤ 1.4 at every radius. Random centres show ratios of
0.5–1.3 (Poisson noise).

Mechanism:
- All states with the same |u − t|² lie on one circle about the x-axis.
- Such a circle carries ≈ r(m_d) r(2^{k+1} − m_d)/64 states, a divisor-type multiplicity.
- A cap of radius ~ε at an arithmetic point swallows the whole circle.
- This is the S² version of the S³ "sphere-section clusters". It is bounded by 2^{O(k/log k)}, as the band theorem
  says.

Conclusion: the volume law can hold only with a 2^{o(k)} factor, not with an absolute constant. This does not affect
the rate.

## Verdict

- The data strongly support square-root cancellation in S_ℓ (θ = 1/2), and hence conjecture SA. There is no sign of
  any power-size excess at k ≤ 16, at either scale.
- The numerics prove nothing.
- What is provable is the mean square over ℓ (pair correlation within a band). It is exactly the L² statement that is
  blind at the critical scale.
- The input actually needed is a pointwise bound with θ > 1/3. That is of subconvexity strength for short shifted
  convolutions of CM Hilbert forms.
