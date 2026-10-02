# Round 11 scratch: the band count B_W = Σ_{m∈W} r(m)r(n−m) at the critical scale (brief item (iii), diagnostic)

Setting (T-count t): n = |u|²+|t|² ≍ 2^{t/2}, N = #states ≍ n² = 2^t. Window W: σ₁(m) ∈ [M₀(1−ε), M₀(1+ε)], M₀ ≍ n, σ₂(m) ∈ [0, n];
|W| := #{m ∈ W} ≍ εn² = εN = 2^{3t/5}. r(m) = #{u ∈ Z[ζ₈] : uū = m} = 8·Σ_{𝔡|m} χ(𝔡) for (m, √2) = 1 (class number 1, 8 units;
χ = quadratic character of Q(ζ₈)/Q(√2)). Expected: B_W = c·|W|·(1 + o(1)) with c = mean of r(m)r(n−m) over the window (≍ 1, log-free
since Σ_{𝔡} χ(𝔡)/N𝔡 converges). Target: an asymptotic with ANY power saving at ε = N^{-2/5}.

## Route 1: hyperbola method + Hardy–Littlewood lattice discrepancy
r(m) = 8(1+χ(m))Σ_{𝔡|m, N𝔡<√N(m)} χ(𝔡) + O(1) (pairing 𝔡 ↔ m/𝔡). So
  B_W = 64 Σ_{N𝔡 ≤ n, N𝔢 ≤ n} χ(𝔡)χ(𝔢) #{m ∈ W : 𝔡|m, 𝔢|n−m, χ-conditions} + (boundary terms).
Each count is the number of points of a coset of the lattice 𝔡𝔢 ⊂ Z[√2] in the box W (sides εn at σ₁, n at σ₂). In the Z²-coordinates
(a,b) ↦ a + b√2 the box is a parallelogram with sides along the lines of slope ±1/√2: badly approximable slopes, so the lattice-point
discrepancy is O(log(area)) uniformly in the shift (Hardy–Littlewood 1922 / Ostrowski; the coset count reduces to Kronecker sequences
{kα + β} with α quadratic). Hence each count = |W|/(c N(𝔡𝔢)) + O(log n). Number of pairs (𝔡,𝔢): ≍ n·n = n² = N.
  Error ≤ N log N   against main |W| = N^{3/5}.   FAILS by N^{2/5}.
Random-sign heuristic for Σ χ(𝔡)χ(𝔢)·O(log): ≈ N^{1/2} log ≪ N^{3/5}: the asymptotic is "true heuristically" but needs cancellation
among N discrepancy terms, which no method provides (this is the Type-I/II obstruction of agent_burgess §4 in its simplest form).

## Route 2: harmonics + Ramanujan (§16(e), recorded here for comparison)
1_W(z) has ≍ 1/ε zonal harmonics, each Weyl sum ≤ k√N·√(2j+1)/√(96π)·(…); Cauchy–Schwarz gives error k√N ε^{-1/2} = N^{7/10}k.
  FAILS by N^{1/10}.  (The whole-sphere pair count succeeds because it is a sum of squares; the band count is not.)

## Route 3: shifted convolution of the weight-(1,1) Eisenstein series (additive divisor problem over Q(√2))
Σ_m r(m)r(n−m)W(m) = shifted convolution of E = θ_{Z[ζ₈]} (Hilbert Eisenstein series of parallel weight 1, character χ) with itself at the
single shift n = 2^{t/2}, with a window short in σ₁ (relative width ε) and full in σ₂.
(3a) Smooth full-range analogue (Kuznetsov/Motohashi over Q(√2)): error ≍ N^{1/2+θ+o} times the Sobolev cost of the window. The window
     has σ₁-frequencies ≤ 1/ε, so the spectral sum runs over Maass forms with t₁ ≲ 1/ε, t₂ ≲ 1: ≍ ε^{-2} = N^{4/5} forms; cost ≥ ε^{-1}.
     Even the optimistic Blomer–Harcos-type shape X^{1/2+θ}(X/Y)^{1/2} gives N^{1/2+θ}·(1/ε)^{1/2} = N^{7/10+θ}.   FAILS by ≥ N^{1/10}.
(3b) Sharp-cutoff full-range analogue (Estermann/Deshouillers–Iwaniec/Motohashi): error x^{2/3+o} for Σ_{m≤x}d(m)d(m+h); here x ↔ N:
     error N^{2/3} against the short window's main term N^{3/5}.   FAILS by N^{1/15}.
     The short window makes it harder, not easier (the interval length y = |W| = x^{3/5} is far below x^{2/3}, the classical threshold).
(3c) δ-method bookkeeping: moduli 𝔮 with N𝔮 ≲ √(X₁X₂) = n = N^{1/2}; Voronoi dual length N𝔮²/(εn²) ≤ 1/ε = N^{2/5} ≫ 1: long dual sums;
     Weil saves (N𝔮)^{1/2} per Kloosterman sum. Trivial accounting reproduces (3b).

## What would suffice
An exponent 3/5 − δ in the binary additive divisor problem over Q(√2) with a σ₁-short window — i.e. beating the Kloosterman/Kuznetsov
exponent 2/3 by 1/15 AND absorbing the window. Over Q the exponent 2/3 (sharp cutoff) has stood since Motohashi; 1/2 is conjectural.
Conclusion: no known method gives the ℓ = 0 band-count asymptotic at the critical scale. Since the ℓ ≠ 0 problems (CM × CM shifted
convolutions, cuspidal, no divisor structure) are strictly harder, every mean-square route inherits at least this gap.

## Numerical check (this round): fluctuations of B_W across bands at k = 12 (sl.c, S₀ = B_W/64 per band)
(filled in below from the run)
Run (cloud container, 4 cores): ./sl k 0.8 → S₀(W_j) = B_W/64 (orbit pairs) per band; bulk |z| ≤ 0.8; local std from neighbour differences.
| k | bands | mean S₀ | mean L | S₀/L | relstd (local) | relstd·√L | max dev / std |
| 12 | 620 | 6.92e5 | 913 | 758 | 0.0247 | 0.745 | 2.97 |
| 13 | 1081 | 1.59e6 | 1917 | 829 | 0.0190 | 0.832 | 2.34 |
| 14 | 1883 | 3.65e6 | 4056 | 900 | 0.0119 | 0.756 | 3.15 |
| 15 | 3276 | 8.39e6 | 8641 | 971 | 0.0082 | 0.761 | 3.16 |
relstd·√L is k-independent (≈ 0.76): B_W = main·(1 + O(L^{-1/2})), i.e. square-root cancellation in the number of norms, Gaussian tails.
Deviation ≍ |W|^{1/2} = N^{3/10}; best provable error ≥ N^{2/3} (Kuznetsov sharp) / N^{7/10} (harmonics). Gap ≥ N^{11/30}.
Proved range of the band asymptotic: harmonics give error k√N ε^{-1/2} < εN iff ε > N^{-1/3}k^{2/3}; critical scale N^{-2/5} is below it.
