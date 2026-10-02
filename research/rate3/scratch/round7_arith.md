# Round 7 scratch: the rich cap in (u,t) coordinates and the exact arithmetic of Z[ζ₈]

Frames V = [[u, −t̄],[t, ū]], u,t ∈ Z[ζ₈], |u|²+|t|² = n = 2^k (both embeddings). State (u:t); fibre coordinate φ = arg σ₁(u) mod π/4;
band coordinate m = |u|² ∈ Z[√2] (window W: σ₁(m) ∈ M₀(1±ε), σ₂(m) ∈ [0,n]); azimuth α = arg σ₁(ūt) = arg t − arg u.
Band frames = ⋃_{m∈W} G_m, G_m = Θ(m) × Θ(n−m) in the (φ, arg t) torus, Θ(m) := {arg u : |u|² = m}, |Θ(m)| = r(m) ≤ 2^{o}.
Cap = strip {α ∈ α₀ ± ε}. N_c = Σ_{m∈W} c_m, c_m := #{(a,b) ∈ Θ(m)×Θ(n−m) : b − a ∈ α₀ ± ε} ≤ 2^{o}.
⇒ N_c ≤ 2^{o}·#{m ∈ W : (Θ(n−m) − Θ(m)) ∩ (α₀ ± ε) ≠ ∅}. A spike M = vol² needs a fraction 1/vol of W (expected fraction ε = 1/vol²).

Unique factorisation: m = Π𝔭_i^{e_i} in Z[√2]; split 𝔭_i = 𝔓_i𝔓̄_i in Z[ζ₈] with Hecke angle θ_i := arg(gen 𝔓_i) (mod π/4);
inert primes need even exponent and contribute angle 0; the ramified prime (√2) contributes a fixed angle.
Θ(m) = { Σ_i (2a_i − e_i) θ_i : 0 ≤ a_i ≤ e_i }.  So the azimuth of a frame is a signed sum of prime Hecke angles of m and n−m:
  α = Σ_{𝔮 | n−m} ±θ_𝔮 − Σ_{𝔭 | m} ±θ_𝔭.
SA (critical scale) ⟺ the multiset A_W := { Σ_{𝔮|n−m}±θ_𝔮 − Σ_{𝔭|m}±θ_𝔭 : m ∈ W, all signs } (size ≈ |W|·2^{o} = B) is
2^{o}-uniform at scale ε = B^{-2/3}: every ε-arc gets ≤ 2^{o}·εB = 2^{o} vol.  Fourier: Σ_{a∈A_W} e^{iℓa} = S_ℓ exactly. Full circle.

Determinant structure in these coordinates: frames of the cap have fibre angles arg u with ≤ 2^{o} per ε-arc. For the u's:
u ∈ Z[ζ₈] in R := {|σ₁u|² ∈ W} × {|σ₂u|² ≤ n} (radial thickness at σ₁: ε√M₀ ≈ 2^{-3k/10} < 1, i.e. R is a union of ≈ |W|
lattice circles), selected by arg u ∈ Θ(n − |u|²) − α₀ ± ε. The spike = correlation between arg u and the Hecke angles of the
complementary norm n − |u|², over u ∈ R. Random model: independent ⇒ M ≈ εB = vol.

Exact identities available: B_W = Σ_{m∈W} r(m)r(n−m) (band count), r(m) = multiplicative (divisor-type), Lagrange identity
|ūu'+t̄t'|² + |ut'−tu'|² = n² (quotients are frames of norm n² near T_z), Bloch vector (2Re ūt, 2Im ūt, |u|²−|t|²) on X²+Y²+Z² = n².
Each rewrites the cap count as a sub-sum of an exact divisor sum selected by an angle condition; none produces a second identity
for the angle-selected part. The angle-selected count is the shifted convolution of CM coefficients S_ℓ, i.e. §10(A).

Verdict: (H4) adds the Hecke-angle-sum description and nothing operational. The problem is equivalent to small-scale
(scale |W|^{-2/3}) equidistribution of signed prime-angle sums over the shifted pair (m, n−m), m in a short σ₁-window: a binary
additive problem with multiplicative (Hecke-angle) data at a single shift n = 2^k, exactly the object agent_burgess §4 calls
"a non-negative incidence kernel with all oscillation in completely multiplicative coefficients".
