> **Historical working notes (2026-09-28).** Theorem numbers below refer to the `main.tex` of
> that date. The summary's statement that the unconditional rate in the full CW model stays
> at 2 is superseded: the paper now proves every rate below 5/2 (Theorem 4), and results A,
> A1, B, D became Theorem 6, Corollary 2, Theorem 7 and Proposition S1. See `README.md` in
> this directory.

# Research notes: beyond the square-root barrier, and the ancilla model

Status: working notes, 2026-09-28. `main.tex` is **not** modified. All numbered
objects (Theorem 3, Conjecture H, Lemma 7, ...) refer to the current `main.tex`.
Scripts in this directory: `zomega.py` (exact Z[ω] arithmetic), `check_lde.py`,
`check_identity.py`, `check_axis.py`, `batch_qrom_cost.py`.

## 0. Summary

| # | Statement | Status | Frames / model |
|---|---|---|---|
| A | Volume law for the tube around arithmetic cosets of bounded height: #{W : t(W) ≤ τ, W ε-near G₁R_zG₂} ≤ 2^{O(τ/log τ)+O(h₀)}(τ+1+ε²2^τ) | proof below, numerically consistent | G₁,G₂ ∈ Γ with t(G₁)+t(G₂) ≤ h₀ (Clifford: h₀=0) |
| A1 | For all but a 2^{-Ω(L)} fraction of angles, every ε-approximation of R_z(θ) has T-count ≥ (3−o(1))L | corollary of A | ancilla-free synthesis |
| B | "Rotation-segmented" CW processes: α ≥ 3 − O(1/log L) unconditionally, and quantization with fee (2−o(1))L | corollary of A | each committed segment ε-close to a Clifford-framed z-rotation |
| C | One-sided cosets R_z·G (G arbitrary in SU(2)): count ≤ 2^{O(τ/log τ)}(τ+1+ε²2^{3τ/2}); hence an unconditional fee (4/3−o(1))L | proof below | r = 2 (and round 1 for any r) in the **full** CW model |
| — | Two-sided frames of large height (the general CW model) | **not achieved**; reduces to small-scale equidistribution of Hecke points, beyond current methods | general |
| D | With clean ancillas + phase-gradient catalyst: committed magic O(L/log L) per share, i.e. α = O(1/log L) → 0. No Ω(m log(1/ε)) law survives ancillas | explicit construction, known primitives | ancilla-assisted |

The square-root barrier is **broken**:
- fully, for tubes around arithmetic cosets of bounded height (Result A);
- below the onset, for one-sided cosets with arbitrary frame (Result C).

It is **not** broken for the two-sided, large-height frames that an adversarial CW process can create. So the unconditional exchange rate in the full CW model stays at 2. What changes:
- the conjectural rate 3 becomes a theorem in a natural submodel (Result B);
- an unconditional quantization of 4L/3 holds for two-round streams (Result C);
- the ancilla question in Sec. 10 of the paper gets a negative answer (Result D).

## 1. Setup and two reductions

Write a single-qubit Clifford+T unitary modulo phase as
W = 2^{-k/2} [[u, −t̄ω^l],[t, ūω^l]], u,t ∈ Z[ω], uū+tt̄ = 2^k, with k minimal (KMM).

**Fact 1 (T-count vs. denominator exponent).** For minimal T-count τ one has
k ∈ {⌈τ/2⌉, ⌈τ/2⌉+1} for τ ≥ 2 (for τ = 1 also k = 0, the words C·T·C').
Checked exactly for all 73 680 words of T-count ≤ 10 (`check_lde.py`). The enumeration also
reproduces 72·2^t − 48 distinct unitaries. We only use k ≤ ⌈τ/2⌉+1, so 4^k ≤ 16·2^τ.

**Fact 2 (tube ⇔ cap).** For g = [[a,−b̄],[b,ā]] ∈ SU(2), min_θ d_proj(g,R_z(θ)) ≤ ε
implies |b| ≤ ε. Since sin(δ/2) = |b| for the angle δ between Ad(g)ẑ and ẑ, it follows that
W is ε-near G₁R_zG₂ ⇒ Ad(W)·p₂ lies within angle 2 arcsin ε of p₁, where
p₂ = Ad(G₂⁻¹)ẑ and p₁ = Ad(G₁)ẑ.

For V ∈ Γ with numerator (u,t) and exponent k, the axis Ad(V)ẑ equals ν/2^k with
ν = (ūt+ut̄, (ūt−ut̄)/i, uū−tt̄) ∈ Z[√2]³ and ⟨ν,ν⟩ = 4^k in *both* real embeddings
of Q(√2). The map V ↦ Ad(V)ẑ has fibres V⟨T⟩ of size ≤ 8. Both facts were verified exactly on the
ball of T-count 10 (`check_axis.py`).

## 2. Result A: arithmetic cosets of bounded height

**Theorem A.** For every ε ≤ 1/8 and τ ≥ 1,
#{W : t_min(W) ≤ τ, ∃θ d_proj(W,R_z(θ)) ≤ ε} ≤ 2^{O(τ/log τ)} (τ + 1 + ε² 2^τ).
If G₁,G₂ ∈ Γ with t(G₁)+t(G₂) ≤ h₀, the same holds for the tube around G₁R_zG₂ with an extra
factor 2^{h₀}. For Clifford G₁,G₂ there is no extra factor.

*Proof.* Fix k ≤ ⌈τ/2⌉+1. The tube condition is |t|²_{σ₁} ≤ ε²2^k, and uū+tt̄ = 2^k with
both terms totally non-negative gives |t|²_{σ₂} ≤ 2^k.

(i) *Counting t.* Write t = a + b i with a, b ∈ ½Z[√2], using ω = (1+i)/√2. Then a and b lie in
the box {|x|_{σ₁} ≤ ε2^{k/2}, |x|_{σ₂} ≤ 2^{k/2}}. Multiplication by a power of the unit 1+√2
is an automorphism of Z[√2] that makes the box square up to a factor 1+√2. So the box holds
≤ C(1 + ε2^k) points (and only 0 if ε2^k < 1, since |N(x)| ≥ 1). Hence #t ≤ C(1+ε2^k)².

(ii) *Counting u given t.* If t = 0 the word is a power of T. Otherwise n := 2^k − tt̄ ≠ 0 is
totally positive, and N_{Q(√2)/Q}(n) ≤ 4^k. Z[ω] is a PID, and its units of relative norm 1
are the 8 powers of ω. So #{u : uū = n} ≤ 8·#{ideals of Z[ω] of relative norm (n)}
≤ 8·d_{Z[√2]}((n)) ≤ 2^{O(k/log k)} by the divisor bound.

(iii) *Summing.* The word is determined by (u, t, l), with l ∈ Z₈. Sum over k ≤ ⌈τ/2⌉+1:
Σ_k (1+ε2^k)² ≤ 2(τ+3) + (4/3)·16·ε² 2^τ.

For G₁,G₂ ∈ Γ apply this to V = G₁⁻¹WG₂⁻¹. This map is injective, and t(V) ≤ τ + h₀. ∎

*Numerics* (`check_identity.py`, all words to T-count 14, L = 3…6):
- The cumulative identity-coset counts track 72ε²2^τ, the Haar value.
- The fibre over one t-class holds at most 64·(2–4) words: 64 from the T-power symmetries, times a divisor factor of 2–4.
- The number of t-classes is far below (1+ε2^k)².

**Corollary A1 (typical z-rotations need 3L, unconditionally).** Distinct grid rotations
have disjoint ε-balls. So #{d ∈ Z_Q : τ_min(d) ≤ τ} ≤ |tube words of cost ≤ τ|, and the fraction of
grid angles (or the Haar measure of θ) with τ_min ≤ (3−δ)L is ≤ 2^{−δL + O(L/log L)}.

For comparison:
- Selinger (§9 of [Selinger15]) says only that this bound is "not a priori clear" for z-rotations.
- Ross–Selinger call 3 log₂(1/ε) "the information-theoretic lower bound" in the typical case. They support it with a heuristic ("it appears", "we expect", §8.3 of arXiv:1403.2975); their rigorous lower bounds are per-instance.
- The paper's eq. (13) derives the same statement from Conjecture H.

A1 makes it unconditional up to o(L). **Novelty check needed:** the diagonal case reduces to sums of two norms over Z[√2]. Something close may exist in the Parzanchevski–Sarnak or Sarnak–letter literature on golden gates, and should be searched before claiming it.

**Result B (rate three in a natural submodel).** Call a CW process
*η-rotation-segmented* if every segment W^{(t)}_j with t < r satisfies
min_θ d_proj(W, C R_z(θ) C') ≤ η for some Cliffords C, C'. Per-rotation synthesis pipelines,
and the fractional-passthrough family with η = ε/r, are of this form. The constraint is on W
itself, so the frame of Corollary 2 plays no role. Theorem A (576 Clifford frames) gives the
profile n(τ) ≤ 2^{O(τ/log τ)}(τ+1+η²2^τ).

The proof of Theorem 5 with this profile in place of (H) gives:
- T̄_t ≥ κ_η (D_t − m a_η) − O(m log(1+T̄_t/m)), with κ_η = 1 + (2log₂(1/η) − O(L/log L))/log₂K and a_η = O(L/log L);
- quantization: every bit beyond O(L/log L) per coordinate is carried by a word of cost ≥ 2log₂(1/η) − O(L/log L).

With η = ε this gives α ≥ 3 − O(1/log L) *unconditionally*. Passthrough is in the model, which gives
α ≤ 3 under the synthesis hypothesis and α ≤ 6 unconditionally (the covering exponent of [PS18]).
So in this model the gap between 2 and 3 in the paper closes. The whole 2-vs-3 gap in the
general CW model comes from the freedom of an adversarial process to emit non-rotation segments
whose frames have large height.

*Caveat.* The O(L/log L) is the worst-case divisor bound. At ε = 10⁻¹⁰, log₂ d(n) for N(n) ≤ 4^{34}
can reach ≈ 18 bits (Nicolas–Robin). So B is an asymptotic statement and gives nothing new at
chemistry scale.

Possible fix: replace the worst case by an average. The quantity to bound is Σ_{t∈box} r(2^k − tt̄). A Nair–Tenenbaum-type bound for multiplicative functions of polynomial values over (unit-balanced) boxes would reduce the slack to polylog(L) above onset. This is not done here.

## 3. Result C: one-sided frames (two-round streams, full CW model)

In round 1 there is no prefix. So the committed word W lies in the tube around R_z·G₂ with
G₂ = (suffix)^† *arbitrary*. By Fact 2 this is the same as W⁻¹ẑ lying in a cap of radius
≈ 2ε around an arbitrary point p. So the count is ≤ 8·#{arithmetic axes ν/2^k, k ≤ ⌈τ/2⌉+1, in the cap}.

**Theorem C.** For every G ∈ SU(2),
#{W : t_min(W) ≤ τ, W ε-near R_z·G} ≤ 2^{O(τ/log τ)}(τ + 1 + ε² 2^{3τ/2}).

*Proof (determinant method, the Q(√2) analogue of Bourgain–Rudnick's cap lemma).* Fix k and
R = 2^k. Take four axes ν₁,…,ν₄ ∈ Z[√2]³ of norm 4^k whose σ₁-images lie in a cap of angular radius ρ.
- In σ₁: the differences have transverse part ≤ 2Rρ and normal part ≤ Rρ²/2, so |det(ν₂−ν₁,ν₃−ν₁,ν₄−ν₁)|_{σ₁} ≤ 6R³ρ⁴.
- In σ₂: the differences have length ≤ 2R, so |det|_{σ₂} ≤ 8R³.
- The determinant lies in Z[√2]. If ρ < 48^{−1/4}R^{−3/2}, its norm is < 1, so det = 0.

Therefore all axes in such a sub-cap lie on one affine Q(√2)-plane, i.e. on a "circle". Completing the
square, the points on a circle are representations of a Z[√2]-integer of height R^{O(1)} by a
totally definite binary form. Their number is at most the number of elements of fixed relative norm
in the maximal order of a CM quadratic extension of Q(√2), which is 2^{O(k/log k)}.

Covering the ε-cap by ≤ C(1+ε²2^{3k}) sub-caps and summing over k ≤ ⌈τ/2⌉+1 proves the claim. ∎

*Consequence.* Combine C with the spectral bound (Theorem 2) as n(τ) := min(C-bound, Theorem-2 bound). Then a committed word of cost
τ ≤ (4/3)L − O(L/log L) carries only O(L/log L) bits about its share. The argument of Theorem 5
then gives, *unconditionally* and for r = 2 (for general r, round 1 only, or cumulatively over the prefix):
- the expected number of coordinates whose committed word costs ≥ (4/3 − o(1))L is ≥ (D_1 − O(mL/log L))/log₂K.

This is a genuine quantization theorem beyond the spectral method, which gives no fee at all. The fee is
4/3 rather than the conjectured 2. The *rate* stays 2, because Theorem C is weaker than the spectral bound above
the onset τ ≥ 2L.

## 4. Why the general case resists

For a round t with both a prefix and a suffix, the frame is two-sided, with G₁, G₂ ∈ Γ of T-counts
h₁ and h₂ that the process controls. Every route we tried loses a factor exponential in h₁+h₂:
1. Reducing to Theorem A via V = G₁⁻¹WG₂⁻¹ loses 2^{h₁+h₂}.
2. Fibring O over the rational subspace G₁Q(ω)G₂ and its orthogonal complement loses the covolume of O ∩ G₁Q(ω)G₂, which is ≍ 4^{h₂} for a primitive G₂. This was checked in both fibration orders.
3. The sphere picture (points Ad(W)p₂ of height τ+h₂ in a cap around p₁) loses 2^{3h₂/2}.
4. Moving the prefix into the word (P = U^{<t}W, a one-sided frame) charges the prefix cost to the word again.

The underlying problem is to count Hecke points of an arithmetic point of large height inside an ε-cap
below the square-root scale, uniformly in the point. This is exactly the content of Conjecture H. It lies
beyond both the spectral method (Remark 2) and the elementary methods above: the analogous rational
problem, lattice points on spheres in small caps, is open at this scale.

So the unconditional rate in the full CW model stays 2 for now.

## 5. Result D: the ancilla-assisted model

1. **Qubit-memory loophole.** If ancillas may survive the cut, a process can copy the share into an
   ancilla register with X gates (Clifford, free) and read it back in round r. That gives T_pre = 0 at
   S = 0. Any ancilla version of the law must therefore charge the qubits that survive the cut as memory.
2. **Clean ancillas still break Ω(mL).** Split the m coordinates into blocks of g. For a block, a
   unary-iteration QROM [Babbush et al. 2018] costs 4(2^g−1) T gates, independently of the word
   length, because the data are written by CNOTs. It loads φ(x_B) = Σ_{j∈B}θ_j x_j into a b-bit
   register. Adding that register into a catalytic phase-gradient register costs 4(b−1) T gates
   [Gidney 2018]. Then the QROM is uncomputed. The ancillas are clean after each round, and the
   catalyst is independent of the shares. The committed cost per share is (8(2^g−1)+4(b−1))/g, which
   is 4L/log₂L·(1+o(1)) for 2^g ≍ L. So α_anc = O(1/log L) → 0. Numbers (`batch_qrom_cost.py`):

   | L | best g | T/share | ancilla-free 3L | T per bit |
   |---|---|---|---|---|
   | 33.2 (ε=10⁻¹⁰) | 4 | 68.0 | 99.7 | 2.05 |
   | 66.4 | 4 | 102.0 | 199.3 | 1.54 |
   | 100 | 5 | 133.6 | 300.0 | 1.34 |
   | 1000 | 7 | 720.0 | 3000 | 0.72 |

   About 3b ≈ 120 qubits suffice at ε = 10⁻¹⁰. This answers the question the paper leaves open in
   Sec. 10, negatively: T_pre = O(m log(1/ε)/log log(1/ε)) at S = O(log(mQr)). It also improves the
   "phase gradient + adder ≈ 4 per bit" row of Table 4.

   Caveats: the construction uses standard primitives and is likely folklore. It is outside the CW
   criterion (block unitaries), so it must be judged against the total-error criterion.
3. **Lower bounds with ancillas.** Counting Clifford+T circuits on N qubits gives only
   T ≥ (mL − S − 2N² − O(N))/(2N + 1). The true exponent lies between that and O(mL/log L); for
   comparison, Gosset–Kothari–Wu reach √(2^n L) for arbitrary n-qubit diagonals. A lower bound of
   the form Ω(mL/log(mL)) is plausible but not proven here.

## 6. What this would change in the paper (suggestions only, not applied)

- **Add Theorem A and Corollary A1** to Sec. 5 or 7. This makes eq. (13) unconditional up to o(L). It also shows Conjecture H is a theorem on identity and Clifford cosets, which leaves only the high-height frames conjectural.
- **Add Result B** as "rate three for rotation-segmented processes". This locates the entire 2-vs-3 gap in adversarial frames, a much sharper story than the current one.
- **Add Theorem C** as an unconditional quantization (fee 4L/3) for two-round streams.
- **Rewrite the ancilla discussion.** The law does not survive clean ancillas, α_anc = O(1/log L); ancillas that survive the cut must be charged as memory. This changes the paper's message: the constant-rate law is specific to the ancilla-free (or bounded-ancilla) setting. It should be said up front.
- All three new bounds carry O(L/log L) divisor slack. None of them improves the numbers of Table 3 at ε = 10⁻¹⁰.

## 7. Attempt on the general case (two-sided frames of arbitrary height)

Outcome: **not solved.** Below: why the full statement (rate 3) is out of reach, and a frame-independent
route that may give a partial result.

### 7.1 The full statement is at least as hard as a known open problem
Take r = 2. The committed word W must satisfy W·S ≈ R_z(grid), where the suffix S = S(σ, Fut) is fixed by the context. By Fact 2 this means W⁻¹ẑ lies in a 2ε-cap around p = Ad(S)ẑ.

A single frame p can serve every context: take S = S₀R_z(φ). So a lower bound of rate 3 for r = 2 needs, uniformly in p ∈ S², the bound
#{ν ∈ Z[√2]³ : ⟨ν,ν⟩ = 4^k, ν/2^k ∈ cap(p, 2ε)} ≤ 2^{o(k)}(k + ε²4^k).

The Γ-orbit of ẑ at height k consists of all primitive points of that "sphere": 2^{2k} points versus 72·2^{2k}/8 words. So this is precisely the Q(√2)-analogue of Linnik's problem for lattice points on spheres in shrinking caps, uniformly in the cap, at cap sizes containing only ~R^{2/3} points (R = 2^k):
- the Ramanujan/Deligne input, and GRH in the rational analogue, reach only the square-root scale;
- known unconditional small-scale results (Duke-type, Bourgain–Rudnick–Sarnak) are either for almost all caps or for caps of size R^{−δ} with small δ;
- the codimension-3 sibling (covering exponent 1 for golden gates, Sarnak's optimal strong approximation) is known only conditionally, e.g. from a twisted Linnik conjecture on Kloosterman sums (Browning–Kumaraswamy–Steiner; Sardari).

Conversely, a cap that violates the bound gives an r = 2 process with a better exchange rate. So for r = 2, rate 3 is essentially *equivalent* to this uniform small-cap conjecture. That makes it a sharper and more classical replacement for Conjecture H.

### 7.2 A frame-independent route: the determinant method on S³ (first step done, second step open)
Work in the quaternion picture: W = w/√2^k, w ∈ O, nrd w = 2^k, R' := 2^{k/2} ≍ 2^{τ/4}. The tube around any great circle G₁R_zG₂ is a thin neighbourhood of a geodesic in S³. Arithmetic of the frame is not used.

**Level 1 (done).** Let V_d be the space of polynomials of degree ≤ d restricted to S³: dim ≈ d³/3, total degree of a monomial basis ≈ d⁴/4. Near the core circle, with coordinates (s, u, v) (along, transverse), the leading jets of V_d are (u,v)-monomials of degree j times trigonometric polynomials in s of degree ≤ d−j. Their jet profile is (j, i) with i ≤ 2(d−j), with multiplicity j+1, so Σj ≈ Σi ≈ d⁴/6.

For D points of X'_k in a box of transverse size ρ and length λ (angular), Bombieri–Pila expansion gives
|det|_{σ₁} ≤ C_d R'^{d⁴/4} ρ^{d⁴/6} λ^{d⁴/6} and |det|_{σ₂} ≤ C_d R'^{d⁴/4}.
The determinant is in Z[√2], so it vanishes once ρλ < c_d R'^{−3−o_d(1)}. **Hence all points of the tube inside such a box lie on one hypersurface section of degree d, for every frame.**

Choose ρ = ε and λ = R'^{−3}/ε. The tube is then covered by 1 + εR'³ = 1 + ε2^{3τ/4} boxes. The volume-law expectation per box is εR' = ε2^{τ/4}, which is ≤ 1 throughout τ ≤ 4L. So the method is automatically tuned to about one point per box.

**Level 2 (open).** It would suffice to prove the following *Uniform section lemma*: points of X'_k on a degree-d section Y inside such a box number ≤ C_d R'^{o(1)} + O(τ). The O(τ) allows the (HT)^j cluster of Proposition 1, which lies on the core circle. Consequences:
- Tube profile n(τ) ≤ 2^{o(τ)}(τ + ε2^{3τ/4}) for **all** frames.
- An unconditional fee (4/3−o(1))L in the full CW model.
- An unconditional exchange rate α ≥ 8/3 − o(1): carrying b bits needs τ ≥ (4/3)(b+L), which is 8L/3 at b = L. This would beat the spectral 2 without reaching the volume-law 3.

Where it stalls: the second determinant step, on Y ∩ box.
- Sections of degree ≥ 3 and curves of degree ≥ 3 look manageable. A 1-dimensional step needs arc length < R'^{−2/e'}, and conics carry only divisor-bound many points.
- The bad case is Y a 2-sphere (degree 2) containing the core circle. That forces the core plane into a rational hyperplane, which is exactly what arithmetic frames can arrange. Then Y ∩ tube is a band around a great circle of a rational 2-sphere, a ternary version of the same problem. The plane method needs λ < r^{1/2}R'^{−2}, which is incompatible with λ = R'^{−3}/ε at τ ≈ 2L.

Handling this band would need either a separate argument for bands on rational 2-spheres (an induction on dimension) or smaller boxes, which weakens the exponent 3/4. I could not close this. Whether the final exponent still beats rate 2 is **unknown**.

Summary of 7.2: the method is rigorous at level 1 and frame-independent, and it is heuristically consistent with the volume law. It is the only route I found that does not hit the GRH/Ramanujan barrier. It is a research project of uncertain outcome, not a result.

## 8. Result E: an elementary bound beyond the spectral barrier for ALL frames (α ≥ 2.22)

Status: new; proof below. **Needs independent checking before any use.** Scripts:
`rigid_plane.py` (optimiser), `certify_plane.py` (explicit certificate), `rigid_AB.py` (variant).

**Setting.** Normalise numerators so that X'_t = {w ∈ O : nrd(w) = n_t} is the set of numerators of
T-count-t words, up to sign. Here n_t is a fixed totally positive element with |n_t|_{σ₁} ≍ |n_t|_{σ₂} ≍ 2^{t/2}
(the unit can be normalised because totally positive units of Z[√2] are squares). Write R' ≍ 2^{t/4}
for the common radius; w ↦ W is 2-to-1. For SU(2) lifts, d_proj(W,V) = min|w/|w| ∓ v|, so the ε-tube of
the coset G₁R_zG₂ is exactly the Euclidean ε-neighbourhood, in σ₁, of the great circle C = G₁R_zG₂ ⊂ S³.
**Arithmetic of G₁, G₂ is never used.** Coordinates of O lie in D₀⁻¹Z[√2]⁴, and the σ₂-image of every
point lies on a sphere of radius ≍ R'.

Exponents are in units of log R'; tube radius ε = R'^{−a₀}, i.e. a₀ = 4L/τ.

**Level 1 (affine dependence on S³).** Let a box B be the set of points of the tube within angular distance
R'^{−a'} of C (a' ≥ a₀) whose projection to C lies in an arc of length R'^{−b'}. Differences of points of
B have along-component ≤ 2λR', transverse components ≤ 2ρR', and radial component ≤ λ²R'
(λ = R'^{−b'}, ρ = R'^{−a'}). So for five points, |det(d₁..d₄)|_{σ₁} ≪ λ³ρ²R'^4 and |det|_{σ₂} ≪ R'^4. The
determinant lies in D₀^{−4}Z[√2], so it vanishes once **2a' + 3b' > 8**. All points of B then lie in one
affine hyperplane H over Q(√2), i.e. on a round 2-sphere Y = S³ ∩ H of σ₁-radius r_Y = R'^{1−g}
(or on a circle).

**Geometry of Y ∩ B.** If 2b' ≥ a', the box lies within ≍ ρR' of a line, and so does H ∩ B inside H.
Hence Y ∩ B lies in ≤ 2 patches of Y, each contained in a geodesic rectangle of intrinsic sides
e_s = min(λR', r_Y, C√(r_Y ρR')) and e_u = min(ρR', r_Y). The square root is the near-tangent
(rim) case of a sphere meeting a thin cylinder.

**Level 2 (coplanarity on Y).** Cover each patch by ≪ (1+e_s/x_L)(1+e_u/x_S) geodesic rectangles of sides
x_L ≥ x_S. Take four points in one rectangle. The 3-vector d₁∧d₂∧d₃ has σ₁-norm
≪ x_L x_S min(x_S, x_L²/r_Y) (sagitta) and σ₂-norm ≪ R'^3. Its coordinates lie in D₀^{−3}Z[√2], so it
vanishes once x_L x_S min(x_S, x_L²/r_Y) R'^3 < c. The points of that rectangle then lie on one circle.

**Level 3 (circles).** Points of X' on a circle S³ ∩ P, with P an affine 2-plane over Q(√2) spanned by
points of height ≤ R', are representations of an element of height R'^{O(1)} by a totally positive binary
form over Z[√2]. They number at most the number of elements of fixed relative norm in O_E, E a CM
quadratic extension of Q(√2), i.e. ≤ R'^{o(1)} by the divisor bound.

**Count.** #tube words of T-count ≤ τ ≤ 2^{o(τ)} R'^{E(a₀)}, where
E(a₀) = min over (a',b') of [ 2(a'−a₀) + b' + max_g min_{x_L,x_S} log(#rectangles) ],
subject to the level-1 and level-2 inequalities. Numerically (`rigid_plane.py`, `certify_plane.py`):

| τ | bits carried ≤ | spectral (Thm 2) | optimal (a', b'), worst-g cost |
|---|---|---|---|
| 1.5L | 0.355L | 0.75L | (2.68, 0.90), 0.02 |
| 1.8L | 0.603L | 0.90L | (2.24, 1.18), 0.13 |
| 2.0L | 0.785L | 1.00L | (2.00, 1.34), 0.23 |
| 2.2L | 0.976L | 1.10L | (1.82, 1.46), 0.31 (certified) |
| 2.23L | ≥ 1.00L | — | — |

**Consequence (unconditional, full CW model, any r, any frames).** Replace (R) and the profile of
Corollary 2 by this profile in Lemma 12 (Gibbs weight λ = 1/α, cutoff τ_* = α log₂K). Since τ/bits(τ)
decreases up to τ* ≈ 2.225L, this gives
T̄_t ≥ α (1 − o(1)) (m log₂K − S) − o(mL), with **α = τ*/L ≈ 2.22.**
This beats the spectral rate 2 of Theorem 3, which has been the only unconditional rate.

The worst case at level 2 is a medium-small sphere section (r_Y ≈ R'^{−0.3…−0.4}) through the box. It
costs R'^{≈0.3} extra rectangles; without it the scheme would give ≈ 5/2.

**Caveats.**
1. As with A–C, the o(τ) is a divisor-bound slack. The result is asymptotic and does not move Table 3 at ε = 10⁻¹⁰.
2. The exponent optimisation is numerical on a 0.01 grid; the value 2.22 is approximate (±0.01).
3. The patch-geometry step (covering Y ∩ B by the stated rectangles, the rim case, and 2b' ≥ a' at the optimum) is the part most in need of a careful written proof.
4. Adding high-degree determinants on the (rigid) spheres (`rigid_AB.py`) improves this only to ≈ 2.24, and requires a level-3 lemma for arbitrary curves on a sphere, which is not proved.
5. 8/3 is not reached. High-degree auxiliary surfaces at level 1 give the better box condition a'+b' > 3. But in anisotropic boxes an auxiliary surface can wiggle at amplitude ε over length λ, which destroys the along-core gain at level 2. I found no way around this.

### 8.1 Errata after independent review (verdict: VALID WITH FIXABLE GAPS; scripts in research/verify/)
- **Level-2 volume bound.** The correct bound is Θ(x_L³x_S/r_Y), not x_L x_S min(x_S, x_L²/r_Y). A counterexample is in `verify/sagitta_check.py`. With the corrected condition 3l + sh + g + 2 < 0 the exponents are unchanged, because at the optimum x_L²/r_Y < x_S.
- **Level-1 constraint.** The radial extent is R'(λ²/8 + ρ²/2), so the scheme needs **a' ≥ b'**. It also needs **2b' ≥ a'**, or else e_s = λR' at g ≈ 0: great spheres containing C give bands of length λR'. The rows at 1.2L and 1.5L in §8 violate this; the corrected values are 0.267L and 0.417L. τ* is unaffected.
- **Worst sphere.** It is g ≈ 1.09–1.16 (r_Y ≈ R'^{−0.1…−0.16}), not R'^{−0.3…−0.4}.
- **Corrected optimum.** On a finer grid, τ* ≈ 2.233L, i.e. **α ≈ 2.23**, and τ/f(τ) decreases on all of [0, τ*].
- **Remaining gaps.**
  - A full proof of the patch geometry. There is a clean argument via the sinusoid d(s) = |A cos(s − s₀) − c|.
  - An exact piecewise-linear certificate in place of the grid search.
  - Restating Lemma 12 for general λ.
  - Writing up level 3.
