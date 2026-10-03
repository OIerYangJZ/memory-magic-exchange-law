# Brief for the cloud session: attack SA (updated 2026-10-03 after rounds 1–13; next session = rounds 14–17)

You are continuing `research/rate3`. Read, in this order: `HANDOFF.md`, `PROGRESS.md` (all rounds, especially 6, 8, 12 and 13),
`NOTES.md` §16, then `scratch/round6_model.md` and `scratch/round12_minimal.md`. Do NOT read everything else up front; pull in earlier sections only
when a specific route needs them. Work on the branch you are checked out on and push to it AND to main
(`git push origin HEAD HEAD:main`); main and prxq-prep are identical at the start of this session.

## Goal
Raise the unconditional exchange rate above 5/2 by proving (any piece of) Conjecture SA for r = 2
(NOTES §2, §7): Clifford+T states of T-count ≤ τ in an ε-cap number at most 2^{o(τ)}(1 + 2^τ ε²), uniformly in
the cap. Equivalent targets, any one suffices (NOTES §14 end, §15(b)):
5-point correlation ⇒ 13/5; in-band pair-count asymptotic with power saving ⇒ > 5/2; RMS_ℓ |S_ℓ| ≤ L^{2/3−δ} ⇒ > 5/2;
centred fourth moment ⇒ 3.

## Rules (from the user; binding)
1. **Direction is fixed: keep attacking SA.** Do not convert the work into write-ups, appendices or deliverables,
   and do not weaken the goal to restricted process classes or conditional results. If a route dies, pick the next
   route and keep attacking.
2. Run **4 rounds** in this session (rounds 14–17). One sub-goal per round. Each round ends with a new section appended to
   `PROGRESS.md` (sub-goal, what was done, intermediate conclusions, where it stops, next sub-goal), then
   `git add -A research/rate3 && git commit && git push origin HEAD`. Push after EVERY round so nothing is lost.
3. **Keep every thinking block short.** Never try to finish a long derivation in one thought. Write intermediate
   conclusions into `PROGRESS.md` (or a scratch file under `research/rate3/scratch/`) and continue from there.
   A session that hits the thinking limit produces nothing; that is the failure mode to avoid above all others.
4. Numerics are allowed when they decide a question cheaply (there is no server access from the cloud; compile
   with gcc, keep runs under ~20 minutes, k ≤ 12–13 for `numerics/sl.c`-style enumeration).
5. Do not touch `main.tex`.

## Where things stand (after rounds 1–13)
Unconditional rate still 5/2. Three orientation results, all rechecked locally (PROGRESS round 13):
(1) black-box insufficiency model (NOTES §16(c)): every input used so far, and eleven further true properties of the orbit
(§16(g)), are satisfied by a configuration with a vol² spike over one cap; (2) the K-determinant constrains only frames within
one fibre ε-arc, so the 5-point criterion is about cross-arc tuples on which no algebraic relation bites (§16(f)); (3) already the
ℓ = 0 statistic (band count; problem (M₀), §16(i)) is a power beyond every analytic method, and spectrally it is the statement
that the *local spectral measure of the Hecke operator at ẑ* has small Fourier coefficients at frequency t ≍ log N.
**(M₀) is a diagnostic only**: it neither implies nor is implied by SA and proving it would not move the rate. Do not make it
the goal. Withdrawn claims (rounds 4, 5) are listed in PROGRESS round 6; do not re-derive them.

## This session: ONE task — look at the spectral data nobody has computed, then use what it shows on SA
Exact spectral form of the SA quantity (NOTES_moments Thm S with x₂ = ẑ): for a cap of radius ε about z₀,
  N_cap(z₀) − main = Σ_{j≥1} ĉ_j(ε) Σ_{F ∈ B_j^C} λ_F(t) F(z₀) F̄(ẑ),
B_j^C = orthonormal Hecke eigenbasis of the Clifford-invariant harmonics of degree j, λ_F(t) = depth-t eigenvalue
(a polynomial in the degree-1 eigenvalue λ_F = 2√2 cos θ_F; Ramanujan |λ_F| ≤ 2√2). SA says sup_{z₀} of the left side is
≤ vol·2^{o}; Ramanujan gives √N. Everything is governed by the complex measures ν_j^{z₀} := Σ_F F(z₀)F̄(ẑ) δ_{θ_F}
(for z₀ = ẑ: the positive local spectral measure Σ_F |F(ẑ)|² δ_{θ_F}) and their Fourier coefficients at frequency t.

Round 14 (build + validate). Construct T₁ := 3·Π_C π_j(T) Π_C on V_j^C for all even j up to at least 400, preferably 1000–2000
(π_j via expm of the tridiagonal Lie-algebra matrices in the |j,m⟩ basis: T is diagonal e^{−imπ/4}, H = exp(−iπ(J_x+J_z)/√2);
Π_C = average over the 24 Clifford rotations; work in the m ≡ 0 mod 4 … subspace if that helps). Diagonalise. Validate hard:
(a) dim V_j^C against the octahedral Molien formula; (b) |λ_F| ≤ 2√2; (c) Z_j computed spectrally equals Z_j computed by direct
enumeration of states (`numerics/sl.c` or `orbit.c`) for t ≤ 20 and several j; (d) a cap count at a random z₀ both ways.
Do not proceed until (a)–(d) pass. Save code under `numerics/` and tables under `numerics/RESULTS.md` (Test F).

Rounds 15–16 (measure; each question gets a number, not an impression).
- Is {θ_F} for fixed j distributed by the Plancherel (Kesten–McKay, q = 2) law, and at what scale does it stop (level spacings)?
- Are the weights |F(ẑ)|² independent of θ_F? Porter–Thomas distributed? Any exceptional families (CM forms from Q(ζ₈), forms with
  F(ẑ) = 0 by symmetry, Clifford-point structure)? The same for F(z₀)F̄(ẑ) at generic z₀ and at the arithmetic points |+⟩, H±.
- Size of the Fourier coefficients ∫U_t(cos θ) dν_j^{z₀} as a function of (j, t), in the SA-relevant regime t ≈ (5/2)log₂ j and
  beyond it: square-root cancellation in dim V_j^C? Where does it fail? Sum over j ≤ 1/ε with the cap weights ĉ_j: how does the
  cancellation *across j* compare with the cancellation *within* each j? (SA needs both; Ramanujan uses neither.)
- Search for z₀ maximising |N_cap(z₀) − main| spectrally (gradient ascent on the explicit trigonometric expression) at t ≈ 20–30
  with ε critical; compare with vol and with 1/ε. This is a direct hunt for a rich cap beyond what enumeration reached.

Round 17 (use it). Whatever structure the data shows (exact identities, positivity, an exceptional family, a rigid relation
between θ_F and F(ẑ)), turn it into an inequality toward N_cap ≤ vol^{2−δ} or record precisely why it cannot be. If the data shows
no structure at all (generic random-matrix behaviour in every statistic), say so with the numbers; then spend the round on the one
place where arithmetic must enter: the *t*-dependence. λ_F(t) = 2^{t/2}U_t(cos θ_F)-type polynomials are the only thing that
distinguishes depth t from a random unitary evolution; look for any statement about Σ_F w_F U_t(cos θ_F) that is false for
random θ_F with the same density but true for Hecke eigenvalues (integrality of traces Tr T_{P^t} on V_j^C, Eichler–Selberg
class-number form of the trace, congruences between degrees j).
