# Brief for the cloud session: attack SA (updated 2026-10-02 after rounds 4–8; next session = rounds 9–12)

You are continuing `research/rate3`. Read, in this order: `HANDOFF.md`, `PROGRESS.md` (all rounds, especially 6 and 8),
`NOTES.md` §15–§16, then `scratch/round6_model.md`. Do NOT read everything else up front; pull in earlier sections only
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
2. Run **4 rounds** in this session (rounds 9–12). One sub-goal per round. Each round ends with a new section appended to
   `PROGRESS.md` (sub-goal, what was done, intermediate conclusions, where it stops, next sub-goal), then
   `git add -A research/rate3 && git commit && git push origin HEAD`. Push after EVERY round so nothing is lost.
3. **Keep every thinking block short.** Never try to finish a long derivation in one thought. Write intermediate
   conclusions into `PROGRESS.md` (or a scratch file under `research/rate3/scratch/`) and continue from there.
   A session that hits the thinking limit produces nothing; that is the failure mode to avoid above all others.
4. Numerics are allowed when they decide a question cheaply (there is no server access from the cloud; compile
   with gcc, keep runs under ~20 minutes, k ≤ 12–13 for `numerics/sl.c`-style enumeration).
5. Do not touch `main.tex`.

## Where things stand (after rounds 1–8)
Every route in NOTES §3–§16 stops at 5/2. The sharpest statement of why is the **black-box insufficiency model** (NOTES §16(c),
`scratch/round6_model.md`): a Poisson-like background plus one spike of 1/ε = vol² frames over a single critical cap, one per
ε-arc of the Hopf fibre, satisfies every input used so far (Ramanujan Weyl sums at all frequencies, the 5/2 tube bound, divisor
bounds on arithmetic circles, Liouville/L1 separation, the band volume law, the exact Hecke recursion of deviation fields, the
σ₂ shadow). So no combination of those inputs by inequalities can give N_c ≤ vol^{2−δ}; a proof must use a property that is
FALSE in that model. Also now known (round 8, NOTES §16(e)): the per-ℓ bound |S_ℓ| ≲ √N for the short window is unproved
(direct harmonics give k√N·ε^{−1/2}, worse than trivial at the critical scale); only the mean square over ℓ is at level √N.
Even the band count B_W = Σ_{m∈W} r(m)r(n−m) has no proved asymptotic at the critical scale, only the divisor upper bound.
Two claims of the previous cloud session were withdrawn in round 6 (prefix splitting does not beat the flat 5/2 tube bound;
anisotropic cells give no single-band asymptotic). Do not re-derive them.

## Suggested round-9 sub-goal (then follow where it points)
The model is killed exactly by the §12–§15 criteria and by nothing weaker, so the work must produce a new *input*, not a new
*rearrangement*. Candidates, in order:
(i) **5-point correlation ⇒ 13/5 via the fibre.** Five frames in a common ε-arc of a Hopf fibre have vanishing 5-point
    K-determinant (this is the engine of the 5/2 proof, Thm rate52 / App. E–H). Σ_c N_c⁵ counts 5-tuples of frames whose
    *states* share an ε-cap but whose fibre coordinates are arbitrary. Work out exactly what the determinant method says about a
    5-tuple of frames in a common cap but in *different* fibre arcs (same cap, five different arcs): is there any algebraic relation
    (a 5-point determinant in the Hopf-lifted coordinates, a Plücker relation, a Cayley–Menger identity) forced to vanish at the
    critical scale? If yes, count the 5-tuples it leaves; if no, record the exact scale where it starts to vanish.
(ii) **What is false in the model.** List concrete true properties of the Hecke orbit not in the black-box list and test each
    against the spike: e.g. exact multiplicativity r_O(P^τ l) = r_O(P^τ)r_O(l)/48 with the induced factorisation of a rich cap's
    frames; the exact identity Σ_W D(W) = whole-sphere pair deviation (sum over bands is a sum of squares with *known* main
    term); unique factorisation of the five frames' quotients. The goal is one property that the spike violates and that is
    provable for the orbit.
(iii) **Diagnostic, not a goal:** the ℓ = 0 case, an asymptotic for the band count B_W at the critical scale (shifted convolution
    of the weight-(1,1) theta series of Z[ζ₈] over Q(√2) at shift n = 2^k in a σ₁-window of relative width ε = N^{−2/5}). If no
    method reaches even this, say so precisely; it bounds from below the difficulty of every mean-square route.
