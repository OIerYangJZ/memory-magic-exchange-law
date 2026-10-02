# Brief for the cloud session (2026-10-02): rounds 4–7, attack SA

You are continuing `research/rate3`. Read, in this order: `HANDOFF.md`, `PROGRESS.md`, `NOTES.md` §12–§15.
Do NOT read everything else up front; pull in earlier sections only when a specific route needs them.

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
2. Run **4 rounds** in this session. One sub-goal per round. Each round ends with a new section appended to
   `PROGRESS.md` (sub-goal, what was done, intermediate conclusions, where it stops, next sub-goal), then
   `git add -A research/rate3 && git commit && git push origin HEAD`. Push after EVERY round so nothing is lost.
3. **Keep every thinking block short.** Never try to finish a long derivation in one thought. Write intermediate
   conclusions into `PROGRESS.md` (or a scratch file under `research/rate3/scratch/`) and continue from there.
   A session that hits the thinking limit produces nothing; that is the failure mode to avoid above all others.
4. Numerics are allowed when they decide a question cheaply (there is no server access from the cloud; compile
   with gcc, keep runs under ~20 minutes, k ≤ 12–13 for `numerics/sl.c`-style enumeration).
5. Do not touch `main.tex`.

## Where things stand (one paragraph)
Every route in NOTES §3–§15 stops at 5/2, for three coinciding reasons: the Burgess endpoint (dual length q^{1/4}),
the norm-window barrier (Archimedean methods see ε^{−1/2} neighbouring norms and the union has one point per cell),
and the L²/ball-non-concentration blindness. Round 3 verified numerically that all moment criteria up to order 8
hold with constant 1 at k ≤ 16. The last unexplored structural fact is NOTES §15(d): the in-band pair-count
deviation is a non-negative sum of squares Σ_{ℓ≠0} ĉ_ℓ |S_ℓ|². Start there (round 4), then move to whatever the
round-4 conclusion points at.

## Suggested round-4 sub-goal
Exploit the sign in §15(d). Candidates: (i) a lower bound for Σ_ℓ |S_ℓ|² from a single rich cell, combined with a
provable upper bound on a *different* functional of the same family {S_ℓ} (e.g. a smoothed mean square over a
longer ℓ-range, which Ramanujan/large sieve may control), to force a contradiction for cells above volume^{1+δ};
(ii) the same with the Hecke-ball (depth-τ arithmetic centre) structure of §15(c) replacing the band.
