# Referee report: integration of Theorem rate52 (α ≥ 5/2 − o(1)) into main.tex (2026-09-30)

**Scope.** I read `git diff main.tex` against `/tmp/main_before_rate52.tex`, and checked it against the reviewed sources: `research/annulus/{appG_draft.tex,NOTES.md,review/REVIEW.md}`, `research/uniform52/{appH_draft.tex,NOTES.md,review/REVIEW.md}` and `research/conditional3/{NOTES.md,review/REVIEW.md}`. I also checked App. E/F of main.tex, README.md and `qip/`.

**Reproduced.**
- **Compilation.** pdflatex ×3 in `/tmp/revbuild` gives 39 pp. There are no undefined references and no overfull boxes. One new warning appears: "A float is stuck". It is benign: all floats land on pp. 13–17.
- **Certificate.** `cert_proj.py 400/249 … 12` gives output byte-identical to `cert_249_100_output.txt`.
- **SMT check.** z3 5.1.0 in `/tmp`, `verify_z3.py`:
  - `1/10 10/33` and `1/10 1/2` are unsat;
  - `1/10 3/10`, `1/10 0.30302` and `1/1000 0.30302` are sat;
  - `1/10 0.30304` is unsat.
  So the threshold is at 10/33.
- **Ball example.** A portable copy of `ballpts.cpp`, run at radius exactly 0.0025 around g/|g|, finds **exactly 126 words**. All have T-count 25, and all lie at distance 0.00224823 on one section. The Haar mean at 0.0025 is 16.02.

## Verdicts

| item | verdict |
|---|---|
| 1a. Thm rate52 vs sources: conversion of the tube bound, A = O_α(1) | **VALID WITH FIXES** (range, F7) |
| 1b. Sharpness sentence (direction of c, "near the limit") | **VALID** (wording nit, F8) |
| 1c. Prop dense: statement and proof | **VALID** (minor precision, F11). Its paraphrases elsewhere are **OVERSTATED/imprecise** (F2) |
| 1d. Remark what52 and the new Discussion sentences: hedging | **VALID WITH FIXES** (F5, F6, F10). "Exact limit of the method" is **OVERSTATED** (F4) |
| 2. App G/H: notation, definitions, cross-refs, numbers, z3 citation | **VALID WITH FIXES** (F9, F12) |
| 3. Global consistency (stale 17/7, numbering, README, abstracts) | **VALID WITH FIXES** (F1, F3, F13) |
| 4. Other overstatements and typos | see F2, F4, F6, F8 |

### 1a. Theorem rate52: the conversion is correct
- **Exponent.** T = E1 + (10/33)μ = 8/5 − 4μ/11. With R' ≤ 2·2^{t/4}, the bound R'^{T+o(1)} is 2^{(2/5 − μ/11)t + o(t)}.
- **Range.** t ≤ 4L/a0 gives ε = 2^{−L} ≤ 2^{a0}R'^{−a0}, which is the appendix hypothesis ε ≤ CR'^{−a0}.
- **Gap.** The gap to 1/α = a0/4 is 15μ/44, which is fixed. Hence Z_λ = O_μ(1) and A = O_α(1), in agreement with the uniform52 review.
- **For α < 40/17 (μ > 1/10).** The "monotonicity in ε" remark covers it.
- **One gap in the range.** The Gibbs step needs the profile up to τ_* = ⌈α log₂K⌉ ≤ αL + 2.7. The main-text range t ≤ 4L/(8/5+μ) = αL stops short of that by up to 3 T-counts.
  - The appendix hypothesis ε ≤ CR'^{−a0} allows any constant C, so t ≤ 4(L+2)/a0 is covered.
  - This matches Theorem dioph and the finite instance, which both state t ≤ …(L+2). See F7.

### 1b. Sharpness
- Put 2/5 − cμ = (E1 + c′μ)/4. Then c = 1/6 − c′/4.
- So "sat (fails) for every c > 1/11" is the same as c′ < 10/33. The direction is correct.
- This is supported by the uniform52 review's exact z3 Optimize (sup = E1 + 10/33·μ at five values of μ) and by my reruns.
- "Near the limit" undersells the result: the exponent is sharp on all of (0, 1/10]. That is harmless.
- "So 5/2 is approached and not attained" can be misread as a statement about the true rate (F8).

### 1c. Proposition dense: the proof is complete
- **Product order.** Correctness uses the product W_u·V. That matches the paper's convention: in Cor. profile, G1 = (U^{<t})†, so earlier gates stand on the left.
- **Correctness.** W_uV ≈ R_z(Δu)GG⁻¹R_z(Δ(ρ+a2)) = R_z(Δx), up to the phase R_z(ΔQ) = −I. The error is ≤ ε by bi-invariance and the triangle inequality.
- **Snapshot.** It is ⌈m log₂γ⌉ bits. With a single cut, no round index or program counter is needed, so S ≤ m(log₂Q − (1−o(1))L) + 1. Since log₂Q − log₂K = O(1/Q), this is the stated bound.
- **Cost.** The process is deterministic and zero-error, and T_pre = T_1.

### 1d. Hedging against the review verdicts
- **Ball heuristic.** It is labelled ("Heuristically") and states the ball-covering conclusion only. That is the conditional3 review's "safe" wording. The exponents check against conditional3 NOTES §3.
- **The two issues in that sentence (F5).**
  - Its **H is N(nrd g)**, the NOTES convention. It is not the paper's height h(g) = |σ1g||σ2g| = N(nrd g)^{1/2}. In the example, N = 7 but h(g) = √7.
  - "**σ1-radius** below 2^{−t/α}" is a normalized radius. In App. E units the σ1-radius carries a factor R'.
- **Exact instance.** Correct, and re-verified.
- **"Compatible with Conjecture H".** It is labelled heuristic, as the review required (it is checked only for ε ≥ r).
- **Quadric paragraph.** The uniform52 NOTES label it "Heuristic (not a theorem)", and its review says it was **not reviewed**. The paper states it as fact ("recovers exactly"), hedged only by "do not obviously help" (F6).
- **Twisted Linnik sentence.** "As far as we can see … in their natural circle-method and spectral implementations" meets the review's requirement ("in this framework", no "cannot"). But it is attributed to Remark what52, which contains nothing about Linnik (F10).
- **"Exact limit of the method" (abstract, intro, What is new, Discussion, all three QIP/arXiv abstracts, extended abstract).**
  - Rigorous content: the method with one auxiliary hyperplane per box (κ = 0) cannot pass 5/2, because the box count R'^{b}, with b > (8−2a0)/3, must stay below R'^{a0}.
  - Heuristic content: κ > 0 fails only in a numerical model. The degree-2 variant (level-1 ceiling 28/11) was not attempted. The quadric obstruction is unreviewed.
  - The uniform52 review says explicitly: "the exponent of this assembly …, not of the determinant method at large."
  - The abstract also attributes "5/2 is where the method stops" to the SMT check. The SMT check proves the rates below 5/2 and the sharp exponent. The stop comes from the E1 box count (F4).

### 2. Appendices G and H
**Numbers.** They match `cert_249_100_output.txt`:
- a0 = 400/249 = 1.60643;
- b = 1192747/747000 = E1 + 1/1000;
- E_* = 319751/199200 = 1.60518;
- worst exponent 532897/332000 = 1.60511;
- E_*/4 = 0.40129 < 100/249.

**z3 citation.** `verify_z3.py 1/10 10/33` → unsat matches `z3_outputs.txt`. The Z3 bibitem (TACAS 2008, LNCS 4963, 337–340) is correct.

**Consistency with App. E/F.** These are consistent: γ, ψ, H_V (F5), ℓ_ε (the draft's "lde" is correctly renamed), X_L < 1 − g, y_j, G_3 and G_3′ (these match the review), and route N/P for d = 2 and d = 3 (the exponents re-derived from G3(a)/(b) agree).

**Route-P prefix.** The paper uses `b` where the draft had `E1`. This is correct: the certificate's "E1" is b, and the output prints E1 = 1192747/747000.

**Cross-references.** They resolve to the right lemmas: G1 = Lemma 24, G2 = 25, G3 = 26; App. G and App. H.

**Remaining gaps (F9, F12).**
- `pm` is not defined in App. F. The paper needs its full definition: max over tangent and crossing classes, offsets at y_j, points per section at y_{j−1}.
- W (d = 3) is not defined as the K-span of v1..v3.
- "A basis of V⊥∩O*" should be "the K-independent 𝔫1, 𝔫2 of F5(b)".
- \bar S should say "Lemma E2 at g_i".
- "Largest exponent" should add "the certificate checks (first admissible threshold)".
- App H defines E_1 only after its first use, writes b = E_1 in one line and b = E_1 + ϖ in the next, and uses T (the paper's E_*) and η (for the height y).
- App H omits two of the uniform52 review's fixes: the last-layer closing (fix 3) and the margin rationale (fix 5).
- δ := εR' in App G and Remark what52 collides with the error probability δ of eq. (rate52).

### 3. Global consistency
- **No stale "17/7 is best" and no "heuristic ceiling 5/2" claim remains**, with one exception at l. 393 (F3).
- **README numbering matches the compiled aux file.** Checked: Thm 8/9/10/15, Prop 4, Cor 3/6/7, Lemma 13, Eq. (33), App. F/G/H, Sec. 9/10.1.
- **Fig. 2(b)'s PNG legend is now wrong.** It says "Thm. 9 (Conj. H, slope κ = 2.82)". After the insertion, main2 is Theorem 10 and Theorem 9 is thm:frontier (F1).
- **Abstracts.**
  - The ASCII, Unicode and arXiv abstracts say the same thing. The arXiv one is 1917 characters.
  - The extended abstract's Theorem 1 matches Thm rate52 (the h₂ term is absorbed in O_α(m)).
  - The shared issues are F2 and F4.
  - `qip/short_abstract_registered_*.txt` and `short_abstract_option_b_unicode.txt` (17:01) still say 17/7. They are fine if they are records of the registered text; otherwise they are stale (F13).
- **Release notes.** These are not errors in main.tex:
  - `research/{annulus,uniform52,conditional3}/` are untracked, but Data availability and the README point to them.
  - `arxiv/v1/src/main.tex` is the old version (not touched).

## Fixes (old → new)

**F1 (must).** Fig. 2(b) legend. In `scripts/make_fig2.py` l. 64, change `'Thm. 9 (Conj. H, …)'` → `'Thm. 10 (Conj. H, …)'` and regenerate `fig2_frontier_twopanel.png`. Also update the docstring "Thm. 8" → "Thm. 9" (the frontier). Several script comments and outputs (`thm_numbers.py`, `mixing_bounds.py`) use old theorem numbers from earlier versions; this is cosmetic.

**F2 (must).** Every paraphrase of Prop dense omits the T-count. Without it the statement is vacuous: every tube carries words at every ε-step at T-count ≈ 3L. "Almost every" also misstates the maximal-gap hypothesis.
- **Abstract:** "a tube carrying words at almost every $\varepsilon$-step of its core would make $\tfrac52$ sharp" → "a tube carrying words of $T$-count $(\tfrac52+o(1))L$ with gaps of at most $2^{o(L)}$ grid steps along its core would make $\tfrac52$ sharp".
- **Same change in:**
  - the intro item (l. 74: "a tube with words at almost every $\varepsilon$-step of its core");
  - the Discussion (l. 737: "tubes that carry words at almost every $\varepsilon$-step of their core");
  - `qip/extended_abstract.tex` Theorem 1.

**F3 (must).** l. 393: "the source of the gap between Theorems~\ref{thm:rate177} and~\ref{thm:main2}" → "…between Theorems~\ref{thm:rate52} and~\ref{thm:main2}".

**F4 (should).** Qualify "the limit of the method", and fix the SMT attribution.
- **Abstract:** "raises it to every $\alpha<\tfrac52$. A continuum form of the resulting case analysis, checked exactly by an SMT solver, shows that $\tfrac52$ is precisely where the method stops, and …" → "raises it to every $\alpha<\tfrac52$, as an SMT solver checks exactly on a continuum form of the resulting case analysis. This is where the method stops: at $\tfrac52$ its boxes shrink to $\varepsilon$-cubes, and …".
- **Intro item:** "This is the limit of the method:" → "This is the limit of the method with one auxiliary hyperplane per box:".
- **What is new:** "which is its exact limit" → "which is the exact limit of this case analysis".
- **Discussion:** "It is the exact limit of that method, because" → "It is the exact limit of that method as set up here, with one auxiliary hyperplane per box, because".
- **QIP/arXiv abstracts:** ", and an exact SMT check of its continuum case analysis shows that the method stops at $5/2$." → ", by an exact SMT check of its continuum case analysis; $5/2$ is where this method stops." The arXiv abstract becomes 1913 characters.
- **Extended abstract:** "the exact limit of the method" → "the exact limit of this method (one hyperplane per box)".
- **Optional:** "We determine the exchange rate" (all abstracts) overstates what an interval [5/2, 3] shows → "We bound the exchange rate".

**F5 (should).** In Remark what52, change "near rational points $g/|g|$ of height $H\asymp2^{t(4/\alpha-1)/3}$ there are sphere sections $\{\dots\}$ of $\sigma_1$-radius below $2^{-t/\alpha}$ that carry about $2^{t(2/3-5/(3\alpha))}$ words each" → "near rational points $g/|g|$, $g\in\mathcal O$ with $N(\mathrm{nrd}\,g)\asymp2^{t(4/\alpha-1)/3}$, there are sphere sections $\{\dots\}$ within distance $2^{-t/\alpha}$ of $g/|g|$ that carry $2^{t(2/3-5/(3\alpha))+o(t)}$ words each (with large constants: $126$ against about $5$ in the instance below)". Add "(heuristically)" to "therefore cannot give more than $\tfrac52$".

**F6 (should).** Quadric paragraph: "Auxiliary surfaces of higher degree do not obviously help. A rational quadric through the points of a longer box … recovers exactly the degree-one box count" → "Heuristically, auxiliary surfaces of higher degree do not help either: a rational quadric through the points of a box longer than a $\delta$-cube … recovers, in a model computation, exactly the degree-one box count $R'^{(8-2a_0)/3}$; we do not prove this." In the same Remark, write "the tube radius is $\delta=\varepsilon R'$ (not the error probability of~\eqref{eq:rate52})".

**F7 (should).** Theorem rate52: "all $t\le4L/(\tfrac85+\mu)$" → "all $t\le4(L+2)/(\tfrac85+\mu)$". In App H: "for $t\le4L/a_0$" → "for $t\le4(L+2)/a_0$ (then $\varepsilon\le2^{a_0+2}R'^{-a_0}$)", and add "this range contains $\tau_*=\lceil\alpha\log_2K\rceil$".

**F8 (nit).** "So $\tfrac52$ is approached and not attained." → "So the method approaches $\tfrac52$ but does not reach it." The paragraph's c and App H's c are different constants (c = 1/6 − c′/4); rename one of them.

**F9 (should).** In the main text: "projecting the numerators of the box along that plane or vector" → "projecting the numerators of the box onto that plane, or along that vector". Lemma G3(b) projects along V⊥ onto V.

**F10 (should).** Discussion: "Two further observations locate the difficulty (Remark~\ref{rem:what52}). Arguments that cover a tube by balls cannot pass $\tfrac52$, heuristically:" → "Two further observations, both heuristic, locate the difficulty. Arguments that cover a tube by balls cannot pass $\tfrac52$ (Remark~\ref{rem:what52}):". The Linnik sentence is not supported anywhere in the paper; leave it hedged as it is.

**F11 (nit).** Prop dense: "Then for $r=2$ there are" → "Then for $r=2$ and $Q=Q_\varepsilon$ (gaps measured cyclically in $\Z_{Q_\varepsilon}$) there are". The hypothesis that G is a Clifford+T word is not used; any G ∈ SU(2) works.

**F12 (should).** App G/H wording.
- **pm:** "Here $\mathrm{pm}$ is the offsets-plus-points-per-section exponent of Appendix~\ref{app:dioph}" → "… of Appendix~\ref{app:dioph} (offsets at $y_j$, points per section at $y_{j-1}$), maximised over the tangent and crossing classes".
- **\bar S:** "the patch exponents of Appendix~\ref{app:dioph}" → "the patch exponents of Lemma~\ref{lem:E2} at $g_i$".
- **Lemma G1 proof:** "some $\mathfrak n_i$ of a basis of $V^\perp\cap\mathcal O^*$ with" → "one of the $K$-independent $\mathfrak n_1,\mathfrak n_2\in V^\perp\cap\mathcal O^*$ of Lemma~\ref{lem:F5}(b), with".
- **Rank-three item:** add "where $W:=Kv_1+Kv_2+Kv_3=\mathfrak m^\perp$ contains $\mathcal K\cap\mathcal O^*$, hence every $n_B$ of the layer".
- **Finite instance:** "The largest exponent is" → "The largest exponent the certificate checks (it stops at the first admissible threshold) is".
- **App H, first line:** "Throughout, $E_1<T<a_0$, $a=a_0$, $b=E_1$" → "Throughout, $E_1:=(8-2a_0)/3<T<a_0$ ($T$ is the target $E_*$ of Appendices~\ref{app:dioph} and~\ref{app:proj}), $a=a_0$, $b=E_1$ up to the shift $\varpi$ below".
- **Margins:** after "with margin $\varpi$" add "(each such inequality only has to beat an absolute constant, which $R'^{\varpi}\to\infty$ does; this replaces the fixed margins of Appendix~\ref{app:elem})".
- **Last layer:** after "… on every layer below it" add "In the discretised assembly the last low layer $[y_{j^*-1},y_{j^*}]$ contains $\eta^*$, and its value is at most $\mathrm{LOW}(y_{j^*-1})+O(\varpi)\le T+O(\varpi)$."
- **Optional:** use y and y* in place of η and η*.

**F13 (should).** Data availability and README.
- **Data availability:** add "The example of Remark~\ref{rem:what52} is reproduced by \texttt{research/conditional3/ballpts.cpp} (\texttt{ballpts 25 5.828427 1 1 -1 0.0025}) and \texttt{research/conditional3/review/check\_sections.py}."
- **README:** add the matching row (output `review/bp_25_0.0024989.txt`, `review/known_output.txt`).
- **Git:** commit the three `research/` directories before release.
- **Short abstracts:** update `qip/short_abstract_*` unless they are frozen records.
