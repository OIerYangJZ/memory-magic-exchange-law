# arXiv replacement, v2

This is a **replacement** of the accepted v1 (17/7), not a new submission: on the arXiv user page choose *Replace* on arXiv:2609.37368, upload the new tarball and update abstract and comments. The arXiv ID stays the same; the new version is public only after it is announced.

Source: repository commit `8bab2d2`, built with `sh arxiv/make_arxiv.sh v2`.

## Upload
- `memory-magic-exchange-law-v2.tar.gz`: only the paper, `main.tex` and its 2 PNG figures (the bibliography is inline, no .bbl needed). Code and data are not uploaded; they are in the GitHub repository cited in the paper.
- No 00README: arXiv detects pdfLaTeX and the single top-level file itself and advises against a hand-written one.
- Compiler: pdfLaTeX. Top-level file: `main.tex`.
- `main_reference.pdf` is the local build (39 pages). arXiv's preview should match it; do **not** upload it.

## Metadata

**Title**
A Memory-Magic Exchange Law in Streaming Clifford+T Compilation

**Authors**
Jinze Yang, Yangyang Li, Xiu-Hao Deng

**Abstract** (1913 characters; limit 1920)
A phase that reaches a fault-tolerant processor in additive pieces can be remembered until the last piece arrives, or executed on arrival: the first option costs classical memory carried across rounds, the second costs magic states committed before the phase is known. We determine the exchange rate $\alpha$, committed $T$ gates per bit of memory forgone, for ancilla-free coordinatewise Clifford+$T$ compilation. The Ramanujan bound gives $\alpha\ge2$ with explicit constants, the square-root barrier of the spectral method. An elementary determinant method, using that quaternion numerators lie on spheres in both real embeddings of $\mathbb{Q}(\sqrt2)$, counts words near any rotation coset below that barrier. With a height dichotomy for the sphere sections it produces and a fibration over rational projections, it gives every $\alpha<5/2$ unconditionally and uniformly over frames, by an exact SMT check of its continuum case analysis; $5/2$ is where this method stops. At Clifford-framed cosets the volume law holds up to subexponential factors: all but a vanishing fraction of $z$-rotations need $T$-count $(3-o(1))\log_2(1/\varepsilon)$, and processes whose committed pieces are near Clifford-framed $z$-rotations, including per-rotation pipelines, have $\alpha\ge3-o(1)$, which a fractional-passthrough family attains under the Ross-Selinger typical-cost hypothesis. Under an equidistribution conjecture supported by exhaustive enumeration to $T$-count 22, $\alpha=3$ in general and memory should be shed in whole rotations. The bounds hold even when the phases cancel to the identity; side information enters through a conditional entropy; probabilistic mixing halves the costs but not the rate. With clean ancillas and a phase-gradient catalyst, table lookups batched over coordinates drive the rate to $O(1/\log\log(1/\varepsilon))$, so the constant-rate law is specific to coordinatewise synthesis.

**Comments**
39 pages, 2 figures. v2: unconditional rate improved from 17/7 to every alpha < 5/2 (Theorem 8, Appendices G-H, exact SMT check). Code, data and exact rational certificates: https://github.com/OIerYangJZ/memory-magic-exchange-law

**Primary category:** quant-ph
**Cross-list (optional):** math.NT
**MSC class (optional):** 81P68, 11P21, 11R52
**License:** CC BY 4.0, unchanged from v1 (the license of PRX Quantum and Quantum; Quantum requires it for the final arXiv version)

## Checks
- [ ] arXiv preview: 39 pages, both figures, affiliations and e-mail footnotes on page 1.
- [ ] Source commit pushed to `main` (the paper's data link points there).

## Timing (QIP 2027)
- Replacements received before 14:00 ET on a weekday are announced at 20:00 ET that day; those received Thu 14:00 ET to Fri 14:00 ET are announced Sunday 20:00 ET. Submit by Fri Oct 2 14:00 ET (Sat Oct 3 02:00 Beijing) at the latest, earlier to leave room for a moderation hold.
- The QIP talk submission (Oct 5 23:59 AoE) needs the versioned link `https://arxiv.org/abs/2609.37368v2` and, as technical manuscript, the v2 PDF exactly as arXiv serves it.

## After v2 is announced
- [ ] `git tag arxiv-v2 8bab2d2 && git push origin arxiv-v2`
- [ ] Check that `https://arxiv.org/abs/2609.37368v2` resolves and its PDF matches `main_reference.pdf`.
- [ ] QIP HotCRP: extended abstract (`qip/extended_abstract.pdf`), technical manuscript (the arXiv v2 PDF), link `https://arxiv.org/abs/2609.37368v2`.

tarball sha256: c49e23affab2ac434c1b611ba469df01dcd262a108359f908660f52a5bf34714
