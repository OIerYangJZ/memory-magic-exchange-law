# arXiv submission, v1

Source: repository commit `bec7c63`, built with `sh arxiv/make_arxiv.sh v1`.

## Upload
- `memory-magic-exchange-law-v1.tar.gz`: only the paper, `main.tex` and its 2 PNG figures (the bibliography is inline, no .bbl needed). Code and data are not uploaded; they are in the GitHub repository cited in the paper.
- No 00README: arXiv detects pdfLaTeX and the single top-level file itself and advises against a hand-written one.
- Compiler: pdfLaTeX. Top-level file: `main.tex`.
- `main_reference.pdf` is the local build (36 pages). arXiv's preview should match it; do **not** upload it.

## Metadata

**Title**
A Memory-Magic Exchange Law in Streaming Clifford+T Compilation

**Authors**
Jinze Yang, Yangyang Li, Xiu-Hao Deng

**Abstract** (1860 characters; limit 1920)
A phase that reaches a fault-tolerant processor in additive pieces can be remembered until the last piece arrives, or executed on arrival: the first option costs classical memory carried across rounds, the second costs magic states committed before the phase is known. We determine the exchange rate $\alpha$, committed $T$ gates per bit of memory forgone, for ancilla-free coordinatewise Clifford+$T$ compilation. The Ramanujan bound for the Clifford+$T$ lattice gives $\alpha\ge2$ with explicit constants, the square-root barrier of the spectral method. An elementary determinant method, using that quaternion numerators are lattice points on spheres in both real embeddings of $\mathbb{Q}(\sqrt2)$, counts words near an arbitrary rotation coset below that barrier and gives $\alpha\ge11/5$ asymptotically and unconditionally, and a height dichotomy for the resulting sphere sections raises this to $\alpha\ge17/7$. At Clifford-framed cosets the volume law holds up to subexponential factors: all but a vanishing fraction of $z$-rotations need $T$-count $(3-o(1))\log_2(1/\varepsilon)$, and processes whose committed pieces are close to Clifford-framed $z$-rotations, including per-rotation pipelines, have $\alpha\ge3-o(1)$, which a fractional-passthrough family attains under the Ross-Selinger typical-cost hypothesis. Under an equidistribution conjecture supported by exhaustive enumeration to $T$-count 22, $\alpha=3$ in general and memory should be shed in whole rotations. The bounds hold even when the phases cancel to the identity; side information enters through a conditional entropy; probabilistic mixing halves the costs but not the rate. With clean ancillas and a phase-gradient catalyst, table lookups batched across coordinates drive the rate to $O(1/\log\log(1/\varepsilon))$, so the constant-rate law is specific to coordinatewise synthesis.

**Comments**
36 pages, 2 figures. Code, data and exact rational certificates: https://github.com/OIerYangJZ/memory-magic-exchange-law

**Primary category:** quant-ph
**Cross-list (optional):** math.NT
**MSC class (optional):** 81P68, 11P21, 11R52
**License:** CC BY 4.0 (the license of PRX Quantum and Quantum; Quantum requires it for the final arXiv version)

## Checks
- [ ] arXiv preview: 36 pages, both figures, affiliations and e-mail footnotes on page 1.
- [ ] Source commit pushed to `main` (the paper's data link points there).

## After the arXiv ID is assigned
- [ ] `git tag arxiv-v1 bec7c63 && git push origin arxiv-v1`
- [ ] Add the arXiv ID to README.md; link it in the QIP 2027 talk submission (due Oct 5).

tarball sha256: 439c6afc9349de2100f6af6b7eb8349dedaf167b7d845af9be360c64ad065e6f
