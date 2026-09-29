#!/bin/sh
# Build the arXiv upload package from the repository root:  sh arxiv/make_arxiv.sh v1
# Copies the source and figures, compiles them in isolation with plain pdflatex (as arXiv does),
# checks the log, and writes arxiv/<tag>/memory-magic-exchange-law-<tag>.tar.gz.
set -e
TAG=${1:-v1}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DIR="$ROOT/arxiv/$TAG"
rm -rf "$DIR"; mkdir -p "$DIR/src" "$DIR/build"
for f in main.tex fig1_tube_counts.png fig2_frontier_twopanel.png; do cp "$ROOT/$f" "$DIR/src/"; done
cp "$DIR/src/"* "$DIR/build/"
cd "$DIR/build"
for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error main.tex > /dev/null; done
if grep -E "^!|undefined|Rerun to get|multiply defined" main.log; then echo "LaTeX problems, see $DIR/build/main.log"; exit 1; fi
cp main.pdf "$DIR/main_reference.pdf"
# ancillary files (arXiv serves anc/ with the paper): code, data and exact certificates
A="$DIR/src/anc"; mkdir -p "$A/scripts" "$A/data" "$A/research"
cp "$ROOT/arxiv/anc_README.md" "$A/README.md"; cp "$ROOT/requirements.txt" "$A/"
cp "$ROOT"/scripts/*.py "$A/scripts/"
for f in "$ROOT"/data/*; do case "$f" in *.bak-*) ;; *) cp "$f" "$A/data/";; esac; done
for d in cert dioph verify4; do mkdir -p "$A/research/$d"; cp "$ROOT/research/$d"/*.py "$ROOT/research/$d"/*.txt "$ROOT/research/$d"/*.md "$A/research/$d/" 2>/dev/null || true; done
if find "$A" -name "*.tex" | grep -q .; then echo "TeX files in anc/"; exit 1; fi
cd "$DIR/src" && COPYFILE_DISABLE=1 tar --uid 0 --gid 0 --uname arxiv --gname arxiv -czf "$DIR/memory-magic-exchange-law-$TAG.tar.gz" main.tex fig1_tube_counts.png fig2_frontier_twopanel.png anc
cd "$DIR" && rm -rf build
if grep -q "Funding and acknowledgments to be supplied" "$DIR/src/main.tex"; then
  echo "WARNING: the acknowledgments placeholder is still in main.tex"
fi
echo "package: $DIR/memory-magic-exchange-law-$TAG.tar.gz"
echo "pages:   $(pdfinfo "$DIR/main_reference.pdf" 2>/dev/null | awk '/^Pages/{print $2}')"

# ---- submission sheet
cd "$ROOT"
PAGES=$(pdfinfo "$DIR/main_reference.pdf" 2>/dev/null | awk '/^Pages/{print $2}')
ABS=$(cat "$ROOT/qip/arxiv_abstract.txt" 2>/dev/null || echo "(qip/arxiv_abstract.txt not found)")
COMMIT=$(git rev-parse --short HEAD)
DIRTY=$(git status --porcelain -- main.tex fig1_tube_counts.png fig2_frontier_twopanel.png | wc -l | tr -d ' ')
cat > "$DIR/SUBMISSION.md" <<EOT
# arXiv submission, $TAG

Source: repository commit \`$COMMIT\`$( [ "$DIRTY" != "0" ] && echo " **plus uncommitted changes to the paper sources**" ), built with \`sh arxiv/make_arxiv.sh $TAG\`.

## Upload
- \`memory-magic-exchange-law-$TAG.tar.gz\`: \`main.tex\` and the 2 PNG figures at the root (the bibliography is inline, no .bbl needed), plus the ancillary directory \`anc/\` (scripts, data, exact certificates; see \`anc/README.md\`).
- No 00README: arXiv detects pdfLaTeX and the single top-level file itself and advises against a hand-written one.
- Compiler: pdfLaTeX. Top-level file: \`main.tex\`.
- \`main_reference.pdf\` is the local build ($PAGES pages). arXiv's preview should match it; do **not** upload it.

## Metadata

**Title**
A Memory-Magic Exchange Law in Streaming Clifford+T Compilation

**Authors**
Jinze Yang, Yangyang Li, Xiu-Hao Deng

**Abstract** ($(printf %s "$ABS" | wc -c | tr -d ' ') characters; limit 1920)
$ABS

**Comments**
$PAGES pages, 2 figures. Code, data and exact rational certificates are included as ancillary files and at https://github.com/OIerYangJZ/memory-magic-exchange-law

**Primary category:** quant-ph
**Cross-list (optional):** math.NT
**MSC class (optional):** 81P68, 11P21, 11R52
**License:** CC BY 4.0 (the license of PRX Quantum and Quantum; Quantum requires it for the final arXiv version)

## Checks
- [ ] arXiv preview: $PAGES pages, both figures, affiliations and e-mail footnotes on page 1.
- [ ] The file list shows \`anc/\` (arXiv lists ancillary files separately on the abstract page).
- [ ] Source commit pushed to \`main\` (the paper's data link points there).

## After the arXiv ID is assigned
- [ ] \`git tag arxiv-$TAG $COMMIT && git push origin arxiv-$TAG\`
- [ ] Add the arXiv ID to README.md; link it in the QIP 2027 talk submission (due Oct 5).

tarball sha256: $(shasum -a 256 "$DIR/memory-magic-exchange-law-$TAG.tar.gz" | cut -d' ' -f1)
EOT
echo "sheet:   $DIR/SUBMISSION.md"
