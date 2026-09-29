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
cd "$DIR/src" && tar --uid 0 --gid 0 --uname arxiv --gname arxiv -czf "$DIR/memory-magic-exchange-law-$TAG.tar.gz" main.tex fig1_tube_counts.png fig2_frontier_twopanel.png
cd "$DIR" && rm -rf build
if grep -q "Funding and acknowledgments to be supplied" "$DIR/src/main.tex"; then
  echo "WARNING: the acknowledgments placeholder is still in main.tex"
fi
echo "package: $DIR/memory-magic-exchange-law-$TAG.tar.gz"
echo "pages:   $(pdfinfo "$DIR/main_reference.pdf" 2>/dev/null | awk '/^Pages/{print $2}')"
