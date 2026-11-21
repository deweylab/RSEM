#!/usr/bin/env bash
# Checks that rsem-generate-data-matrix rejects result files whose rows differ.
# Usage: tests/test_generate_data_matrix.sh   (run from the repository root)
set -u

script="$(cd "$(dirname "$0")/.." && pwd)/rsem-generate-data-matrix"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

header=$'gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count'
row() { printf '%s\tT\t1\t1\t%s\n' "$1" "$2"; }
{ echo "$header"; row G1 1.00; row G2 2.00; row G3 3.00; } > three.genes.results
{ echo "$header"; row G1 4.00; row G2 5.00; row G3 6.00; row G4 7.00; } > four.genes.results
{ echo "$header"; row G1 4.00; row G3 6.00; row G2 5.00; } > shuffled.genes.results
: > empty.genes.results

fail=0
expect_fail() {  # name, file args...
  local name="$1"; shift
  if perl "$script" "$@" > out.tsv 2> err.txt; then
    echo "FAIL: $name: exited 0 but the inputs differ"; fail=1
  else
    echo "ok:   $name: $(head -n 1 err.txt)"
  fi
}

perl "$script" three.genes.results three.genes.results > out.tsv 2> err.txt \
  && [[ "$(wc -l < out.tsv | tr -d ' ')" == 4 ]] \
  && echo "ok:   identical row sets" || { echo "FAIL: identical row sets"; fail=1; }

expect_fail "extra row in a later file" three.genes.results four.genes.results
expect_fail "fewer rows in a later file" four.genes.results three.genes.results
expect_fail "same rows in another order" three.genes.results shuffled.genes.results
expect_fail "empty file" three.genes.results empty.genes.results
expect_fail "missing file" three.genes.results missing.genes.results

exit "$fail"
