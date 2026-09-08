#!/usr/bin/env bash
#
# Histogram of CpG methylation rates inside one HMM state.
#
# Behaviour is the manuscript's
# workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh, which
# intersects the methylation bedGraph with the state and counts how many CpGs
# carry each rounded rate. Output is `count rate`, one line per rate, the same
# shape `sort | uniq -c` produces and the same shape the compile step reads.
#
# Two differences:
#
#   * `-u`. Without it a CpG that falls across the boundary of two segments of
#     the same state is counted twice. Ours is a 1 kbp grid, so boundaries are
#     everywhere.
#
#   * the histogram is built in a single awk pass rather than
#     `sort | uniq -c`. The Zebra finch methylation file is genome wide, so
#     sorting it once per state means ~20 M records x 77 states through sort
#     for a table with at most a thousand rows in it.
#
# One scale difference to be aware of downstream: the manuscript's methylation
# bedGraph held fractions (0-1) and this one holds percentages (0-100). The
# rates are printed as they are found; Compile_methylation.R bins on 0-100.
#
# Usage:
#   methylation_in_state.sh -s state.bed -m methylation.bed -o out.txt

set -euo pipefail

usage() {
    echo "Usage: $0 -s <state.bed> -m <methylation.bed> -o <out.txt>" >&2
    exit 1
}

STATE=""; METHYLATION=""; OUT=""

while getopts "s:m:o:h" opt; do
    case "$opt" in
        s) STATE="$OPTARG" ;;
        m) METHYLATION="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$STATE" && -n "$METHYLATION" && -n "$OUT" ]] || usage

mkdir -p "$(dirname "$OUT")"

# Staged and renamed for the same reason as Enrichment.sh: a job killed at the
# wall clock must leave no output rather than a half-written histogram.
STAGING="${OUT}.partial.$$"
cleanup() {
    rm -f "$STAGING" 2>/dev/null || true
}
trap cleanup EXIT

# An unoccupied state has no segments; an empty histogram is the honest answer
# and Compile_methylation.R reports it as zero CpGs.
if [[ ! -s "$STATE" ]]; then
    echo "WARNING: $STATE is empty -- this state is unoccupied in the fitted model." >&2
    : > "$OUT"
    exit 0
fi

bedtools intersect -a "$METHYLATION" -b "$STATE" -u \
    | awk '{ rate = sprintf("%.1f", $4); count[rate]++ }
           END { for (rate in count) printf "%d %s\n", count[rate], rate }' \
    > "$STAGING"

mv -f "$STAGING" "$OUT"

echo "$(wc -l < "$OUT" | tr -d "[:space:]") distinct methylation rates in $STATE"
