#!/usr/bin/env bash
#
# Enrichment of a discrete annotation (BED) in one HMM state.
#
# Behaviour is the manuscript's
# workflow/scripts/07_PhyloHGMP_model/Enrichment.sh:
#
#   repeat <iterations> times:
#       subsample <n> annotations at random
#       Observed = bp of intersection with the state
#       shuffle the subsample to random coordinates on the same chromosome
#       Null     = bp of intersection with the state
#
# and enrichment is mean(Observed) / mean(Null), computed by the compile step.
#
# Differences from the manuscript script, all mechanical:
#
#   * the chrom.sizes path is an argument instead of the hardcoded
#     resources/genomes/Hsapien.chrom.sizes;
#   * temp files go to a per-job mktemp -d under a caller-supplied temp dir,
#     so 1200 of these can run concurrently without colliding on
#     "$output.temp";
#   * awk prints `sum + 0`, so a zero-overlap iteration writes 0 rather than
#     an empty field that would break the R compile step;
#   * an unoccupied state short circuits to a table of zeros instead of doing
#     100 iterations of work against an empty file;
#   * the record count of the annotation is carried in a third column. The
#     methods pick the significance test on that number (Wilcoxon above
#     10,000 loci, empirical below), so it has to reach the compile step.
#
# The shuffle is genome wide and NOT restricted to the segmented genome, which
# is the manuscript's behaviour. See README section 06 for what that means for
# reading the ratios.
#
# Usage:
#   Enrichment.sh -s state.bed -a annotation.bed -g chrom.sizes -o out.csv \
#                 [-n subsamples] [-i iterations] [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -s <state.bed> -a <annotation.bed> -g <chrom.sizes> -o <out.csv> [-n <subsamples>] [-i <iterations>] [-d <tmpdir>]" >&2
    exit 1
}

STATE=""; ANNOTATION=""; GENOME=""; OUT=""
SUBSAMPLES="10000"; ITERATIONS="100"; TMP_DIR="temp"

while getopts "s:a:g:o:n:i:d:h" opt; do
    case "$opt" in
        s) STATE="$OPTARG" ;;
        a) ANNOTATION="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        n) SUBSAMPLES="$OPTARG" ;;
        i) ITERATIONS="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$STATE" && -n "$ANNOTATION" && -n "$GENOME" && -n "$OUT" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")"

WORK="$(mktemp -d "${TMP_DIR}/enrichment_XXXXXXXX")"

# Rows are appended one iteration at a time, so a job killed at the SLURM wall
# clock would leave a short file behind. Snakemake only re-runs an output it
# recorded as incomplete, and a SIGKILLed snakemake never gets to record
# anything -- the truncated CSV would then be read as a finished result with
# fewer iterations in it. Staging next to the final path (same directory, so
# the rename is atomic) means a killed job leaves nothing behind and the next
# run simply redoes the job. This matters because Phase 4 is expected to be
# restarted several times against a billing-minutes budget.
STAGING="${OUT}.partial.$$"

cleanup() {
    rm -rf "$WORK" 2>/dev/null || true
    rm -f "$STAGING" 2>/dev/null || true
}
trap cleanup EXIT

N_ANNOTATION=$(wc -l < "$ANNOTATION" | tr -d "[:space:]")
N_STATE=$(wc -l < "$STATE" | tr -d "[:space:]")

if [[ "$N_ANNOTATION" -eq 0 ]]; then
    echo "ERROR: $ANNOTATION has no records." >&2
    exit 1
fi

echo "Observed,Null,N_annotation" > "$STAGING"

# An unoccupied state has no segments, so every intersection is empty by
# construction. Writing the zeros keeps the output shape uniform.
if [[ "$N_STATE" -eq 0 ]]; then
    echo "WARNING: $STATE is empty -- this state is unoccupied in the fitted model." >&2
    for _ in $(seq 1 "$ITERATIONS"); do
        echo "0,0,${N_ANNOTATION}" >> "$STAGING"
    done
    mv -f "$STAGING" "$OUT"
    exit 0
fi

SUBSAMPLE="${WORK}/subsample.bed"
SHUFFLED="${WORK}/shuffled.bed"

for _ in $(seq 1 "$ITERATIONS"); do

    # shuf -n returns the whole file when it holds fewer than n records, which
    # is the methods' "all data were retained if the number of annotations was
    # <= 10,000". No special case needed.
    shuf -n "$SUBSAMPLES" "$ANNOTATION" > "$SUBSAMPLE"

    OBSERVED=$(bedtools intersect -a "$STATE" -b "$SUBSAMPLE" \
        | awk '{ sum += $3 - $2 + 1 } END { print sum + 0 }')

    # -chrom keeps each annotation on the chromosome it came from, so the null
    # controls for chromosome-level composition.
    bedtools shuffle -i "$SUBSAMPLE" -g "$GENOME" -chrom > "$SHUFFLED"

    NULL=$(bedtools intersect -a "$STATE" -b "$SHUFFLED" \
        | awk '{ sum += $3 - $2 + 1 } END { print sum + 0 }')

    echo "${OBSERVED},${NULL},${N_ANNOTATION}" >> "$STAGING"
done

# Only now does the output exist, and it is complete.
mv -f "$STAGING" "$OUT"
