#!/usr/bin/env bash
#
# Enrichment of a continuous signal (bedGraph) in one HMM state.
#
# Behaviour is the manuscript's
# workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh: subsample the
# bedGraph, sum column 4 over the records that hit the state, shuffle the
# subsample, sum again. Enrichment is mean(Observed) / mean(Null).
#
# Differences from the manuscript script:
#
#   * chrom.sizes is an argument, temp files are per job, awk prints `sum + 0`,
#     an unoccupied state short circuits, and the record count is carried in a
#     third column -- all as in Enrichment.sh, for the same reasons;
#
#   * `bedtools intersect -u`. Without it a bedGraph record straddling two
#     segments of the SAME state is emitted twice and its value counted twice.
#     The manuscript's input was a raw genomecov bedGraph whose records are a
#     few bp long, so straddling was rare; ours is on the 1 kbp window grid the
#     states are built from, where a shuffled window landing across a 1 kbp gap
#     between two segments of one state is common. -u makes it "the sum of the
#     observations that fall in this state", counted once each.
#
# The input bedGraph is RNA-Seq coverage summed over the same 1 kbp windows the
# HMM was fit to, not a raw `bedtools genomecov -bg`. A genomecov bedGraph for
# a 464 M pair library is on the order of 10^8 records and `shuf` reads its
# input once per iteration; at 100 iterations x 154 jobs that is a petabyte of
# I/O for a signal that is being read at 1 kbp resolution anyway.
#
# Usage:
#   BG_enrichment.sh -s state.bed -a signal.bg -g chrom.sizes -o out.csv \
#                    [-n subsamples] [-i iterations] [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -s <state.bed> -a <signal.bg> -g <chrom.sizes> -o <out.csv> [-n <subsamples>] [-i <iterations>] [-d <tmpdir>]" >&2
    exit 1
}

STATE=""; SIGNAL=""; GENOME=""; OUT=""
SUBSAMPLES="100000"; ITERATIONS="100"; TMP_DIR="temp"

while getopts "s:a:g:o:n:i:d:h" opt; do
    case "$opt" in
        s) STATE="$OPTARG" ;;
        a) SIGNAL="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        n) SUBSAMPLES="$OPTARG" ;;
        i) ITERATIONS="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$STATE" && -n "$SIGNAL" && -n "$GENOME" && -n "$OUT" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")"

WORK="$(mktemp -d "${TMP_DIR}/bg_enrichment_XXXXXXXX")"

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

N_SIGNAL=$(wc -l < "$SIGNAL" | tr -d "[:space:]")
N_STATE=$(wc -l < "$STATE" | tr -d "[:space:]")

if [[ "$N_SIGNAL" -eq 0 ]]; then
    echo "ERROR: $SIGNAL has no records." >&2
    exit 1
fi

echo "Observed,Null,N_annotation" > "$STAGING"

if [[ "$N_STATE" -eq 0 ]]; then
    echo "WARNING: $STATE is empty -- this state is unoccupied in the fitted model." >&2
    for _ in $(seq 1 "$ITERATIONS"); do
        echo "0,0,${N_SIGNAL}" >> "$STAGING"
    done
    mv -f "$STAGING" "$OUT"
    exit 0
fi

SUBSAMPLE="${WORK}/subsample.bg"
SHUFFLED="${WORK}/shuffled.bg"

for _ in $(seq 1 "$ITERATIONS"); do

    shuf -n "$SUBSAMPLES" "$SIGNAL" > "$SUBSAMPLE"

    OBSERVED=$(bedtools intersect -a "$SUBSAMPLE" -b "$STATE" -u \
        | awk '{ sum += $4 } END { print sum + 0 }')

    bedtools shuffle -i "$SUBSAMPLE" -g "$GENOME" -chrom > "$SHUFFLED"

    NULL=$(bedtools intersect -a "$SHUFFLED" -b "$STATE" -u \
        | awk '{ sum += $4 } END { print sum + 0 }')

    echo "${OBSERVED},${NULL},${N_SIGNAL}" >> "$STAGING"
done

# Only now does the output exist, and it is complete.
mv -f "$STAGING" "$OUT"
