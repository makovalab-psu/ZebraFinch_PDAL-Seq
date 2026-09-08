#!/usr/bin/env bash
#
# Read enrichment of one sequencing library in one repeat class.
#
# Behaviour is the manuscript's
# workflow/scripts/17_PDALseq_reads_in_CenSat/Read_enrichment.sh:
#
#   repeat <iterations> times:
#       subsample <n> of the window-clipped annotation fragments
#       Observed = summed per-base read coverage over the subsample
#       shuffle the subsample to random coordinates on the same chromosome
#       Null     = summed per-base read coverage over the shuffle
#
# and enrichment is mean(Observed) / mean(Null), computed by the compile step.
#
# ONE substantive difference from the manuscript script: the coverage comes
# from the sample's bigWig instead of from its BAM.
#
#   megadepth reads a BAM front to back -- there is no index seek -- so one
#   iteration costs a full scan of a 5 GB alignment file. The loop makes two
#   megadepth calls, and Phase 5 has 498 of these jobs, so the manuscript's
#   version is ~100,000 whole-BAM scans and would not finish. A bigWig is
#   indexed, so the same query is a range read.
#
#   The number is the same one. Phase 1 built data/phase1/bigwig/<sample>.bigwig
#   from `bedtools genomecov -bg -split -ibam` on exactly the All_reads BAM
#   this phase is asking about, so `megadepth <bw> --op sum` and
#   `megadepth <bam> --op sum` are both "summed per-base read depth over these
#   intervals". It is a faster route to the statistic, not a different
#   statistic.
#
# Two mechanical differences, matching phase4_state_function/Enrichment.sh:
# temp files go to a per-job mktemp -d so hundreds of these can run at once,
# and the output is staged and renamed atomically so a job killed at the wall
# clock leaves nothing behind that could be mistaken for a finished result.
#
# Plain `bedtools sort` -- lexicographic, not .fai order -- as in the
# manuscript script. That is also the order the bigWig itself is in:
# bedGraphToBigWig refuses anything else, so Phase 1's bg_to_bigwig.sh sorted
# the bedGraph that way before building the track.
#
# Depth does not need normalising. Observed and Null are both sums over the
# same library, so the library size cancels in the ratio -- which is the whole
# point of running the PCR-free control through the identical procedure.
#
# Usage:
#   Read_enrichment.sh -c coverage.bigwig -a fragments.bed -g chrom.sizes \
#                      -o out.csv [-n subsamples] [-i iterations] [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -c <coverage.bigwig> -a <fragments.bed> -g <chrom.sizes> -o <out.csv> [-n <subsamples>] [-i <iterations>] [-d <tmpdir>]" >&2
    exit 1
}

COVERAGE=""; ANNOTATION=""; GENOME=""; OUT=""
SUBSAMPLES="1000"; ITERATIONS="100"; TMP_DIR="temp"

while getopts "c:a:g:o:n:i:d:h" opt; do
    case "$opt" in
        c) COVERAGE="$OPTARG" ;;
        a) ANNOTATION="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        n) SUBSAMPLES="$OPTARG" ;;
        i) ITERATIONS="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$COVERAGE" && -n "$ANNOTATION" && -n "$GENOME" && -n "$OUT" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")"

WORK="$(mktemp -d "${TMP_DIR}/read_enrichment_XXXXXXXX")"
STAGING="${OUT}.partial.$$"

cleanup() {
    rm -rf "$WORK" 2>/dev/null || true
    rm -f "$STAGING" 2>/dev/null || true
}
trap cleanup EXIT

N_ANNOTATION=$(wc -l < "$ANNOTATION" | tr -d "[:space:]")

# Every class in the Snakefile's lists has at least one annotation, and
# split_repeats.py already fails if one of them ends up empty, so an empty
# fragment file here means the window intersection produced nothing -- a real
# problem worth stopping for rather than papering over with a row of zeros.
if [[ "$N_ANNOTATION" -eq 0 ]]; then
    echo "ERROR: $ANNOTATION has no records." >&2
    exit 1
fi

SUBSAMPLE="${WORK}/subsample.bed"
SHUFFLED="${WORK}/shuffled.bed"

# megadepth prints the annotation with the summary statistic appended, one
# line per interval, so `$NF` is the sum whether the coverage source is a
# bigWig or a BAM. Empty stdout is never a legitimate answer -- an interval
# with no reads under it still prints a line with a 0 -- so it means megadepth
# failed, and the whole job should stop rather than report an enrichment of
# 0/0 for this class.
sum_coverage() {
    local bed="$1"
    local raw
    raw="$(megadepth "$COVERAGE" --annotation "$bed" --op sum)"
    if [[ -z "$raw" ]]; then
        echo "ERROR: megadepth returned nothing for $COVERAGE over $bed." >&2
        echo "       Check that the bigWig is readable and that its contig" >&2
        echo "       names match $GENOME." >&2
        exit 1
    fi
    printf '%s\n' "$raw" | awk '{ s += $NF } END { printf "%.0f\n", s + 0 }'
}

echo "Observed,Null,N_annotation" > "$STAGING"

for _ in $(seq 1 "$ITERATIONS"); do

    # shuf -n returns the whole file when it holds fewer than n records, which
    # is the methods' "all data were retained if the number of annotations was
    # <= the subsample size". Many satellite families have only a handful of
    # fragments, so this is the normal case here, not the exception -- see the
    # note on constant Observed in Compile_repeat_enrichment.R.
    shuf -n "$SUBSAMPLES" "$ANNOTATION" | bedtools sort > "$SUBSAMPLE"

    OBSERVED="$(sum_coverage "$SUBSAMPLE")"

    # -chrom keeps each fragment on the chromosome it came from, and shuffle
    # preserves interval length, so Observed and Null are sums over identical
    # length distributions on identical chromosomes.
    bedtools shuffle -i "$SUBSAMPLE" -g "$GENOME" -chrom \
        | bedtools sort > "$SHUFFLED"

    NULL="$(sum_coverage "$SHUFFLED")"

    echo "${OBSERVED},${NULL},${N_ANNOTATION}" >> "$STAGING"
done

# Only now does the output exist, and it is complete.
mv -f "$STAGING" "$OUT"
