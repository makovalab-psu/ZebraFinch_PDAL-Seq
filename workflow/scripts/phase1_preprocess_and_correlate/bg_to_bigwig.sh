#!/usr/bin/env bash
#
# Convert a bedGraph to a bigWig.
#
# Adapted from workflow/scripts/01_preprocess_data/bg_to_bigwig.sh in the
# PDAL-Seq manuscript pipeline. The chrom.sizes path is passed in rather
# than parsed out of the sample name, and the bedGraph is sorted first.
#
# The sort is NOT optional. bedtools genomecov emits records in BAM header
# order, which for this assembly is the .fai order
# (chr1, chr1A, chr2, ... chr10, ...). bedGraphToBigWig requires ASCII
# order and aborts with "is not case-sensitive sorted" otherwise.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -b <in.bg> -c <chrom.sizes> -o <out.bigwig> [-d tmpdir]

  -b  Input bedGraph
  -c  UCSC chrom.sizes for the reference
  -o  Output bigWig
  -d  Directory for sort temp files (default: temp)
  -h  Show this message
EOF
}

BEDGRAPH=""
CHROM_SIZES=""
OUTPUT=""
TMP_DIR="temp"

while getopts "b:c:o:d:h" opt; do
    case "$opt" in
        b) BEDGRAPH="$OPTARG" ;;
        c) CHROM_SIZES="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "$BEDGRAPH" || -z "$CHROM_SIZES" || -z "$OUTPUT" ]]; then
    echo "ERROR: -b, -c and -o are all required" >&2
    usage
    exit 1
fi

mkdir -p "$TMP_DIR" "$(dirname "$OUTPUT")"

SORTED="${TMP_DIR}/$(basename "$OUTPUT" .bigwig)_$$.sorted.bg"
cleanup() {
    [[ -f "$SORTED" ]] && rm -f "$SORTED"
}
trap cleanup EXIT

LC_ALL=C sort -k1,1 -k2,2n -S 4G -T "$TMP_DIR" "$BEDGRAPH" > "$SORTED"

bedGraphToBigWig "$SORTED" "$CHROM_SIZES" "$OUTPUT"
