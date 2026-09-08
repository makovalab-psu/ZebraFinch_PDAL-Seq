#!/usr/bin/env bash
#
# Concatenate the cleaned non-B motif classes into the "all" set.
#
# Concatenated, NOT merged -- this reproduces the manuscript's all.sh.
# The classes overlap one another (and six of the seven self-overlap, per
# the merge tests in resources/README.md), so motif COUNT in "all" double
# counts. That is why the compile step also reports base_density, which is
# computed from bases covered and is unaffected by overlap.
#
# Usage: combine_nonb.sh -g chrom.sizes -o all.bed [-d tmpdir] in1.bed in2.bed ...

set -euo pipefail

usage() {
    echo "Usage: $0 -g <chrom.sizes> -o <out.bed> [-d <tmpdir>] <in.bed> [in.bed ...]" >&2
    exit 1
}

GENOME=""; OUT=""; TMP_DIR="temp"

while getopts "g:o:d:h" opt; do
    case "$opt" in
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done
shift $((OPTIND - 1))

[[ -n "$GENOME" && -n "$OUT" && $# -gt 0 ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")"

# Same rank sort as clean_nonb.sh -- the result feeds bedtools coverage
# -sorted -g, which demands .fai order.
LC_ALL=C awk '
    NR == FNR { rank[$1] = FNR; next }
    ($1 in rank) { print rank[$1] "\t" $1 "\t" $2 "\t" $3 }
' "$GENOME" "$@" \
    | LC_ALL=C sort -k1,1n -k3,3n -S 4G -T "$TMP_DIR" \
    | cut -f 2- \
    > "$OUT"

echo "all: $(wc -l < "$OUT" | tr -d "[:space:]") records from $# classes"
