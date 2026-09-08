#!/usr/bin/env bash
#
# Filter one functional annotation to the contigs in our genome, reduce it to
# three columns, clamp it to the contig bounds, and sort it into .fai order.
#
# Same job as phase2_merge_and_nonb/clean_nonb.sh, with one addition: the
# coordinates are clamped to the contig length and zero-length records are
# dropped. `bedtools shuffle` aborts with
#
#   Error: Interval ... is larger than the length of chromosome ...
#
# on a single out-of-range record, which would take out one of ~600 enrichment
# jobs for a reason that has nothing to do with the analysis. These annotations
# were derived from a GFF against Linnea's matZ extraction rather than against
# our make_PDAL-Seq_fasta.sh extraction, so an off-by-one at a contig end is
# cheap to guard against and expensive to debug at hour 30 of a run.
#
# Sorting is by .fai RANK, not ASCII: the .fai order is chr1_mat, chr1A_mat,
# chr2_mat, ... which is neither lexicographic nor numeric.
#
# Usage:
#   clean_annotation.sh -i in.bed -g chrom.sizes -o out.bed -r report.txt \
#                       -n NAME [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -i <in.bed> -g <chrom.sizes> -o <out.bed> -r <report.txt> -n <name> [-d <tmpdir>]" >&2
    exit 1
}

IN=""; GENOME=""; OUT=""; REPORT=""; NAME=""; TMP_DIR="temp"

while getopts "i:g:o:r:n:d:h" opt; do
    case "$opt" in
        i) IN="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        r) REPORT="$OPTARG" ;;
        n) NAME="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$IN" && -n "$GENOME" && -n "$OUT" && -n "$REPORT" && -n "$NAME" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")" "$(dirname "$REPORT")"

# NR == FNR reads chrom.sizes: rank[] for the sort key, len[] for the clamp.
# The numeric guard on $2/$3 skips a header or track line rather than letting
# awk coerce it to 0 and emit a record at the start of the contig.
LC_ALL=C awk '
    NR == FNR { rank[$1] = FNR; len[$1] = $2 + 0; next }
    !($1 in rank) { next }
    ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) { next }
    {
        s = $2 + 0; e = $3 + 0;
        if (s < 0) s = 0;
        if (e > len[$1]) e = len[$1];
        if (e > s) print rank[$1] "\t" $1 "\t" s "\t" e;
    }
' "$GENOME" "$IN" \
    | LC_ALL=C sort -k1,1n -k3,3n -S 2G -T "$TMP_DIR" \
    | cut -f 2- \
    > "$OUT"

IN_N=$(wc -l < "$IN" | tr -d "[:space:]")
OUT_N=$(wc -l < "$OUT" | tr -d "[:space:]")

if [[ "$OUT_N" -eq 0 ]]; then
    echo "ERROR: no records in $IN survived the contig filter." >&2
    echo "       Not one chromosome name matched $GENOME." >&2
    echo "       Chromosome names in the annotation:" >&2
    cut -f 1 "$IN" | LC_ALL=C sort -u | head -20 >&2
    echo "       Chromosome names in the genome:" >&2
    cut -f 1 "$GENOME" | head -20 >&2
    exit 1
fi

CONTIGS=$(cut -f 1 "$OUT" | LC_ALL=C sort -u | wc -l | tr -d "[:space:]")

printf '%s\t%s\t%s\t%s\t%s\n' \
    "$NAME" "$IN_N" "$OUT_N" \
    "$(awk -v a="$OUT_N" -v b="$IN_N" 'BEGIN { printf "%.1f", 100 * a / b }')" \
    "$CONTIGS" \
    > "$REPORT"

echo "$NAME: kept $OUT_N of $IN_N records on $CONTIGS contigs"
