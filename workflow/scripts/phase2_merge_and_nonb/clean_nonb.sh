#!/usr/bin/env bash
#
# Filter one non-B DNA motif annotation to the contigs that exist in our
# genome, reduce it to three columns, and sort it into .fai order.
#
# The Zenodo release annotates the FULL bTaeGut7v0.4 assembly, so every
# autosome is present twice (chr10_mat AND chr10_pat). Our genome is the
# single haplotype extracted by make_PDAL-Seq_fasta.sh: chr1-37 + chrW
# maternal, chrZ paternal, chrMT, 3 rDNA morphs. So the _pat autosomes are
# dropped here and roughly half of each file should survive.
#
# Sorting is by .fai RANK, not ASCII. The .fai order is chr1_mat,
# chr1A_mat, chr2_mat, ... which sorts neither lexicographically nor
# numerically, and `bedtools coverage -sorted -g` requires exactly this
# order. Same class of trap as the bedGraphToBigWig sort in Phase 1.
#
# Usage:
#   clean_nonb.sh -i in.bed -g chrom.sizes -o out.bed -r report.txt \
#                 -n CLASS [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -i <in.bed> -g <chrom.sizes> -o <out.bed> -r <report.txt> -n <class> [-d <tmpdir>]" >&2
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

# Rank every contig by its position in chrom.sizes, keep only records on a
# contig that has a rank, and sort by (rank, start).
LC_ALL=C awk '
    NR == FNR { rank[$1] = FNR; next }
    ($1 in rank) { print rank[$1] "\t" $1 "\t" $2 "\t" $3 }
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
