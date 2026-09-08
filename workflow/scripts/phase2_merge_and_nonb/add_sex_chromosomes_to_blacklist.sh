#!/usr/bin/env bash
#
# Add the sex chromosomes to the blacklist in their entirety, then put the
# file back into .fai order.
#
# Adapted from 03_.../add_sex_chromosomes_to_blacklist.sh, which hard-codes
# chrX and chrY. The Zebra finch assembly is ZW: chrW_mat and chrZ_pat.
# Any partial sex-chromosome intervals the mappability rule produced are
# dropped first so they are not duplicated by the whole-chromosome entry.
#
# Usage:
#   add_sex_chromosomes_to_blacklist.sh -b merged.bed -g chrom.sizes \
#       -s chrW_mat,chrZ_pat -o out.bed [-d tmpdir]

set -euo pipefail

usage() {
    echo "Usage: $0 -b <blacklist.bed> -g <chrom.sizes> -s <chr1,chr2> -o <out.bed> [-d <tmpdir>]" >&2
    exit 1
}

BLACKLIST=""; GENOME=""; SEX=""; OUT=""; TMP_DIR="temp"

while getopts "b:g:s:o:d:h" opt; do
    case "$opt" in
        b) BLACKLIST="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        s) SEX="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$BLACKLIST" && -n "$GENOME" && -n "$SEX" && -n "$OUT" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")"

COMBINED="${TMP_DIR}/blacklist_with_sex_$$.bed"
trap 'rm -f "$COMBINED"' EXIT

# Records on autosomes, plus one full-length record per sex chromosome.
# awk rather than `grep -v`, which exits 1 when it filters everything out
# and would kill the script under `set -e`.
LC_ALL=C awk -v sex="$SEX" '
    BEGIN {
        n = split(sex, a, ",")
        for (i = 1; i <= n; i++) drop[a[i]] = 1
    }
    !($1 in drop) { print $1 "\t" $2 "\t" $3 }
' "$BLACKLIST" > "$COMBINED"

LC_ALL=C awk -v sex="$SEX" '
    BEGIN {
        n = split(sex, a, ",")
        for (i = 1; i <= n; i++) want[a[i]] = 1
    }
    ($1 in want) { print $1 "\t0\t" $2 }
' "$GENOME" >> "$COMBINED"

LC_ALL=C awk '
    NR == FNR { rank[$1] = FNR; next }
    ($1 in rank) { print rank[$1] "\t" $1 "\t" $2 "\t" $3 }
' "$GENOME" "$COMBINED" \
    | LC_ALL=C sort -k1,1n -k3,3n -S 2G -T "$TMP_DIR" \
    | cut -f 2- \
    > "$OUT"

echo "blacklist: $(wc -l < "$OUT" | tr -d "[:space:]") intervals covering $(awk '{ s += $3 - $2 } END { printf "%.1f", s / 1e6 }' "$OUT") Mbp"
