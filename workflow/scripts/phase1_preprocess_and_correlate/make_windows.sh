#!/usr/bin/env bash
#
# Make a BED file of evenly spaced genome windows.
# Verbatim behaviour from
# workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh
# in the PDAL-Seq manuscript pipeline, including its $2+1 start shift.
#
# The +1 is kept deliberately: it means window coordinates in this project
# line up with every window coordinate in the manuscript, at the cost of
# skipping the first base of each window. Do not "fix" it without
# regenerating every downstream coverage file.
#
# $1 a UCSC chrom.sizes file
# $2 the output bed file
# $3 the window size in base pairs

set -euo pipefail

mkdir -p "$(dirname "$2")"

bedtools makewindows -g "$1" -w "$3" \
    | cut -f 1,2,3 \
    | awk 'BEGIN {OFS="\t"} {print $1, $2 + 1, $3}' > "$2"
