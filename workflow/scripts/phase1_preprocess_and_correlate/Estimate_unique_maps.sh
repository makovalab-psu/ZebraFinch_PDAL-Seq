#!/usr/bin/env bash
#
# Histogram of (primary alignment score - suboptimal alignment score),
# i.e. AS:i minus XS:i. A large positive difference means the read has one
# clearly best placement; 0 means the aligner found an equally good
# alternative.
#
# Adapted from
# workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh
# in the PDAL-Seq manuscript pipeline.
#
# DEVIATION FROM THE MANUSCRIPT: the original awk never reset as_value /
# xs_value between records, so a read with no XS tag (i.e. a uniquely
# placed read, which is exactly the case being counted) silently inherited
# the previous read's XS. Both are reset to 0 for each record here, so a
# read with AS and no XS scores AS - 0 and lands in the unique bin where
# it belongs. Counts from this script are therefore not directly
# comparable to the manuscript's Table S3.
#
# $1 a mapped bam file
# $2 sort temp directory (optional, default: temp)

set -euo pipefail

TMP_DIR="${2:-temp}"
mkdir -p "$TMP_DIR"

echo "count Primary.score.minus.suplemental.score"

samtools view "$1" \
    | awk 'BEGIN {OFS="\t"}
        {
            as_value = 0
            xs_value = 0
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^AS:i:/) as_value = substr($i, 6) + 0
                if ($i ~ /^XS:i:/) xs_value = substr($i, 6) + 0
            }
            print as_value - xs_value
        }' \
    | sort -T "$TMP_DIR" -n \
    | uniq -c
