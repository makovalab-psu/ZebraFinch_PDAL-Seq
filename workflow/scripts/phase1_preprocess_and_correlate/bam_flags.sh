#!/usr/bin/env bash
#
# Count reads carrying each SAM flag.
# Adapted from workflow/scripts/02_preprocessing_statistics/bam_flags.sh
# in the PDAL-Seq manuscript pipeline; the only change is an explicit
# sort temp directory so large BAMs do not fill /tmp on a compute node.
#
# $1 a mapped bam file
# $2 output: column 1 = read count, column 2 = SAM flag
# $3 sort temp directory (optional, default: temp)
#
# See https://broadinstitute.github.io/picard/explain-flags.html

set -euo pipefail

TMP_DIR="${3:-temp}"
mkdir -p "$TMP_DIR" "$(dirname "$2")"

samtools view "$1" | cut -f 2 | sort -T "$TMP_DIR" | uniq -c > "$2"
