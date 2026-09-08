#!/usr/bin/env bash
#
# Count reads in a fastq.gz.
# Verbatim behaviour from
# workflow/scripts/02_preprocessing_statistics/read_count_fastq.gz.sh
# in the PDAL-Seq manuscript pipeline.
#
# $1 the fastq.gz to count
# $2 the text file to write the answer to

set -euo pipefail

echo $(( $(zcat -f "$1" | wc -l) / 4 )) > "$2"
