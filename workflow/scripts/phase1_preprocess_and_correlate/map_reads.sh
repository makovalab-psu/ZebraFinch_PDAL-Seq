#!/usr/bin/env bash
#
# Map single-end (R1) reads to the reference genome.
#
# Adapted from workflow/scripts/01_preprocess_data/map_reads.sh in the
# PDAL-Seq manuscript pipeline. Two changes:
#   * The genome is passed in explicitly. The manuscript version derived
#     it from the first underscore-delimited field of the sample name,
#     which has nothing to infer here (one genome, and the fastq stems
#     start with "Tcas_" regardless of what they map to).
#   * bwa mem streams into samtools sort instead of writing a temp SAM.
#     At 119 M reads a SAM is ~40 GB on scratch for no benefit.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -f <fastq.gz> -g <genome.fa> -o <out.bam> [-t threads] [-d tmpdir]

  -f  Trimmed single-end fastq.gz (R1)
  -g  Reference fasta, with bwa indices alongside it
  -o  Output coordinate-sorted BAM
  -t  Threads for bwa mem (default 1)
  -d  Directory for samtools sort temp files (default: temp)
  -h  Show this message
EOF
}

FASTQ=""
GENOME=""
OUTPUT=""
THREADS=1
TMP_DIR="temp"

while getopts "f:g:o:t:d:h" opt; do
    case "$opt" in
        f) FASTQ="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "$FASTQ" || -z "$GENOME" || -z "$OUTPUT" ]]; then
    echo "ERROR: -f, -g and -o are all required" >&2
    usage
    exit 1
fi

mkdir -p "$TMP_DIR" "$(dirname "$OUTPUT")"

# Sort threads are kept well below bwa's so the two stages do not fight
# over the cores snakemake allocated to this job.
SORT_THREADS=4
if [[ "$THREADS" -lt 4 ]]; then
    SORT_THREADS="$THREADS"
fi

SORT_PREFIX="${TMP_DIR}/sort_$(basename "$OUTPUT" .bam)_$$"

echo "genome:  $GENOME"
echo "fastq:   $FASTQ"
echo "output:  $OUTPUT"
echo "threads: bwa=${THREADS} sort=${SORT_THREADS}"

bwa mem -t "$THREADS" "$GENOME" "$FASTQ" \
    | samtools sort -@ "$SORT_THREADS" -m 2G -T "$SORT_PREFIX" -O bam -o "$OUTPUT"
