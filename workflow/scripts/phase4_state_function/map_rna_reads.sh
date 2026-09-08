#!/usr/bin/env bash
#
# Splice-aware alignment of one paired RNA-Seq library.
#
# This is the mapping half of the manuscript's 18_ATACseq_analysis.smk with
# bwa mem replaced by STAR. ATAC-Seq fragments are genomic, so bwa mem is
# correct there; RNA-Seq reads cross splice junctions and bwa mem would
# soft-clip every one of them.
#
# STAR streams unsorted BAM on stdout straight into samtools sort. The
# manuscript writes a SAM to disk between the two steps (ATACseq_bwa_mem ->
# ATACseq_sort_sam); at 464 M read pairs that intermediate is several hundred
# GB for no benefit. Phase 1's map_reads.sh already made the same call.
#
# Junctions are discovered de novo unless -a is given. resources/ holds the
# BED files Linnea derived from the annotation but not the GFF they came from,
# so there is nothing to hand --sjdbGTFfile by default. At the 1 kbp window
# resolution this signal is read at, annotation-guided junctions do not move
# the numbers.
#
# MAPQ note for anyone reading these BAMs later: STAR encodes 255 = unique,
# 3 = two loci, 1 = three or four, 0 = more. That is NOT bwa's 0-60 scale, so
# the `--min-MQ 20` idiom from Phase 1 means something different here.
#
# Memory: STAR holds the index (~12-14 GB for this genome) for the whole run
# and samtools sort runs concurrently in the same pipe. The sort is therefore
# pinned to 6 threads x 1 GB rather than inheriting -t, which keeps the pair
# inside the ~28 GB the calling rule reserves. Raising -t alone will not
# overrun that; raising the sort budget below will.
#
# Usage:
#   map_rna_reads.sh -1 R1.fastq.gz -2 R2.fastq.gz -x star_index_dir \
#                    -o out.bam -p log_prefix -t threads [-d tmpdir] [-a ann.gff]

set -euo pipefail

usage() {
    echo "Usage: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -x <star_index_dir> -o <out.bam> -p <log_prefix> -t <threads> [-d <tmpdir>] [-a <annotation.gff>]" >&2
    exit 1
}

R1=""; R2=""; INDEX=""; OUT=""; PREFIX=""; THREADS="1"; TMP_DIR="temp"; ANNOTATION=""

while getopts "1:2:x:o:p:t:d:a:h" opt; do
    case "$opt" in
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        x) INDEX="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        a) ANNOTATION="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -n "$R1" && -n "$R2" && -n "$INDEX" && -n "$OUT" && -n "$PREFIX" ]] || usage

mkdir -p "$TMP_DIR" "$(dirname "$OUT")" "$(dirname "$PREFIX")"

SORT_THREADS=6
SORT_MEMORY="1G"

# An empty array expanded under `set -u` is an unbound variable on older bash,
# hence the ${x[@]+"${x[@]}"} guard rather than a bare "${STAR_ARGS[@]}".
STAR_ARGS=()
if [[ -n "$ANNOTATION" ]]; then
    # --sjdbGTFtagExonParentTranscript Parent is what makes STAR read GFF3
    # rather than GTF attribute names.
    STAR_ARGS+=(--sjdbGTFfile "$ANNOTATION" --sjdbGTFtagExonParentTranscript Parent)
    echo "using annotation-guided junctions from $ANNOTATION"
else
    echo "discovering junctions de novo (no annotation given)"
fi

SORT_TMP="${TMP_DIR}/$(basename "${OUT%.bam}")_sort"

# STAR refuses to start if its scratch directory already exists, which is
# exactly the state a killed or failed run leaves behind:
#
#   EXITING because of fatal ERROR: could not make temporary directory
#
# Snakemake clears the declared outputs before re-running a rule but knows
# nothing about this directory, so clear it here.
rm -rf "${PREFIX}_STARtmp"

STAR \
    --genomeDir "$INDEX" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --runThreadN "$THREADS" \
    --outSAMtype BAM Unsorted \
    --outStd BAM_Unsorted \
    --outSAMattributes NH HI AS nM \
    --outFileNamePrefix "$PREFIX" \
    ${STAR_ARGS[@]+"${STAR_ARGS[@]}"} \
    | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEMORY" -T "$SORT_TMP" -o "$OUT"

echo "wrote $OUT"
