#!/bin/bash
#
# js4025 Phase 5 -- PDAL-Seq versus control read density in satellites,
# centromeres and transposable elements.
#
# Submit from the project root with:
#
#     mkdir -p slurm_logs && sbatch submit_phase5.sh
#
# The mkdir is NOT optional and cannot move into this script: SLURM opens the
# --output file before the first line of the body runs, so a missing
# slurm_logs/ kills the job with no log explaining why.
#
#SBATCH --job-name=js4025_phase5
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/js4025_phase5_%j.out
#SBATCH --error=slurm_logs/js4025_phase5_%j.err
#
# Cores: 20, as requested and as in Phases 3 and 4. This phase is 498
# independent single-threaded enrichment jobs at 2 GB each, so it is purely a
# question of how many slots fit; 20 cores on a 4 GB/core partition fits all
# 20 at once.
#
# Wall time: 24 h. The whole phase is
#
#   3 split jobs                     seconds
#   83 window-intersect jobs         ~1 min each, 20-way parallel
#   498 enrichment jobs              200 megadepth range queries each
#   compile and 2 plots              minutes
#
# The enrichment jobs are the only real cost and each is 100 iterations x
# (1 subsample + 1 shuffle) bigWig queries. At a second or two per query that
# is 3-7 min per job, so ~2-3 h at 20-way parallelism; 24 h leaves room for
# megadepth being slower than that on a 1.14 Gbp bigWig without risking the
# long queue wait a bigger request would buy.
#
# Nothing here re-reads a BAM. Phase 1 already wrote data/phase1/bigwig/, and
# that is the coverage source -- see the header of
# workflow/scripts/phase5_repeat_enrichment/Read_enrichment.sh for why a
# bigWig and not the BAM the manuscript used.
#
# Partition: basic (no --account line, exactly as in submit_phase1-4.sh).
# basic is 4 GB/core, so --mem-per-cpu must stay at 4G.
#
# Disk: negligible next to Phase 4. The per-class BEDs and their window
# clippings are ~1-2 GB in total, and the 498 enrichment CSVs are 101 lines
# each.
#
# Restarting is safe. Read_enrichment.sh stages its output and renames it
# atomically, so a job killed at the wall clock leaves no half-written CSV
# that a later run would read as a finished result. Resubmit until
# `snakemake -n` reports nothing to do.

set -euo pipefail

# Same driver environment as Phases 1-4. Created once with:
#     module load anaconda
#     conda env create -f workflow/env/snakemake.yml
CONDA_ENV="/storage/work/jus841/.conda/envs/js4025_snakemake"

echo "=================================================================="
echo "job id      : ${SLURM_JOB_ID}"
echo "job name    : ${SLURM_JOB_NAME}"
echo "partition   : ${SLURM_JOB_PARTITION}"
echo "nodelist    : ${SLURM_NODELIST}"
echo "ntasks      : ${SLURM_NTASKS}"
echo "submit dir  : ${SLURM_SUBMIT_DIR}"
echo "started     : $(date)"
echo "=================================================================="

cd "${SLURM_SUBMIT_DIR}"

module load anaconda
conda activate "${CONDA_ENV}"

# Snakemake resolves conda in a child shell, and `module load anaconda`
# installs `conda` as an exported shell FUNCTION that ignores PATH and calls
# $CONDA_EXE. `conda activate` alone therefore does not stop snakemake from
# finding the module's conda 4.11.0 and refusing to build rule environments.
# Fix PATH, repoint CONDA_EXE, drop the function. See submit_phase4.sh.
export PATH="${CONDA_ENV}/bin:${PATH}"
export CONDA_EXE="${CONDA_ENV}/bin/conda"
unset -f conda 2>/dev/null || true

if [[ ! -x "${CONDA_EXE}" ]]; then
    echo "ERROR: no conda binary at ${CONDA_EXE}." >&2
    echo "       Create the driver environment first:" >&2
    echo "         module load anaconda && conda env create -f workflow/env/snakemake.yml" >&2
    exit 1
fi

# Probe the way snakemake does -- through a child non-interactive bash, not
# through this shell.
CONDA_VERSION="$(bash -c 'conda --version' 2>/dev/null | awk '{print $2}' || true)"
echo "conda      : ${CONDA_EXE} $("${CONDA_EXE}" --version | awk '{print $2}')"
echo "             child shell resolves conda ${CONDA_VERSION:-<none>}  (what snakemake sees)"

if [[ -z "${CONDA_VERSION}" ]]; then
    echo "ERROR: a child shell could not run conda at all." >&2
    exit 1
fi

CONDA_MIN="24.7.1"
if [[ "$(printf '%s\n%s\n' "${CONDA_VERSION}" "${CONDA_MIN}" | sort -V | head -1)" != "${CONDA_MIN}" \
      && "${CONDA_VERSION}" != "${CONDA_MIN}" ]]; then
    echo "ERROR: snakemake needs conda >= ${CONDA_MIN} to build rule environments, but a" >&2
    echo "       child shell resolves conda ${CONDA_VERSION} -- Roar's module conda, not the" >&2
    echo "       driver environment's ($("${CONDA_EXE}" --version | awk '{print $2}'))." >&2
    echo "       See README section 02, 'If snakemake refuses to build the rule environments'." >&2
    exit 1
fi

mkdir -p temp

# Each enrichment job makes a mktemp -d under temp/ and removes it on exit, but
# a job killed by the wall clock cannot run its trap. Clear the strays so they
# do not accumulate across resubmissions.
find temp -maxdepth 1 -type d -name 'read_enrichment_*' \
    -mmin +60 -exec rm -rf {} + 2>/dev/null || true

# --resources mem_mb is a GLOBAL cap on concurrently running jobs, not a
# per-job limit, and it is DERIVED from the allocation rather than hardcoded.
# A cap above what SLURM granted makes snakemake oversubscribe and the kernel
# OOM-kills mid-run.
MEM_PER_CPU_MB=4096
TOTAL_MEM_MB=$(( SLURM_NTASKS * MEM_PER_CPU_MB ))
SNAKE_MEM_MB=$(( TOTAL_MEM_MB * 95 / 100 ))

echo "memory     : ${SLURM_NTASKS} cores x ${MEM_PER_CPU_MB} MB = ${TOTAL_MEM_MB} MB, capping snakemake at ${SNAKE_MEM_MB} MB"

# The largest single-rule request in Phase 5 is 8000 MB (the split rules and
# the compile). Fail now rather than after the queue wait.
PHASE5_PEAK_MB=8000
if [[ "${SNAKE_MEM_MB}" -lt "${PHASE5_PEAK_MB}" ]]; then
    echo "ERROR: ${SNAKE_MEM_MB} MB is below the ${PHASE5_PEAK_MB} MB that the split and" >&2
    echo "       compile rules request. Increase --ntasks (at 4 GB/core, 3 is the floor)." >&2
    exit 1
fi

# Phase 5 targets by name rather than `rule all`, for the same reason
# submit_phase2_including_repetative.sh and submit_phase4.sh do it: `rule all`
# spans every phase, so anything snakemake decides is stale in Phases 1-4 -- a
# re-hashed environment, a touched resource file -- would restart mapping, a
# 4 h HMM fit or 15 h of STAR inside this job. Naming the targets makes that
# impossible.
#
# The two plots pull in both CSVs, which pull in all 498 enrichment jobs, so
# this list is the whole phase.
TARGETS=(
    "results/phase5/Repeat_annotation_report.txt"
    "plots/phase5/Repeat_enrichment_heatmap.svg"
    "plots/phase5/PDALSeq_versus_control.svg"
)

echo "targets    : ${#TARGETS[@]}"

snakemake \
    --snakefile workflow/Snakefile \
    --cores "${SLURM_NTASKS}" \
    --use-conda \
    --resources mem_mb="${SNAKE_MEM_MB}" \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going \
    "${TARGETS[@]}"

echo "=================================================================="
echo "finished    : $(date)"
echo "=================================================================="

sacct -j "${SLURM_JOB_ID}" --format=JobID,JobName,MaxRSS,Elapsed,TotalCPU,State
