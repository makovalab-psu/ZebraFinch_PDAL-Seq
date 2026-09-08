#!/bin/bash
#
# js4025 Phase 2 -- merge the well-correlated Zebra finch PDAL-Seq
# datasets and compare read density to non-B DNA motif density.
#
# Submit from the project root with:
#
#     mkdir -p slurm_logs && sbatch submit_phase2.sh
#
# The mkdir is NOT optional and cannot move into this script: SLURM opens
# the --output file before the first line of the body runs, so a missing
# slurm_logs/ kills the job with no log explaining why.
#
#SBATCH --job-name=js4025_phase2
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/js4025_phase2_%j.out
#SBATCH --error=slurm_logs/js4025_phase2_%j.err
#
# Cores: 10, not Phase 1's 20. The long pole here is samtools merge, which
# is I/O bound rather than CPU bound -- more cores buy nothing and lengthen
# the queue wait. See the note on the reserved_resources reservation in
# submit_phase1.sh and README.md before assuming a pending job is stuck.
#
# Partition: basic (there is no `open` partition on this cluster; see
# submit_phase1.sh). basic is 4 GB/core, so --mem-per-cpu must stay at 4G.
# basic consumes allocation credit.
#
# Wall time: three real merges (~117 M, ~114 M and ~231 M reads) plus two
# single-library copies, then megadepth over 5 experiments x 3 window
# sizes, bedtools coverage over 8 motif classes x 3 window sizes, and the
# R compile. Expect 3-4 h against a 12 h request.
#
# Disk: the merged BAMs add roughly 65-70 GB across both MapQ branches, on
# top of Phase 1's ~29 GB. Check quota before submitting:
#     check_storage_quotas

set -euo pipefail

# Same driver environment as Phase 1 -- snakemake plus a conda new enough
# to build rule environments. Created once with:
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

# Snakemake resolves conda in a child shell. Three things have to be true for
# it to find the driver environment's conda rather than the module's 4.11.0,
# and job 55276541 proved that `conda activate` alone delivers none of them:
# it activated correctly (CONDA_PREFIX was right, the env's binary reported
# 26.7.1) and snakemake still refused with
#
#   CreateCondaEnvironmentException: Conda must be version 24.7.1 or later,
#   found version 4.11.0.
#
# `module load anaconda` installs `conda` as a shell FUNCTION that ignores PATH
# and calls $CONDA_EXE, and Lmod exports that function, so it survives into the
# non-interactive child shell -- where $CONDA_EXE still points at the module.
# So: fix PATH, repoint CONDA_EXE, and drop the exported function.
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
# through this shell. If the exported function were still in play this would
# report the module's version while the line above reports the environment's,
# which is exactly the failure being guarded against.
CONDA_VERSION="$(bash -c 'conda --version' 2>/dev/null | awk '{print $2}' || true)"
echo "conda      : ${CONDA_EXE} $("${CONDA_EXE}" --version | awk '{print $2}')"
echo "             child shell resolves conda ${CONDA_VERSION:-<none>}  (what snakemake sees)"

if [[ -z "${CONDA_VERSION}" ]]; then
    echo "ERROR: a child shell could not run conda at all." >&2
    exit 1
fi

# Fail here with an explanation rather than 30 seconds later inside snakemake
# with a stack trace. sort -V puts the older version first, so if the lowest of
# the two is not the requirement, what we have is older than the requirement.
CONDA_MIN="24.7.1"
if [[ "$(printf '%s\n%s\n' "${CONDA_VERSION}" "${CONDA_MIN}" | sort -V | head -1)" != "${CONDA_MIN}" \
      && "${CONDA_VERSION}" != "${CONDA_MIN}" ]]; then
    echo "ERROR: snakemake needs conda >= ${CONDA_MIN} to build rule environments, but a" >&2
    echo "       child shell resolves conda ${CONDA_VERSION} -- Roar's module conda, not the" >&2
    echo "       driver environment's ($("${CONDA_EXE}" --version | awk '{print $2}'))." >&2
    echo "       If the driver env's own conda is also old, refresh it (this does not" >&2
    echo "       touch the module):" >&2
    echo "         module load anaconda" >&2
    echo "         conda install -n js4025_snakemake -c conda-forge 'conda>=${CONDA_MIN}' -y" >&2
    echo "       See README section 02, 'If snakemake refuses to build the rule environments'." >&2
    exit 1
fi

# The rank sorts in clean_nonb.sh / combine_nonb.sh and the blacklist sort
# all write here.
mkdir -p temp

# --resources mem_mb is a GLOBAL cap on concurrently running jobs, not a
# per-job limit, and it is DERIVED from the allocation rather than
# hardcoded. A cap above what SLURM granted makes snakemake oversubscribe
# and the kernel OOM-kills mid-pipe; the visible error is usually a
# downstream "truncated file", with the real "Killed" line further up.
MEM_PER_CPU_MB=4096
TOTAL_MEM_MB=$(( SLURM_NTASKS * MEM_PER_CPU_MB ))
SNAKE_MEM_MB=$(( TOTAL_MEM_MB * 95 / 100 ))

echo "memory     : ${SLURM_NTASKS} cores x ${MEM_PER_CPU_MB} MB = ${TOTAL_MEM_MB} MB, capping snakemake at ${SNAKE_MEM_MB} MB"

# The largest single rule request in Phase 2 is compile_reads_vs_nonb at
# mem_mb=16000, so the allocation must clear that or it can never run.
if [[ "${SNAKE_MEM_MB}" -lt 16000 ]]; then
    echo "ERROR: ${SNAKE_MEM_MB} MB is below the 16000 MB that compile_reads_vs_nonb requests." >&2
    echo "       Increase --ntasks (at 4 GB/core, 5 cores is the floor)." >&2
    exit 1
fi

snakemake \
    --snakefile workflow/Snakefile \
    --cores "${SLURM_NTASKS}" \
    --use-conda \
    --resources mem_mb="${SNAKE_MEM_MB}" \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going

echo "=================================================================="
echo "finished    : $(date)"
echo "=================================================================="

sacct -j "${SLURM_JOB_ID}" --format=JobID,JobName,MaxRSS,Elapsed,TotalCPU,State
