#!/bin/bash
#
# js4025 Phase 1 -- preprocess Zebra finch PDAL-Seq / WGS reads and
# correlate them in genome windows.
#
# Submit from the project root with:
#
#     mkdir -p slurm_logs && sbatch submit_phase1.sh
#
# The mkdir is NOT optional and cannot move into this script: SLURM opens
# the --output file before the first line of the body runs, so a missing
# slurm_logs/ kills the job with no log explaining why.
#
#SBATCH --job-name=js4025_phase1
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_logs/js4025_phase1_%j.out
#SBATCH --error=slurm_logs/js4025_phase1_%j.err
#
# Partition: this was --partition=open until 2026-09-01, when SLURM
# answered with
#     slurm_job_submit: Open partition does not exist. Using basic
#     partition as default.
# and silently rewrote it. It is set to basic explicitly so the header
# says what actually runs. NOTE: unlike the open queue, basic consumes
# allocation credit. To use a paid allocation deliberately:
#     #SBATCH --partition=sla-prio
#     #SBATCH --account=<your_allocation_id>
# basic takes no --account line.
#
# Memory: basic is 4 GB/core, so 20 cores x 4 GB = 80 GB, and
# --mem-per-cpu must not exceed that rate. Estimated peak is ~30 GB --
# two concurrent map_reads jobs (bwa mem on a 1.14 Gbp index is ~6 GB,
# plus samtools sort at 4 threads x 2 GB). If
# benchmarks/phase1/map_reads/*.tsv shows max_rss above 4 GB/core, ask
# for more cores rather than more memory per core.
#
# Wall time: bwa index is ~1-1.5 h and blocks everything; then 6 samples
# at threads:10 run two at a time, ~2-3 h each, so 3 waves; then
# megadepth, statistics and bigwigs. Expect 12-16 h of work.
#
# On a pending job, squeue reason "(Reservation)" on this cluster is
# normally NOT a maintenance window. There is a standing ACTIVE MAGNETIC
# reservation ("reserved_resources", ~156 nodes, preempt/reserved QOS,
# runs to 2030) that permanently removes those nodes from the pool, so
# ordinary jobs wait for a node outside it. Check the predicted start
# with `squeue -j <id> --start` before changing anything.
#
# A shorter --time improves the odds of landing in a backfill gap, which
# is the main lever on queue wait here. 24:00:00 still leaves ~8 h of
# headroom over the 12-16 h estimate, and overrunning the limit is
# recoverable: every completed rule's output is on disk, so resubmitting
# resumes rather than restarting (--rerun-incomplete is already passed).

set -euo pipefail

# The driver environment (snakemake + a modern conda/mamba), NOT the
# per-rule envs. Those are built by snakemake --use-conda from
# workflow/env/*.yml. Create it once with:
#     module load anaconda
#     conda env create -f workflow/env/snakemake.yml
#
# The env pins conda>=24.7.1 because Roar's anaconda module ships conda
# 4.11.0 and snakemake 8 refuses to build rule environments with it. See
# the README section "If snakemake refuses to build the rule environments".
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

# samtools sort, the uniq -c sorts and the bedGraph sort all write here.
mkdir -p temp

# --resources mem_mb is a GLOBAL cap on concurrently running jobs, not a
# per-job limit. It is DERIVED from the allocation, never hardcoded: if it
# exceeds what SLURM actually granted, snakemake oversubscribes and the
# kernel OOM-kills mid-pipe (the visible error is usually a downstream
# samtools "truncated file", with the real "Killed" line further up). If
# it is set too low, snakemake serialises the run while cores sit idle.
#
# 4096 MB/core is the basic partition rate and must match --mem-per-cpu
# above. 5% is held back for the snakemake process itself.
MEM_PER_CPU_MB=4096
TOTAL_MEM_MB=$(( SLURM_NTASKS * MEM_PER_CPU_MB ))
SNAKE_MEM_MB=$(( TOTAL_MEM_MB * 95 / 100 ))

echo "memory     : ${SLURM_NTASKS} cores x ${MEM_PER_CPU_MB} MB = ${TOTAL_MEM_MB} MB, capping snakemake at ${SNAKE_MEM_MB} MB"

# The largest single rule request is map_reads at mem_mb=16000, so the
# allocation must clear that or nothing can be scheduled at all.
if [[ "${SNAKE_MEM_MB}" -lt 16000 ]]; then
    echo "ERROR: ${SNAKE_MEM_MB} MB is below the 16000 MB that map_reads requests." >&2
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
