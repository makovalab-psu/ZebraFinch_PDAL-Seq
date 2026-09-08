#!/bin/bash
#
# js4025 Phase 4 -- intersect the Zebra finch PDAL-Seq HMM states with non-B
# DNA motifs, 5mC rates, RNA-Seq and functional genomic annotations, for every
# k in PHASE4_STATES (8 through 14).
#
# Submit from the project root with:
#
#     mkdir -p slurm_logs && sbatch submit_phase4.sh
#
# The mkdir is NOT optional and cannot move into this script: SLURM opens the
# --output file before the first line of the body runs, so a missing
# slurm_logs/ kills the job with no log explaining why.
#
#SBATCH --job-name=js4025_phase4
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/js4025_phase4_%j.out
#SBATCH --error=slurm_logs/js4025_phase4_%j.err
#
# Cores: 20, as in Phase 3. This phase has two very different shapes of work.
# The RNA-Seq mapping is three long single-job stages that want every core;
# the enrichment sweep is ~1,500 independent two-minute jobs that want as many
# slots as memory allows. 20 cores serves both -- STAR takes the whole node
# while it runs, and afterwards the enrichment jobs at 2 GB each fill it.
#
# Wall time: 48 h, and the whole header is deliberately identical to
# submit_phase3.sh -- same partition, same nodes/ntasks/mem-per-cpu, same 48 h,
# and no --account line. Phase 3 submitted and ran with exactly this; Phase 4
# at 96 h did not, so the header is not the place to be inventive.
#
# Phase 4 is ~35-40 h of work, so 48 h should cover it in one go. If it does
# hit the wall clock, just resubmit: every stage is a separate snakemake
# target, the STAR BAMs are kept rather than temp(), and the enrichment scripts
# stage their output and rename it atomically, so a killed job loses only the
# handful of jobs that were actually in flight. Re-run until
# `snakemake -n` reports nothing to do.
#
# Budget:
#   STAR index                    ~0.5 h
#   fastp, 2 x ~460 M pairs       ~4 h
#   STAR + sort, 2 x ~460 M pairs ~15-20 h   <- the long pole
#   megadepth                     ~2 h
#   ~1,500 enrichment jobs        ~6-10 h at 20-way parallelism
#   compile and plot              minutes
# That is ~35-40 h of work. If you would rather cut it down than risk the wall
# clock, set RNA_SEQ_READ_LIMIT in workflow/Snakefile -- see README section 06.
#
# Partition: basic (no --account line, exactly as in submit_phase1/2/3.sh).
# basic is 4 GB/core, so --mem-per-cpu must stay at 4G, and 20 cores is what
# makes star_index's 40 GB request satisfiable at all.
#
# Disk: this is the heaviest phase for scratch. Peak usage is roughly
#   STAR index                       ~11 GB   (kept, resources/genomes/)
#   trimmed fastqs                  ~160 GB   (temp)
#   sorted RNA-Seq BAMs             ~140 GB   (kept, data/phase4/RNA_seq/sorted)
#   enrichment CSVs, states, logs     ~2 GB
# There is no duplicate-removal step, so the sorted BAM is the analysis BAM and
# is not temp(). Snakemake deletes the trimmed fastqs as soon as STAR finishes
# with them, but the two RNA-Seq libraries are processed independently and can
# overlap, so plan for ~320 GB free before submitting.

set -euo pipefail

# Same driver environment as Phases 1-3. Created once with:
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

mkdir -p temp

# The enrichment jobs each make a mktemp -d under temp/ and remove it on exit,
# but a job killed by the wall clock cannot run its trap. Clear the strays now
# so they do not accumulate across resubmissions.
find temp -maxdepth 1 -type d \( -name 'enrichment_*' -o -name 'bg_enrichment_*' \) \
    -mmin +60 -exec rm -rf {} + 2>/dev/null || true

# --resources mem_mb is a GLOBAL cap on concurrently running jobs, not a
# per-job limit, and it is DERIVED from the allocation rather than
# hardcoded. A cap above what SLURM granted makes snakemake oversubscribe
# and the kernel OOM-kills mid-run; the visible error is usually something
# misleading further downstream, with the real "Killed" line further up.
MEM_PER_CPU_MB=4096
TOTAL_MEM_MB=$(( SLURM_NTASKS * MEM_PER_CPU_MB ))
SNAKE_MEM_MB=$(( TOTAL_MEM_MB * 95 / 100 ))

echo "memory     : ${SLURM_NTASKS} cores x ${MEM_PER_CPU_MB} MB = ${TOTAL_MEM_MB} MB, capping snakemake at ${SNAKE_MEM_MB} MB"

# The largest single-rule request in Phase 4 is star_index at 40000 MB. If the
# cap is below that the job can never run, so fail now rather than after the
# queue wait.
PHASE4_PEAK_MB=40000
if [[ "${SNAKE_MEM_MB}" -lt "${PHASE4_PEAK_MB}" ]]; then
    echo "ERROR: ${SNAKE_MEM_MB} MB is below the ${PHASE4_PEAK_MB} MB that star_index requests." >&2
    echo "       Increase --ntasks (at 4 GB/core, 11 cores is the floor)." >&2
    exit 1
fi

# Phase 4 targets by name rather than `rule all`, for the same reason
# submit_phase2_including_repetative.sh does it: `rule all` spans every phase,
# so anything snakemake decides is stale in Phases 1-3 -- a re-hashed
# environment, a touched resource file -- would restart mapping or a 4 h HMM
# fit inside this job. Naming the targets makes that impossible.
#
# The two plots per k pull in all six CSVs for that k as dependencies, so this
# list is the whole phase. Keep PHASE4_STATES in step with workflow/Snakefile.
PHASE4_STATES=(8 9 10 11 12 13 14)

TARGETS=(
    "results/phase4/Annotation_report.txt"
    "results/phase4/Segmented_genome.bed"
    "results/phase4/RNA_seq_mapping_statistics.txt"
)
for K in "${PHASE4_STATES[@]}"; do
    TARGETS+=(
        "plots/phase4/${K}/HMM_functional_summary.svg"
        "plots/phase4/${K}/State_distribution_by_chromosome.svg"
    )
done

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
