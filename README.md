# Overview

Code for JPSieg's analysis of Zebra finch cell line PDAL-Seq data (notebook page number).

Contains final plots and a Snakemake pipeline in "workflow".

Ran in five phases:

Phase 1: Data preprocessing and correlation

Phase 2: Merge replicates, construct blacklist, and examine correspondence to non-B DNA motifs

Phase 3: Choose optimum HMM states

Phase 4: PDAL-Seq HMM correspondence to genome function

Phase 5: PDAL-Seq read enrichment in repeats

The "workflow/" directory contains conda environments in "env", rules for each step in "rules", scripts for each step in "scripts", and the overall workflow configuration in "Snakefile".

Requires large "resources" and "data/raw_reads/" directories. Contact authors for access. Will create large intermediate "data" and "results" folders, along with smaller "logs" and "benchmarks". Requires 648G of total space following implementation. Recommend 1 TB disk total for intermediate files.

To run:

1. Place "resources" and "data/raw_reads/" in a directory.
2. Install the driver environment

```bash
module load anaconda
mamba env create -f workflow/env/snakemake.yml
mamba activate js4025_snakemake
```
4. Test the dag:

```bash
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda
```

5. Implement the workflow directly

```bash
snakemake \
    --snakefile workflow/Snakefile \
    --cores "${CORES}" \
    --use-conda \
    --resources mem_mb="${MEM_MB}" \
    --rerun-incomplete \
    --printshellcmds \
    --keep-goingh
```

Or modify the submit slurm scripts to submit batch jobs on your sever. Example:

```bash
mkdir -p slurm_logs && sbatch submit_phase1.sh
```

See my step by step configuration and implementation below:

# 01 Upload and check data

#### Download data from the Huck portal to data/raw_reads. Check the data:

#Analysis_1
ls -l data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_1.zip 
-rw-r--r-- 1 jus841 kdm16_collab 83432425233 Aug 31 13:35 data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_1.zip
md5sum data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_1.zip
6858c72aaf82e6ab35e1ae70ed9640ac  data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_1.zip

#Analysis_2
ls -l data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_2.zip 
-rw-r--r-- 1 jus841 kdm16_collab 78191207072 Aug 31 15:50 data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_2.zip
md5sum data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_2.zip 
3a00b583d24fc03e7c1e677ddc527a2c  data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_2.zip


cd data/raw_reads/
unzip data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_1.zip
unzip data/raw_reads/260826_VH00707_339_AAJJGHNM5_Analysis_2.zip

# 02 Phase 1 -- Preprocess reads and compare datasets in genome windows

Trim, map, deduplicate and quantify every R1 library in 10 kbp (plus 100 kbp and
1 Mbp) windows across bTaeGut7v0.4, then correlate every dataset against every
other one. The output answers the one question Phase 2 needs settled first:
**which libraries are technically sound, and which of them are similar enough to
merge as replicates.**

Adapted from steps 01, 02, 03 and 21 of the PDAL-Seq manuscript pipeline
(`Manuscripts/PDALseq_manuscript/PNAS_revised_version_2/PDALSeq_manuscript`).

## Inputs this phase expects

| Path | Where it came from |
|---|---|
| `Zebrafinch_samples.txt` | in this repo -- sample sheet, 6 R1 libraries |
| `data/raw_reads/fastq/Tcas_CFS414_*_R1_001.fastq.gz` | section 01 above. R2 files are ignored |
| `resources/genomes/bTaeGut7v0.4_mat_Z_MT_rDNA.fa` + `.fa.fai` | `resources/README.md`, `make_PDAL-Seq_fasta.sh` |

BWA indices are **not** expected -- `rule bwa_index` builds them into
`resources/genomes/`.

## Install the environments (once)

```bash
module load anaconda

# Driver environment: snakemake itself. This is the only env you create by
# hand; snakemake builds the per-rule envs from workflow/env/*.yml on first run.
conda env create -f workflow/env/snakemake.yml
```

### If snakemake refuses to build the rule environments

```
CreateCondaEnvironmentException:
Conda must be version 24.7.1 or later, found version 4.11.0.
```

Roar's `module load anaconda` provides conda 4.11.0. Snakemake 8 will not build
rule environments with it. `workflow/env/snakemake.yml` now pins
`conda>=24.7.1` and `mamba` inside the driver environment, which is the fix
snakemake itself recommends: once `js4025_snakemake` is activated its `bin/`
comes first on `PATH`, so snakemake finds the new conda and the module's old
one is never modified.

If the environment was created before that pin was added, add the new conda to
it in place -- no need to rebuild:

```bash
module load anaconda
conda install -n js4025_snakemake -c conda-forge "conda>=24.7.1" mamba -y

# "conda is not recommended in non-base environments" is expected here; that
# is exactly what is being done, deliberately.

conda activate js4025_snakemake

# Check the BINARY, not `conda --version`. `conda` is a shell function
# installed by `module load anaconda`; it ignores PATH and calls $CONDA_EXE,
# which still points at the module's 4.11.0. So `conda --version` reports
# 4.11.0 forever, even with 26.x correctly installed in the env.
"$CONDA_PREFIX/bin/conda" --version     # this is the one that matters
```

Snakemake shells out to `conda` from a non-interactive bash where that
function does not exist, so it resolves through `PATH`. The real test is simply
that `snakemake -n --use-conda` stops raising `CreateCondaEnvironmentException`.

**Inside a batch job, `conda activate` is not enough, and `PATH` is not the
mechanism.** Job 55276541 activated the driver environment correctly --
`CONDA_PREFIX` was right and `$CONDA_PREFIX/bin/conda` reported 26.7.1 -- and
snakemake still refused, having found 4.11.0.

The reason is the shell function. `module load anaconda` defines `conda` as a
function that ignores `PATH` and calls `$CONDA_EXE`, and Lmod **exports** it, so
it survives into the non-interactive child shell snakemake spawns -- where
`$CONDA_EXE` still points at the module. Reproduced in isolation:

| Probe | Reports |
|---|---|
| `$CONDA_PREFIX/bin/conda --version` | 26.7.1 |
| `$(type -P conda) --version` | 26.7.1 |
| `bash -c 'conda --version'` | **4.11.0** <- what snakemake sees |

So a `PATH` fix alone does not work, and neither does checking `type -P` --
both report the right answer while snakemake gets the wrong one. All four
`submit_*.sh` scripts do this after activation:

```bash
export PATH="${CONDA_ENV}/bin:${PATH}"
export CONDA_EXE="${CONDA_ENV}/bin/conda"
unset -f conda 2>/dev/null || true
```

and then verify through a child shell -- `bash -c 'conda --version'`, the same
resolution snakemake performs -- failing immediately with an explanation if that
is older than 24.7.1, instead of a `CreateCondaEnvironmentException` stack trace
thirty seconds later.

`--conda-frontend mamba` was also dropped from all four scripts: the installed
snakemake now reports `Support for alternative conda frontends has been
deprecated ... Ignoring the alternative conda frontend setting (mamba)`.

Before any of that, it is worth seeing whether Roar has a newer module, which
would make this unnecessary:

```bash
module avail anaconda
```

## Run it

```bash
conda activate js4025_snakemake

# Check the DAG first. Should report 135 jobs: bwa_index 1, trim_fastq 6,
# map_reads 6, genome_window_coverage 36, index_bam 12, bam_flags 18.
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda

# Submit. The mkdir is required: SLURM opens the --output file before the
# job body runs, so a missing slurm_logs/ kills the job with no log.
mkdir -p slurm_logs && sbatch submit_phase1.sh
```

`submit_phase1.sh` requests 1 node, 20 cores, 4 GB/core (80 GB) on
`--partition=basic` for 48 h.

**There is no `open` partition on this cluster.** The script asked for it on
2026-09-01 and SLURM answered `slurm_job_submit: Open partition does not exist.
Using basic partition as default.`, silently rewriting the request. The header now
says `basic` so it matches what runs. Unlike an open queue, **basic consumes
allocation credit**; to use a paid allocation deliberately, swap the partition line
for `--partition=sla-prio` plus `--account=<allocation_id>` (basic itself takes no
account line). `basic` is 4 GB/core, so `--mem-per-cpu` must not exceed `4G` --
to get more memory, ask for more cores.

If `squeue` shows the job pending with reason `(Reservation)`, check the predicted
start before assuming anything:

```bash
squeue -j <jobid> --start
scontrol show reservation
```

On this cluster that reason is usually **not** a maintenance window. A standing
`ACTIVE MAGNETIC` reservation (`reserved_resources`, ~156 nodes, `preempt,reserved`
QOS, running to 2030) permanently removes those nodes from the general pool, so an
ordinary job simply waits for a node outside it. Observed 2026-09-01: a 20-core,
48 h job queued at ~21:00 was scheduled to start the next day at 11:56.

A shorter `--time` is the main lever -- it improves the odds of fitting a backfill
gap. 24 h leaves ~8 h of headroom over the 12-16 h estimate, and overrunning is
recoverable: completed rule outputs stay on disk, so resubmitting resumes rather
than restarting.

### Do not trade cores for a shorter queue wait

`map_reads` is the bottleneck and it is thread-bound (~481 M reads total). At 20
cores two `bwa mem` jobs run at 10 threads each, ~7 h; at 8 cores only one fits
(snakemake caps `threads: 10` down to the available cores) and mapping stretches
to ~17 h, pushing the run to ~22-26 h. Dropping cores buys queue position and
pays for it twice over in runtime.

Ask SLURM instead of guessing -- `--test-only` predicts a start without queuing
anything:

```bash
sbatch --test-only --ntasks=20 --time=24:00:00 submit_phase1.sh
sbatch --test-only --ntasks=8  --time=24:00:00 submit_phase1.sh
```

If you do change `--ntasks`, **change nothing else**: the script derives
snakemake's global `--resources mem_mb` from `$SLURM_NTASKS` x 4096 MB and aborts
if the result falls below the 16000 MB that `map_reads` requests (5 cores is the
floor). A hardcoded cap larger than the allocation makes snakemake oversubscribe
and the kernel OOM-kills mid-pipe -- which surfaces as a misleading samtools
"truncated file" error, with the real `Killed` line further up the log.

Restarting after a kill needs `--rerun-incomplete`, which the submit script
already passes; run `snakemake --unlock` first if snakemake died holding the lock.

### If only the heatmap failed

Observed on the first run (job 55234228, 2026-09-01): 133 of 135 steps finished
and `plot_correlation` died with

```
Error in dev(...): The package "svglite" is required to save as SVG.
```

`ggsave()` dispatches on the file extension, and `.svg` needs `svglite`, which
`Rutility.yml` does not carry. The fix is a **separate** `workflow/env/Rplot.yml`
holding `r-svglite`; `plot_correlation` now points at it. Adding svglite to
`Rutility.yml` instead would change that file's content hash, and snakemake would
rerun both completed `compile_*` jobs to match -- which is why every plotting rule
from here on should use `Rplot.yml`, and why neither file gets edited for a later
phase.

Only one job remains, so resubmit small rather than asking for the full node
again. 5 cores is the floor the script's memory guard allows:

```bash
sbatch --ntasks=5 --time=1:00:00 submit_phase1.sh
```

The first few minutes are mamba solving `Rplot.yml`; nothing else reruns.

## Outputs

| File | What it answers |
|---|---|
| `results/phase1/Preprocessing_statistics.csv` | one row per library: raw and trimmed reads, % mapped, duplicates removed, % MapQ >= 20, median MapQ, % uniquely placed |
| `results/phase1/fastp_reports/*.html` | per-library trimming QC |
| `results/phase1/Window_coverage_correlation.csv` | every sample pair x 3 window sizes x 2 MapQ branches, Pearson and Spearman, with `Comparison` and `Same_experiment` columns |
| `results/phase1/Window_coverage_10000_CPM.csv` | wide 10 kbp coverage matrix, normalized to coverage per million. **This is the Phase 2 input.** |
| `plots/phase1/Correlation_heatmap.svg` | the correlation table as a heatmap, faceted by window size and MapQ branch |
| `data/phase1/All_reads/*.bam`, `data/phase1/High_MapQ/*.bam` (+ `.bai`) | deduplicated alignments, with and without a MapQ >= 20 filter |
| `data/phase1/bigwig/*.bigwig` | per-library browser tracks |
| `benchmarks/phase1/<rule>/*.tsv` | runtime and `max_rss` per job -- the memory budget for Phase 2 |

## Sanity checks after the run

1. `Raw_reads` in `Preprocessing_statistics.csv` must match Table js4025.1
   (119, 60.5, 76.2, 75.9, 70.7, 78.4 million).
2. `Percent_duplicates_removed` should land near the MultiQC `Dups` column
   (0.109 for the control, 0.189-0.236 for PDAL-Seq).
3. Every self-pair in `Window_coverage_correlation.csv` must have
   `Spearman == 1`. That is the cheapest proof the pairing logic is right.
4. In the heatmap, 20mM-1/20mM-2 and 40mM-1/40mM-2 should be the brightest
   off-diagonal cells if those replicates are mergeable in Phase 2.

## Decisions recorded for this phase

- **R1 only.** No rule anywhere reads an R2 file.
- **The PCR-free WGS control is deduplicated** with `samtools rmdup -s`, exactly
  like the PDAL-Seq libraries. A PCR-free library has no PCR duplicates, so this
  removes some genuinely independent fragments, but it keeps the Phase 2
  denominator processed identically to the numerator, which is what the
  manuscript did with its TruSeq controls.
- **`bwa mem` streams into `samtools sort`** rather than writing a temp SAM as
  the manuscript does. A 119 M read SAM is ~40 GB on scratch for no benefit.
- **`make_windows.sh` keeps the manuscript's `$2+1` start shift.** It skips the
  first base of every window, but it means window coordinates here line up with
  every window coordinate in the manuscript. Do not "fix" it without
  regenerating all downstream coverage.
- **`Estimate_unique_maps.sh` fixes a bug in the manuscript version.** The
  original awk never reset `as_value`/`xs_value` between records, so a read with
  no `XS` tag -- i.e. a uniquely placed read, the exact case being counted --
  inherited the previous read's `XS`. Both are reset per record here, so the
  `Percent_uniquely_placed` column is not comparable to the manuscript's Table S3.
- **Coverage is summed per-base coverage, not read counts.** `megadepth --op sum`
  reports base-pairs of coverage per window; the normalized matrix is labelled
  CPM (coverage per million) for that reason.
- **Blacklist generation and per-experiment BAM merging are deliberately not in
  this phase.** Both depend on knowing which datasets correlate, so they belong
  to Phase 2.
- Chromosome ordering: `bedGraphToBigWig` needs ASCII-sorted chromosomes, and the
  `.fai` order is `chr1, chr1A, chr2, ...`, so `bg_to_bigwig.sh` sorts the
  bedGraph first. The sort is not optional.

# 03 Phase 2 -- Merge datasets and compare read density to non-B DNA motif density

Phase 1 showed the PDAL-Seq libraries correlate tightly with one another and
much less well with the PCR-free WGS control (Figure js4025.1), so the
replicates are mergeable. This phase merges them, builds the blacklist that
Phase 1 deferred, and asks the actual question: **does PDAL-Seq read density
track non-B DNA motif density in a way the control does not?**

Adapted from steps 03 (blacklist), 04 (non-B comparison) and 21 (the
per-experiment merge pattern) of the PDAL-Seq manuscript pipeline.

## Experiments

| Merged output | Libraries | Dedup reads |
|---|---|---|
| `Tguttata_CFS414_WGS` | PCRfree | ~106 M |
| `Tguttata_CFS414_PDALSeq_0mM` | 0mM | ~49 M |
| `Tguttata_CFS414_PDALSeq_20mM` | 20mM-1 + 20mM-2 | ~117 M |
| `Tguttata_CFS414_PDALSeq_40mM` | 40mM-1 + 40mM-2 | ~114 M |
| `Tguttata_CFS414_PDALSeq` | all four 20 mM + 40 mM libraries | ~231 M |

Grouping comes from the `Experiment` column of `Zebrafinch_samples.txt`; the
combined dataset is assembled in the Snakefile. 0 mM is deliberately excluded
from the combined set -- it is the untreated control, not a replicate. Both
MapQ branches are merged; Phase 3's HMM will want `High_MapQ`.

The two single-library experiments still go through `samtools merge`, where it
is effectively a copy. That is intentional: uniform output paths mean every
downstream rule has one input pattern. It costs ~12 GB.

## Inputs this phase expects

| Path | Where it came from |
|---|---|
| `data/phase1/{All_reads,High_MapQ}/*.bam` | Phase 1 |
| `data/phase1/genome_windows/{1000000,100000,10000}_nucleotides/*.bed` | Phase 1 |
| `data/phase1/genome_window_coverage/{All_reads,High_MapQ}/10000_nucleotides/*.bed` | Phase 1 -- the blacklist reads these, nothing is recomputed |
| `resources/nonB_DNA_annotations/final_nonB/bTaeGut7v0.4_MT_rDNA.{APR,DR,G4,IR,STR,TRI,Z}.bed` | Zenodo 10.5281/zenodo.19225047, see `resources/README.md` |

### Seven motif classes, not the manuscript's ten

The Zebra finch release annotates `APR, DR, G4, IR, STR, TRI, Z`. `CRU`, `MR`
and `SLS` were not annotated for this genome, and what the manuscript calls
`GQ` is called `G4` here. The pipeline uses the release's names.

### The annotations cover both haplotypes

They were built against the full `bTaeGut7v0.4_MT_rDNA` assembly, so every
autosome appears twice -- `chr10_mat` **and** `chr10_pat` -- while our genome is
the single haplotype `bTaeGut7v0.4_mat_Z_MT_rDNA`. `clean_nonb.sh` filters each
file to the 45 contigs in our `.fai`, so **roughly half of every class is
expected to drop**. Check that in `results/phase2/nonB_annotation_report.txt`;
the rule hard-fails, printing the observed chromosome names, if a class loses
everything.

Sorting is by `.fai` **rank**, not ASCII: the order is `chr1_mat, chr1A_mat,
chr2_mat, ...`, which is neither lexicographic nor numeric, and
`bedtools coverage -sorted -g` requires exactly it. Same class of trap as the
`bedGraphToBigWig` sort in Phase 1. Do not replace the rank sort with
`sort -k1,1`.

## Run it

Nothing new to install by hand -- snakemake builds `workflow/env/phase2_nonb.yml`
on first use. `Rutility.yml` and `Rplot.yml` are reused **unmodified**, so no
Phase 1 job reruns.

```bash
conda activate js4025_snakemake

# Should report 77 new jobs and nothing from Phase 1.
snakemake --snakefile workflow/Snakefile -n --cores 10 --use-conda

mkdir -p slurm_logs && sbatch submit_phase2.sh
```

`submit_phase2.sh` requests 1 node, 10 cores, 4 GB/core on `basic` for 12 h.
Ten rather than Phase 1's twenty: `samtools merge` is the long pole and it is
I/O bound, so extra cores buy nothing and lengthen the queue wait. As in
Phase 1, `--resources mem_mb` is derived from `$SLURM_NTASKS`, so change
`--ntasks` and nothing else.

**Disk**: the merged BAMs add roughly 65-70 GB on top of Phase 1's ~29 GB.
Run `check_storage_quotas` before submitting.

## Outputs

| File | What it answers |
|---|---|
| `results/phase2/nonB_correlation_coefficient.csv` | 864 rows: experiment x 8 classes x 3 window sizes x 2 motif metrics x 2 blacklist states, Pearson and Spearman |
| `results/phase2/nonB_coverage_versus_read.csv` | per-window table at 1 Mbp behind the scatter plot |
| `results/phase2/Blacklist.bed` | unmappable windows plus whole `chrW_mat` and `chrZ_pat` |
| `results/phase2/nonB_annotation_report.txt` | records in, records retained, % retained, contigs, per class |
| `results/phase2/nonb_annotations/*.bed` | cleaned, `.fai`-sorted motif annotations |
| `plots/phase2/nonB_correlation_heatmap.svg` | the correlation table as a diverging heatmap |
| `plots/phase2/nonB_scatter.svg` | PDAL-Seq vs WGS against motif density, faceted by class |
| `data/phase2/{All_reads,High_MapQ}/*.bam` (+ `.bai`) | merged experiment BAMs -- **the Phase 3 input** |

### Reading the correlation table

Four columns decide which number you want:

- `Signal` -- `CPM` (per-experiment density) or `Ratio_vs_WGS` (PDAL-Seq divided
  by the control, present only for PDAL-Seq experiments)
- `Metric` -- `count` (motifs per window, the manuscript's quantity) or
  `base_density` (fraction of the window covered). Six of seven classes
  self-overlap and the classes overlap each other, so `count` double counts,
  worst of all in `all`. `base_density` does not.
- `Blacklist_filtered` -- `TRUE` drops unmappable windows
- `Window` -- 1 Mbp, 100 kbp or 10 kbp

The heatmap plots `Signal == CPM`, `Blacklist_filtered == TRUE`.

## Sanity checks after the run

1. `nonB_annotation_report.txt` -- `Percent_retained` should be near 50 for
   every class (the `_pat` autosomes dropping out). A number near 100 means the
   annotation was single-haplotype after all; near 0 means the rule should have
   failed and something is wrong.
2. `Blacklist.bed` -- must contain full-length `chrW_mat` and `chrZ_pat` rows.
   The mappability portion ran at 7% of windows on synthetic data; a real value
   above ~30% suggests the threshold or the WGS control needs a look.
3. In `nonB_correlation_coefficient.csv`, the WGS control should sit near zero
   for most classes. It is the negative control for this whole comparison -- if
   WGS tracks motif density as strongly as PDAL-Seq does, the signal is
   compositional (GC, mappability), not permanganate reactivity.
4. `count` and `base_density` should broadly agree per class. Where they
   disagree sharply, trust `base_density`.

## Decisions recorded for this phase

- **Both merges, per the manuscript.** Concentration-level BAMs *and* a single
  combined PDAL-Seq BAM, matching `K562_merge_concentrations` in rule 21.
- **The blacklist reuses Phase 1 coverage.** The manuscript recomputes 10 kbp
  window coverage for this; Phase 1 already wrote exactly those files for both
  MapQ branches, so `make_raw_blacklist` reads them directly.
- **Sex chromosomes are blacklisted whole**, as the manuscript does for
  chrX/chrY. `chrMT` and the three `rDNA_morph_*` contigs are not blacklisted
  but are excluded from every correlation -- the rDNA morphs are alternative
  representations of one locus and so are guaranteed multi-mappers.
- **`all.bed` is concatenated, not merged**, reproducing the manuscript's
  `all.sh`. It therefore double counts overlapping motifs by design; that is
  what `base_density` is for.
- **Correlations are reported both ways on every axis** -- two metrics, two
  signals, two blacklist states -- because none of those choices has a single
  defensible default and all four are cheap once the coverage is in memory.
- **Experiment window coverage is computed on `All_reads` only**, matching the
  manuscript convention that downstream analyses use deduplicated BAMs with no
  MapQ filter. The merged `High_MapQ` BAMs exist for Phase 3.
- Heatmap label contrast is computed from the **actual WCAG luminance** of each
  cell's fill, not from a fraction of the ramp. Phase 1 could use "55% along the
  ramp" because viridis is monotonic in lightness; a diverging ramp is lightest
  in the middle, so that rule would put white text on pale cells.

# 04 Phase 3 -- Determine the optimum number of HMM states

Segment the genome with a Gaussian HMM in 1 kbp windows and answer one
question: **how many hidden states does the PDAL-Seq signal support?** Phase 3
stops at the diagnostics. Which *k* to use is a judgement made by eye from the
three figures below; interpreting the states functionally is Phase 4.

Adapted from `09_Human_HMM.smk` and the model scripts it borrows from
`07_PhyloHGMP_model` in the manuscript pipeline.

## The structural difference from the manuscript

The manuscript fit four human cell lines at once, so its Gaussian emission was
four dimensional. **We have one cell line, so ours is univariate.** Every rule
and script is written over `HMM_FEATURES`, a list, rather than a fixed column,
so adding the WGS control or the 0 mM library as a second dimension later is a
one-line Snakefile change and no rule rewrites.

## Configuration

Set at the top of `workflow/Snakefile`:

| Constant | Value | Note |
|---|---|---|
| `HMM_WINDOW` | `"1000"` | separate from `WINDOWS`; see below |
| `HMM_FEATURES` | `[COMBINED_EXPERIMENT]` | the merged 20 mM + 40 mM BAM |
| `HMM_STATES` | 2, 4, …, 20, 25, 30 | the manuscript stops at 20 |
| `HMM_TRANSFORM` | `"none"` | or `"log1p"` |

`HMM_WINDOW` is deliberately **not** appended to `WINDOWS`. Doing that would
re-run Phase 1 coverage for 6 samples x 2 MapQ branches and change the Phase 2
compile outputs. `make_genome_windows` and `experiment_window_coverage` both
take `{size}` as a free wildcard, so asking for 1000 here costs one extra
windows file and one extra megadepth pass and touches nothing else. A dry run
confirms it: `genome_window_coverage` stays at 36 jobs, not 42.

## Run it

```
rsync -Pav ./ submit.hpc.psu.edu:/storage/group/kdm16/default/jus841/js4-2026/js4025_Analyze_Zebra_finch_PDAL-Seq_data

module load anaconda
conda activate js4025_snakemake

# Should report 44 new jobs and nothing from Phase 1 or Phase 2.
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda

mkdir -p slurm_logs && sbatch submit_phase3.sh
```

`submit_phase3.sh` requests 20 cores and **48 h**. The twelve fits are
independent single-core jobs, so cores buy real parallelism here, unlike Phase
2 where the long pole was I/O-bound `samtools merge`. Measured on 10^6 windows
x 1 feature, `k = 30` costs ~14.5 s per EM iteration, so a fit that runs the
full `n_iter = 1000` without converging early is ~4 h; cost scales as `k^2`, so
`k = 30` sets the wall clock and everything smaller runs alongside it. The
fits are separate targets, so a job killed at the wall clock loses only what
had not finished and a re-submission resumes.

Disk is small compared with Phases 1-2: one extra windows file, one megadepth
pass, and twelve gzipped state assignments -- on the order of 1 GB.

## Outputs

| File | What it answers |
|---|---|
| `plots/phase3/HMM_BIC_scree.svg` | how many states BIC prefers, and whether that fit is trustworthy |
| `plots/phase3/Kmeans_scree.svg` | the model-free elbow, as a cross-check |
| `plots/phase3/gaussian_fit/<k>_states.svg` | is the model fitting the data? |
| `plots/phase3/HMM_model_vs_empirical_means.svg` | does each state's Gaussian match the windows it was given? |
| `results/phase3/HMM_BIC.csv` | BIC, log likelihood, parameter count, occupied states |
| `results/phase3/Kmeans_scree.csv` | wss and fraction of variance explained per *k* |
| `results/phase3/HMM_state_summary.csv` | every state of every model: fitted mean/variance, weight, empirical mean/sd |
| `data/phase3/GHMM/<k>/GHMM_states.csv.gz` | the segmentation itself, with `chr/start/end` -- Phase 4's input |

### Reading the BIC scree

Three panels, and they must be read together.

- **BIC** -- the manuscript's plot. The minimum is circled.
- **log likelihood (must not dip)** -- every model saw the same windows, so a
  model with more states cannot legitimately explain them worse. Where this
  curve dips, that *k*'s fit fell into a bad optimum and its BIC describes the
  fit rather than the data. Those points are marked with a red **x**, and
  `calc_BIC.py` prints a warning naming them.
- **States the model actually uses** -- a variational HMM is free to leave
  states empty. Where this curve leaves the grey diagonal, extra states are
  being fit but not used, which is independent evidence about the optimum.

**Bad fits at particular *k* are expected and are not a bug to tune away.**
`hmmlearn` initialises the emissions from k-means, and k-means splits wide
clusters and merges narrow ones, so a mixture with unequal variances will
mislead it at some *k* no matter how many times it is restarted (k-means
already runs `n_init=10` internally, and the seed is fixed at `random_state=42`
so results reproduce). On the synthetic three-state test data below, `k = 3`
merged two well-separated states and split a third, while `k = 4` recovered the
truth exactly. Read the three panels together and confirm the chosen *k*
against its `gaussian_fit` figure.

### Reading the Gaussian fit plots

The first panel is **all windows** with the weighted mixture
`sum_k w_k * dnorm(x; mu_k, sd_k)` drawn over the histogram; the rest are one
panel per occupied state, each against its own Gaussian. The manuscript plots
only the per-state panels, but a mixture can look fine state by state and still
miss the overall distribution, and the first panel is what catches that. Bars
and curves are on the same footing: per-state bars integrate to 1 over that
state, the "all windows" bars integrate to 1 over the genome.

If the Gaussians clearly miss the histograms, set `HMM_TRANSFORM = "log1p"` in
`workflow/Snakefile` and re-run. PDAL-Seq coverage in 1 kbp windows is zero
inflated with a long right tail, which a Gaussian emission may not describe
well. It is a `params` value and snakemake treats params as a rerun trigger, so
the whole phase re-runs from `normalize_hmm_coverage` and nothing in Phase 1 or
Phase 2 is recomputed.

## Sanity checks after the run

1. `wc -l data/phase3/filtered_coverage/*.bg` -- expect roughly 1.0 M windows
   (1.14 Gbp minus chrZ/chrW/chrMT/rDNA, minus the blacklist). A number near
   1.14 M means the blacklist intersect silently did nothing.
2. `results/phase3/HMM_BIC.csv` -- `occupied_states` should track `States` at
   small *k* and fall behind it once the model saturates. Check the log for the
   `WARNING: log likelihood dips at ...` line before trusting the BIC minimum.
3. `plots/phase3/Kmeans_scree.svg` -- the elbow should be in the same
   neighbourhood as the BIC minimum. Two independent methods disagreeing
   wildly is worth understanding before Phase 4.
4. `plots/phase3/HMM_model_vs_empirical_means.svg` -- points should sit on the
   `y = x` line. A point off the line means that model is not describing its
   own states, which is a reason to reject a *k* whatever BIC says.

## Verification done before this phase shipped

Verified locally on synthetic data with a planted three-state structure
(means 1000/3000/9000 summed coverage, sticky transitions, 1570 windows across
`chr1_mat`, `chr2_mat`, `chrMT` and `chrW_mat`):

- the blacklist + contig filter kept exactly 1450 of 1570 windows -- 1000 + 500
  autosomal, minus the 50 overlapping a planted unmappable interval -- matching
  an independent reimplementation row for row, with `chrMT` and `chrW_mat` gone
- normalisation returned mean 0, sd 1
- **`k = 4` recovered the planted states with purity 1.0000**: 562 / 559 / 329
  windows, every one assigned correctly, with the fourth state left empty
- BIC selected `k = 4` and the occupancy panel showed it using only 3 states
- `k = 3`, `5` and `6` were flagged by the log-likelihood check, and the `k = 3`
  fit plot shows the failure plainly: its Gaussian sits in the valley between
  two humps
- the k-means scree elbowed at 3
- `log1p` and the error paths (unknown transform, missing feature file) behave

Not verifiable locally, and stated as such: the conda solve for
`workflow/env/hmm.yml`, and behaviour at real genome scale.

## Decisions recorded for this phase

- **Univariate**, because there is one cell line. Written over a feature list
  so a second dimension is a one-line change.
- **1 kbp windows are a separate constant**, not an addition to `WINDOWS`, so
  Phase 1 and Phase 2 outputs are untouched.
- **`All_reads` only**, matching the manuscript's `window_coverage` rule. The
  blacklist is what handles multi-mappers.
- **Coordinates are carried through the whole chain.** The manuscript drops
  `chr/start/end` from the feature matrix and re-attaches them later by row
  position (`Annotate_state.R`). Carrying them removes that failure mode and
  hands Phase 4 ready-to-use state assignments.
- **Normalisation divides by `sd()`, not the mean.** The manuscript's
  `Normalize_human_window_coverage.R` sets `SD = mean(...)`. With a single
  feature this is only a scale factor -- it shifts every log likelihood by the
  same `N*log(scale)`, leaving the argmin of BIC unchanged -- but the plots are
  read on a z-scale, so it is fixed here.
- **`random_state` is fixed at 42**, so a re-run reproduces the same
  segmentation. The manuscript leaves it unset.
- **The parameter count in `calc_BIC.py` is the manuscript's exactly** (means +
  full covariances + transitions, omitting the initial-state probabilities), so
  the BIC values stay comparable with the human analysis. The `params` column
  is written out so BIC can be recomputed under another convention without
  refitting.
- **Signal is pre-binned in Python, not re-read in R.** `extract_GHMM_parameters.py`
  already holds the model and the data, so it writes the histogram counts
  directly; `Plot_gaussian_fit.R` then reads a few hundred rows instead of 10^6
  windows twelve times.
- **No new R environment.** Every Phase 3 figure is built from `tidyverse` +
  `svglite` + `viridis`, which `Rplot.yml` already has. The manuscript's
  `plot_to_test_Guassian.R` uses nested `gridExtra::grid.arrange`; that is
  rewritten as ggplot facets rather than adding `gridExtra` to `Rplot.yml`,
  which would rehash the env and re-run every finished Phase 1 and Phase 2
  figure.
- **`fit_ghmm` memory is measured, not guessed**: 1.19 GB peak RSS at `k = 10`
  and 2.62 GB at `k = 30` on 10^6 windows, so `mem_mb = 3000 + 250 * k` leaves
  roughly a 4x margin. `mem_mb` is a global concurrency cap, so every GB
  over-reserved is a fit that cannot start alongside.

# 05 Phase 2 including repetitive -- the same comparison with nothing filtered out

Belongs with section 03; it is numbered 05 only so the existing sections keep
their numbers. Phase 2 is **not modified** and both results stand side by side.

Phase 2 draws its figures from windows that survive two filters: the mappability
blacklist, and the `EXCLUDE_CHROMOSOMES` contig list. Regions land on the
blacklist *because* they are repetitive, and `bwa mem` assigns their
multi-mapping reads to one of the alternatives at random rather than discarding
them -- so the reads are there, the motifs are there, and the filter removes
exactly the repeat content this comparison is about.

## The manuscript does not filter here either

`04_.../compile_reads_vs_nonb.R` never opens `blacklist_files/`. Its only filter
is `!chr %in% c("chrX","chrY")`. The blacklist first appears at step 09, the
HMM (`rule Remove_human_black_list`, `bedtools intersect -v`) -- which is where
Phase 3 applies it here. So the blacklist belongs at the HMM, not at this
correlation, and Phase 2 was over-filtering.

## Why the boolean filter is worse than its window size suggests

The blacklist is built on **10 kbp** windows in both pipelines (`BLACKLIST_WINDOW`
here, hardcoded `10000` in the manuscript's `window_coverage_for_blacklist`),
with the same rule: ratio = `All_reads / High_MapQ` per window, `Inf`/`NA` set to
2, capped at 2, keep ratio >= 1.1. Then `bedtools merge -d 1` merges adjacent
flagged windows into long intervals, and the sex chromosomes are appended whole.

That merge is the part that bites. Measured on the manuscript's own
`Hsapien.bed` against its `Hsapien.chrom.sizes`: 9,500 intervals, mean length
65 kbp, 403 Mbp on autosomes -- about 14% of autosomal sequence. Autosomal
windows touched by at least one interval, i.e. dropped by Phase 2's
`blacklist_window_flag` boolean:

| Window | Touched | Retained by Phase 2 |
|---|---|---|
| 1 Mbp | 93.0% | **7%** |
| 100 kbp | 37.0% | 63% |
| 10 kbp | 13.9% | 86% |

So `plots/phase2/nonB_scatter.svg`, which draws 1 Mbp windows with
`!Blacklisted`, is fitting the cleanest few percent of the genome, and the 1 Mbp
column of `plots/phase2/nonB_correlation_heatmap.svg` has the same problem. This
phase scores each window by the **fraction** of its bases that are blacklisted
instead, which degrades gracefully where the boolean does not.

## Nothing expensive is recomputed

Phase 2's filtering happens inside its R scripts, so the files underneath it
already cover the whole genome. This phase reuses them untouched:

| Reused as-is | Rule that built it |
|---|---|
| `data/phase2/experiment_window_coverage/All_reads/{size}_nucleotides/{experiment}.bed` | `experiment_window_coverage` |
| `data/phase2/nonb_density/{size}_nucleotides/{nonb}.bg` | `nonb_density` |
| `results/phase2/Blacklist.bed` | `add_sex_chromosomes_to_blacklist` |
| `data/phase1/genome_windows/{size}_nucleotides/<genome>.bed` | `make_genome_windows` |

No BAM is read, no `megadepth` pass, no motif intersect. **Seven jobs**: one
`bedtools coverage` per window size, one R compile, three plots. `Rutility.yml`
and `Rplot.yml` are reused **unmodified** and there is no new conda env, so
nothing already finished rehashes and reruns.

## Two filters become two columns, not two decisions

`Region_set`, from the blacklisted fraction *f* of each window:

| Level | Definition |
|---|---|
| `All_windows` | every window -- **what the figures show** |
| `Blacklist_free` | *f* == 0, i.e. exactly Phase 2's `Blacklist_filtered == TRUE` |
| `Blacklist_light` | *f* <= `REPEAT_LIGHT_FRACTION` (0.05) |
| `Blacklist_heavy` | *f* >= `REPEAT_HEAVY_FRACTION` (0.5) |

`Chromosome_set`:

| Level | Definition |
|---|---|
| `All_contigs` | everything: sex chromosomes, `chrMT` and the rDNA morphs included -- **what the figures show** |
| `Phase2_contigs` | `EXCLUDE_CHROMOSOMES` removed, i.e. Phase 2's contig set |

`All_windows` x `Phase2_contigs` reproduces the manuscript's step 04 convention
exactly, so these numbers stay directly comparable with the published human ones.

The figures take "all regions of the genome" literally, which means the three
`rDNA_morph_*` contigs -- three assemblies of one locus -- each contribute. If
you want them out of the figures, that is one `filter()` line in the plot
scripts, and `Phase2_contigs` is already in the CSV either way.

CPM is normalised genome wide over every window *before* any subsetting, exactly
as Phase 2 does, so CPM means the same number in both phases.

## Run it

Phase 2 must have finished. Nothing new to install.

```bash
conda activate js4025_snakemake

# Should report 7 new jobs and nothing from Phase 1, 2 or 3.
snakemake --snakefile workflow/Snakefile -n --cores 5 --use-conda

mkdir -p slurm_logs && sbatch submit_phase2_including_repetative.sh
```

`submit_phase2_including_repetative.sh` requests 1 node, 5 cores, 4 GB/core on
`basic` for 2 h. Five cores is the **floor**, not a preference: at 4 GB/core it
is the smallest allocation whose derived `mem_mb` cap clears the 16000 MB that
`compile_reads_vs_nonb_all_regions` requests. More cores buy nothing.

Unlike the other submit scripts this one names its six targets explicitly rather
than building `rule all`, so submitting it cannot start a 48 h Phase 3 HMM fit by
accident.

## Outputs

| File | What it answers |
|---|---|
| `results/phase2_including_repetative/nonB_correlation_coefficient.csv` | 3456 rows: experiment x 8 classes x 3 windows x 2 metrics x 2 signals x 4 region sets x 2 chromosome sets |
| `results/phase2_including_repetative/nonB_coverage_versus_read.csv` | per-window table at 1 Mbp, **every** window, carrying `blacklist_fraction` and `Phase2_excluded_contig` |
| `results/phase2_including_repetative/Region_composition.csv` | windows and Mbp per region set x chromosome set x window size -- how much of the genome Phase 2 was discarding |
| `plots/phase2_including_repetative/nonB_correlation_heatmap.svg` | the Phase 2 heatmap, on all windows and all contigs |
| `plots/phase2_including_repetative/nonB_correlation_delta.svg` | **the figure for this question**: what including the repeats does to each correlation |
| `plots/phase2_including_repetative/nonB_scatter.svg` | per-window scatter with the blacklisted fraction encoded rather than filtered |

### Reading the delta heatmap

Each cell is rho(`All_windows`) minus rho(`Blacklist_light`), so:

- **near zero** -- the blacklist was irrelevant to that correlation, and Phase 2's
  number stands
- **away from zero** -- the repeats were carrying signal Phase 2 removed, and the
  sign says which way. Positive means including them raises the correlation.
- **grey `n/a`** -- fewer than 30 windows survive the filtered side, so there is
  nothing to difference against. Expect this in the 1 Mbp column.

`Blacklist_light` rather than `Blacklist_free` is the comparison because
`Blacklist_free` is a tiny, strongly biased slice at coarse window sizes (see the
table above). `Blacklist_free` is still in the CSV, where it is exactly Phase 2's
filtered branch.

### Reading the scatter

One row per series, one column per motif class, every window drawn:

- **colour** -- fraction of the window that is blacklisted, on viridis
- **shape** -- triangles are contigs Phase 2 excluded outright (`chrW_mat`,
  `chrZ_pat`, `chrMT`, the rDNA morphs), so they are identifiable as outliers
  rather than anonymous
- **solid line** -- `lm` over all windows in the panel
- **dashed line** -- `lm` over windows at most `REPEAT_LIGHT_FRACTION` blacklisted

If the dashed line sits on the solid one, the repeats are not driving the
correlation. If it swings away, they are.

**A reading note on the y axis.** `megadepth --op sum` reports summed per-base
coverage, which is not length-normalised, so a short window has a proportionally
small CPM. Every contig's final window is short, and `chrMT` is one 16 kbp
window, so a handful of points sit low and stretch the axis downward. Phase 2 and
the manuscript share this property; it is not introduced here, and CPM was left
alone so the two phases stay comparable.

## Sanity checks after the run

1. `Region_composition.csv` first. At `Phase2_contigs` x `Blacklist_free` it says
   what fraction of windows Phase 2 was keeping. If that is not far below 100% at
   1 Mbp, either the blacklist is much smaller than the human one or the
   fraction rule did not run.
2. `Pct_windows` and `Pct_Mbp` should diverge -- windows are counted whole, bases
   are not. If they are identical everywhere, the fraction column is boolean and
   `bedtools coverage` wrote the wrong thing.
3. In `nonB_correlation_coefficient.csv`, `All_windows` should sit **between**
   `Blacklist_free` and `Blacklist_heavy` for any class where the two differ.
   It is a mixture of them; landing outside both means something is misaligned.
4. The WGS control should still sit near zero on most classes in the
   `All_windows` branch. It is the negative control for the whole comparison. If
   it now tracks motif density as strongly as PDAL-Seq does, the signal that
   appeared when the repeats came back is compositional (GC, mappability), not
   permanganate reactivity -- and the delta heatmap's WGS row is where that shows.
5. Compare `plots/phase2_including_repetative/nonB_correlation_heatmap.svg`
   against `plots/phase2/nonB_correlation_heatmap.svg` cell by cell. The delta
   figure is that comparison done arithmetically.

## Verification done before this phase shipped

- `bash -n` on the submit script, R `parse()` on all four scripts.
- Dry run against a mock tree holding only the four Phase 2 / Phase 1 inputs
  above: resolved to **exactly 7 jobs**, reaching for no `megadepth`, no
  `samtools`, no `nonb_density`, and nothing in Phase 3.
- End-to-end on synthetic data with a **planted contrast** -- 53 Mbp over seven
  contigs, a blacklist of 62 merged intervals, and a G4/coverage relationship
  deliberately made positive outside the blacklist, flat in between, and negative
  inside it. Recovered: `Blacklist_free` rho = **+0.997**, `Blacklist_heavy` rho =
  **-0.998**, `All_windows` **+0.392** (the mixture, between them, as it must be).
  Only the planted class moved (all others |rho| < 0.035); WGS and 0mM stayed
  flat; `Ratio_vs_WGS` appeared only for PDAL-Seq and dropped exactly the 55
  windows where the control had zero coverage. All 16 rows of
  `Region_composition.csv` matched an independent `awk` recount, window counts
  and Mbp. Row counts came out at the expected product exactly.
- All three figures rendered and inspected, plus four edge cases: a filtered
  branch too thin at one window size (draws grey `n/a`, keeps the rest), too thin
  at every window size (fails loudly rather than emitting a blank figure), NA
  cells in the heatmap (drawn `n/a`, not dropped), and no window qualifying for
  the scatter's comparison line (draws, reports `0 of 16 panels`).

Not checkable locally: behaviour at real genome scale.

## Decisions recorded for this phase

- **Phase 2 is untouched.** No existing rule, script, env or constant was
  modified, so nothing already finished reruns and both analyses remain
  available. The cost is four near-duplicate R scripts, which is the right
  trade -- editing the Phase 2 scripts in place would have rehashed their rules.
- **Fraction, not boolean.** Phase 2 asks "does any blacklist interval touch this
  window?". After `merge -d 1` that question drops 93% of 1 Mbp windows to
  remove 14% of the sequence. This phase asks what fraction of the window is
  blacklisted, which is the same question at 10 kbp and a much better one above.
- **Both filters are columns, not decisions.** Four region sets and two
  chromosome sets, so the CSV answers every reading of "all regions" -- including
  the manuscript's exact step 04 convention, which is
  `All_windows` x `Phase2_contigs`.
- **The figures take "all regions" literally**, rDNA morphs included, rather than
  quietly reinstating a contig filter the request did not ask for.
- **CPM is unchanged.** Normalised genome wide before subsetting, not
  length-normalised, exactly as in Phase 2 and the manuscript. Fixing the
  short-window property would have broken comparability with both, so it is
  documented instead.
- **The delta compares against `Blacklist_light`, not `Blacklist_free`.** A
  difference against a branch that is 7% of the genome at 1 Mbp is not
  interpretable. `Blacklist_free` stays in the CSV as the exact Phase 2 branch.
- **`Blacklist_heavy` exists** so "do the repeats behave differently?" can be
  answered directly rather than inferred from the difference between two
  overlapping sets.
- **The submit script names its targets** instead of building `rule all`, so it
  cannot start Phase 3's 48 h HMM fits by accident.

# 06 Phase 4 -- Intersect HMM states with non-B DNA motifs, 5mC rates, RNA-Seq and functional genomic annotations

Take the Gaussian HMMs Phase 3 fit and ask **what the states mean**. Nothing is
re-fit here: this phase reads `data/phase3/GHMM/<k>/GHMM_states.csv.gz` and
`Parameters.csv` and measures, for every state of every model, how enriched it
is for non-B DNA motifs, functional genomic elements and RNA-Seq reads, and how
methylated it is.

Adapted from `07_PhyloHGMP_model.smk` (state annotation, `Enrichment.sh`,
`BG_enrichment.sh`, the compile pattern), `10_Methylation_analysis.smk` (5mC in
states), `18_ATACseq_analysis.smk` (the read mapping and filtering protocol) and
`Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R` (the plot layout).

**Seven models are analysed, not one.** `PHASE4_STATES = 8 … 14` brackets the
elbow of the BIC scree, so every candidate *k* gets the full treatment and the
choice is made from what the states turn out to mean rather than from BIC
alone. That is what makes this phase 1,539 jobs instead of 220.

## Inputs this phase expects

Everything below lives on Roar and is described in `resources/README.md`.

| Path | What it is |
|---|---|
| `resources/RNA-Seq/SRR17849680_{1,2}.fastq.gz` | CFS414 RNA-Seq, passage 50 |
| `resources/RNA-Seq/SRR17849681_{1,2}.fastq.gz` | CFS414 RNA-Seq, passage 24 |
| `resources/gene_annotation/bTaeGut7v0.4_MT_rDNA.matZ.{promoter,UTR5,CDS,introns,UTR3,lncrna,intergenic}.bed` | functional elements, longest transcript |
| `resources/centromere/bTaeGut7v0.4_MT_rDNA.matZ.CEN.bed` | centromeres, one per chromosome |
| `resources/methylation_data/bTaeGut7v0.4_MT_rDNA.matZ.PBmethylation.v0.1.bed` | PacBio 5mC, `chr start end percent` |
| `results/phase2/nonb_annotations/*.bed` | the seven cleaned non-B classes plus `all`, from Phase 2 |

All the annotation files are on the single haplotype already -- `chr*_mat` plus
`chrZ_pat` -- so their contig names match our genome and `clean_annotation.sh`
only has to drop `chrMT` / rDNA and clamp to the contig ends.

The repeat annotations (EDTA2, Satellites, TRF_withMers) and the A/B compartment
track that Linnéa listed are **not** used here. `resources/README.md` only
records the gene annotations and the centromeres being copied, and the
satellites are Phase 5's subject.

## RNA-Seq is mapped with STAR, not bwa mem

`18_ATACseq_analysis.smk` uses `bwa mem`, which is correct for ATAC-Seq --
those fragments are genomic. RNA-Seq reads cross splice junctions, and bwa mem
would soft-clip every junction-spanning read. Everything else in that protocol
is carried over unchanged: the same `fastp` flags (`-p -P 50 -y -x -3`), the
same "count reads in intervals" downstream.

The one other change is that the deduplication step is dropped; see below.

Three consequences of the aligner swap worth knowing:

- **STAR encodes MAPQ as 255 / 3 / 1 / 0**, not bwa's 0-60. Nothing in Phase 4
  filters on MAPQ (the ATAC protocol does not either), but the `--min-MQ 20`
  idiom from Phase 1 means something different on these BAMs.
- **STAR streams unsorted BAM into `samtools sort`** rather than writing the
  intermediate SAM that `ATACseq_bwa_mem` -> `ATACseq_sort_sam` writes. At
  ~460 M read pairs that file would be several hundred GB. Phase 1's
  `map_reads.sh` already made the same call.
- **Junctions are found de novo.** `resources/` has the BED files Linnéa
  derived from the gene annotation but not the GFF itself. To use
  annotation-guided junctions instead:

```
cp /storage/group/kdm16/default/lbs5874/ZebraFinch/ref/bTaeGut7v0.4_MT_rDNA.matZ.gff \
   resources/gene_annotation/
```

  then point `RNA_SEQ_GTF` at it in `workflow/Snakefile`. At the 1 kbp window
  resolution this signal is read at, it will not move the numbers.

### There is no duplicate-removal step

`ATACseq_remove_dup` runs `samtools rmdup -s`, which discards reads sharing an
alignment start. On genomic fragments those are PCR duplicates. On RNA-Seq they
are mostly **not**: a highly expressed gene generates many independent
fragments from the same start position, so rmdup strips the most reads from
exactly the regions this analysis is meant to detect, compressing the dynamic
range of the expression signal in a depth-dependent way.

The step was in the first version of this phase, on the grounds that the
protocol being followed does it. It was removed on 2026-09-04. `rna_signal` and
`rna_flagstat` now read `data/phase4/RNA_seq/sorted/` directly, and that BAM is
no longer `temp()` -- with nothing downstream of it, it is the analysis BAM.

To put it back: restore a `rna_dedup` rule running
`samtools rmdup -s` from `sorted/` to `dedup/`, point those two rules at
`dedup/`, and re-`temp()` the `rna_map` output.

## The enrichment calculation

Exactly the manuscript's, from the methods:

> each annotation type was first subsampled at random to 10,000 annotations …
> the length in base pairs of the intersection between the HMM states segments
> and the subsampled annotations was recorded … the subsample was shuffled to
> random genomic coordinates -- but keeping annotations on the same chromosome
> … This process was repeated 100 times, and enrichment was calculated by
> dividing the mean of the first intersect by the mean of the second, shuffled
> intersect.

**The shuffle is genome wide.** `bedtools shuffle -g chrom.sizes -chrom`, with
no `-incl`, as in the manuscript's `Enrichment.sh`. Our states only cover the
segmented genome (autosomes, blacklist removed), so shuffled annotations can
land in regions no state covers and the null is diluted by that fraction. The
consequence is that these ratios are enrichment **relative to a whole-genome
expectation**, not relative to a within-segmentation expectation, and every
state will look slightly more enriched than it would under the latter.
Comparisons *between* states and *between* annotations are unaffected, because
the dilution is the same for all of them. `results/phase4/Segmented_genome.bed`
records the covered fraction if you want to correct for it after the fact.

### Adjacent windows are merged into segments

The manuscript's states were alignment blocks and were already segment-like.
Ours are a fixed 1 kbp grid, so `Annotate_states.R` merges adjacent windows in
the same state into a segment before anything intersects them. Without it a
single state is ~10⁵ one-window records and `sum += $3 - $2 + 1` double counts
at every window boundary.

### Both significance tests are reported

The manuscript's compile scripts call `wilcox.test` unconditionally, but the
methods say a Wilcoxon test only applies above 10,000 loci and that smaller
annotations get an empirical p-value from the null distribution. That
distinction matters here: **the centromere annotation has 41 records**, so
`shuf -n 10000` returns all of them every iteration, `Observed` is constant,
and a Wilcoxon test is comparing a point mass against a distribution.

`Compile_enrichment.R` therefore reports `Wilcox.p`, `Empirical.p`,
`N_annotation`, and a `P_value` column carrying whichever one the methods
select, with `Test` saying which. The plots mark significance off `P_value`.

### Two `-u` flags the manuscript scripts do not have

`BG_enrichment.sh` and `methylation_in_state.sh` add `bedtools intersect -u`.
Without it a record straddling the boundary between two segments of the *same*
state is emitted twice and counted twice. The manuscript's inputs were raw
`genomecov` bedGraphs with few-bp records, so straddling was rare; ours are on
the 1 kbp grid the states are built from, where it is common. This is a
counting fix, not a change of analysis.

### The RNA-Seq signal is binned to 1 kbp before enrichment

`BG_enrichment.sh` re-reads its input once per iteration. A raw
`bedtools genomecov -bg` for a 464 M pair library is on the order of 10⁸
records, so 100 iterations x 154 jobs would be roughly a petabyte of I/O for a
signal being read at 1 kbp resolution anyway. `rna_signal` instead runs
`megadepth --op sum` over the same window BED the HMM was fit to, giving
~1.15 M records. `RNA_ENRICHMENT_SUBSAMPLES` is 100,000 rather than the
manuscript's 1,000,000 for the same reason: 1,000,000 out of 1.15 M would draw
87% of the data every iteration and leave the null with almost no variance.

## Configuration

Set at the top of `workflow/Snakefile`:

| Constant | Value | Note |
|---|---|---|
| `PHASE4_STATES` | 8 … 14 | one full analysis per *k* |
| `RNA_SRA` | `SRR17849680`, `SRR17849681` | passage 50 and passage 24 |
| `RNA_SEQ_GTF` | `""` | de novo junctions; see above |
| `RNA_SEQ_READ_LIMIT` | `"0"` | `fastp --reads_to_process`; the wall-clock lever |
| `ANNOTATIONS` | 7 gene features + `CEN` | |
| `ENRICHMENT_SUBSAMPLES` | `"10000"` | manuscript value |
| `ENRICHMENT_ITERATIONS` | `"100"` | manuscript value, deliberately not higher |
| `RNA_ENRICHMENT_SUBSAMPLES` | `"100000"` | see above |

The methods are explicit that 100 iterations is a ceiling, not a budget: *"we
used 100 iterations of subsampling because considering a larger number of
subsamples resulted in all enrichments being significant due to the large size
of the datasets."* Do not raise it.

## Run it

```
rsync -Pav ./ submit.hpc.psu.edu:/storage/group/kdm16/default/jus841/js4-2026/js4025_Analyze_Zebra_finch_PDAL-Seq_data

module load anaconda
conda activate js4025_snakemake

# Should report 1539 new jobs and nothing from Phase 1, 2 or 3.
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda

mkdir -p slurm_logs && sbatch submit_phase4.sh
```

`submit_phase4.sh` uses **exactly the header submit_phase3.sh uses** -- same
partition, nodes, ntasks, mem-per-cpu and 48 h, and no `--account` line. Phase 3
submitted and ran with it. A 96 h version of this job was rejected twice with

```
JOBID     PARTITION NAME     USER   ST TIME NODES NODELIST(REASON)
55341276  basic     js4025_p jus841 PD 0:00 1     (AssocGrpBillingMinutes)
```

and the cause was never established -- `open`/`jus841` showed 6,882,469 of
9,300,096 billing minutes still free at the time, which should have covered the
115,200 the request needed. If a future phase hits the same wall, the thing
that is known to work is this header, unchanged.

It names its Phase 4 targets rather than building `rule all` -- `rule all`
spans every phase, so anything snakemake decided was stale in Phases 1-3 would
restart mapping or a 4 h HMM fit inside this job.

### Restarting is safe

Phase 4 is ~35-40 h of work in a 48 h reservation, so it should finish in one
go. If it does hit the wall clock, just resubmit. Every stage is a separate
snakemake target; the STAR BAMs are kept rather than `temp()`; and
`Enrichment.sh`, `BG_enrichment.sh` and `methylation_in_state.sh` write to
`<output>.partial.$$` and rename into place only once complete, so a killed job
leaves **no** output file rather than a truncated one. That matters because a
`SIGKILL`ed snakemake never records the job as incomplete, so
`--rerun-incomplete` would not catch a short CSV -- it would be read as a
finished result with fewer iterations in it.

```
sbatch submit_phase4.sh
# ... when it ends:
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda
# resubmit until that reports nothing to do
```

### Disk

This is the heaviest phase for scratch. Plan for **~320 GB free**:

| | |
|---|---|
| STAR index | ~11 GB, kept in `resources/genomes/STAR_index/` |
| trimmed fastqs | ~160 GB, `temp()` |
| sorted RNA-Seq BAMs | ~140 GB, kept -- the analysis BAMs |
| everything else | ~2 GB |

Snakemake deletes the temp files as soon as their consumers finish, but the two
libraries are processed independently and can overlap.

## Outputs

Per *k* in `results/phase4/<k>/` and `plots/phase4/<k>/`:

| File | What it answers |
|---|---|
| `HMM_functional_summary.svg` | the Figure 3 composite: signal, size, non-B, functional, RNA-Seq, 5mC |
| `State_distribution_by_chromosome.svg` | where each state sits in the genome |
| `State_summary.csv` | per state: segments, bp, percent of segmented genome, fitted and empirical mean |
| `State_by_chromosome.csv` | bp and percent per chromosome per state, plus the unsegmented remainder |
| `nonB_enrichment.csv` | 8 motif classes x *k* states, both p-values |
| `Functional_enrichment.csv` | 8 annotations x *k* states, both p-values |
| `RNA_enrichment.csv` | 2 libraries x *k* states |
| `Methylation_by_state.csv` | CpG counts and percents in five 5mC rate bins |

Once, not per *k*:

| File | What it answers |
|---|---|
| `results/phase4/Segmented_genome.bed` | what the HMM could have covered -- the denominator |
| `results/phase4/Annotation_report.txt` | how many records of each annotation survived the contig filter |
| `results/phase4/RNA_seq_mapping_statistics.txt` | STAR `Log.final.out` and flagstat per library |
| `data/phase4/states/<k>/state_<i>.bed` | the state segments themselves |

### Reading the composite figure

Six panels sharing one state axis, ordered by mean PDAL-Seq signal with the
highest at the top.

- **States are labelled `s<i>` and the ordering is recomputed per model.** `s3`
  in the *k* = 9 panel is not `s3` in the *k* = 12 panel. Read the state id, not
  the row position.
- **The heatmap colour is clipped at ±log2(3); the printed number is not.** A
  tile can read `12.4` on a scale that stops at 3, exactly as in the manuscript
  figure.
- **`ns` / `*` / `**` come from `P_value`**, which is the empirical test for
  `CEN` and the Wilcoxon test for everything else.
- **A grey tile** is a non-finite enrichment -- an unoccupied state, or a null
  that came out zero.
- **The `unsegmented` band** in the distribution plot is the blacklist plus the
  contigs Phase 3 excluded outright, which is why `chrW`, `chrZ`, `chrMT` and
  the rDNA morphs are solid grey.

### Unoccupied states

A variational HMM is free to leave a state empty, and Phase 3's fit plots show
that some *k* do. Phase 4 does not treat that as an error: `Annotate_states.R`
writes an empty BED, the enrichment scripts short circuit to a table of zeros,
and the compile steps report `NaN` enrichment and `NA` for the Wilcoxon p. A
model with many empty states is telling you *k* is too high -- that is a result,
not a failure.

## Sanity checks after the run

```
# 1. Did every annotation survive the contig filter? The gene annotations
#    should keep ~100% (they are already on our haplotype); CEN should keep 41.
cat results/phase4/Annotation_report.txt

# 2. Did STAR map a sensible fraction? Look for "Uniquely mapped reads %".
#    Anything below ~70% on a matched cell line means the index or the reads
#    are wrong, not the biology.
grep -A2 "Uniquely mapped reads" results/phase4/RNA_seq_mapping_statistics.txt

# 3. Do the state sizes add up to the segmented genome?
awk -F, 'NR>1 {s += $4} END {print s}' results/phase4/8/State_summary.csv
awk '{s += $3 - $2 + 1} END {print s}' results/phase4/Segmented_genome.bed

# 4. Do the per-chromosome percents reach 100?
awk -F, 'NR>1 {s[$1] += $6} END {for (c in s) printf "%s\t%.1f\n", c, s[c]}' \
    results/phase4/8/State_by_chromosome.csv

# 5. How many states did each model actually occupy?
for K in 8 9 10 11 12 13 14; do
    echo -n "k=$K occupied: "
    awk -F, 'NR>1 && $4 > 0' results/phase4/$K/State_summary.csv | wc -l
done

# 6. Did any enrichment job leave a temp directory behind (i.e. was killed)?
ls -d temp/enrichment_* temp/bg_enrichment_* 2>/dev/null | wc -l
```

## Verification done before this phase shipped

The whole chain was run end to end on synthetic fixtures -- a 3-contig genome,
a *k* = 3 model, a 20,000-record annotation and a 2-record one, a methylation
BED and an RNA-Seq bedGraph -- with the real scripts and real `bedtools`:

- segment merging collapses 146 windows into 23 segments and the per-state bp
  sums back to the segmented total;
- `clean_annotation.sh` clamps a record running past the contig end and drops a
  record on a contig that is not in the genome;
- the 2-record annotation produces a constant `Observed`, is correctly routed
  to the empirical p-value, and the 20,000-record one to the Wilcoxon;
- an empty state BED short circuits to a table of zeros instead of failing;
- every enrichment job removes its temp directory on exit;
- the per-chromosome percents sum to exactly 100 including the unsegmented row;
- both figures render.

A full `--dry-run` against snakemake 8 resolves the DAG with no ambiguity or
missing input and reports the expected 1,539 Phase 4 jobs: 616 non-B, 616
functional, 154 RNA-Seq and 77 methylation intersections, 7 of each compile and
plot, 8 annotation cleans, 7 state annotations and the RNA-Seq mapping chain.

## Decisions recorded for this phase

- **Seven models, not one.** The functional analysis is what distinguishes a
  *k* that splits a real biological state from a *k* that splits noise, so it
  runs for every candidate rather than after a choice has been made.
- **STAR, not bwa mem**, for RNA-Seq -- see above. Everything else in the ATAC
  protocol is unchanged.
- **No duplicate removal.** The ATAC protocol's `rmdup -s` is wrong for
  RNA-Seq -- it strips the most reads from the most transcribed regions. It was
  in the first version of this phase and was removed on 2026-09-04; the sorted
  BAM became the analysis BAM and is no longer `temp()`.
- **The shuffle is genome wide**, as the manuscript does it, with the dilution
  it implies written down rather than corrected for.
- **Adjacent windows are merged into segments** before intersecting. The
  manuscript did not need to; a 1 kbp grid does.
- **Both p-values are reported.** The manuscript's compile scripts only ever
  computed the Wilcoxon; the methods describe both, and `CEN` needs the other.
- **Coordinates are carried through, not re-joined by row position.** The
  manuscript's `Annotate_state.R` pairs a BED and a state CSV by row order;
  Phase 3 already carries `chr/start/end` into `GHMM_states.csv.gz`.
- **`Compile_enrichment.R` is one script, not three.** The manuscript's three
  compile scripts differ only in a directory and a column name.
- **Two new environment files** (`phase4_functional.yml`, `Rphase4.yml`). No
  existing yml was edited -- snakemake hashes them, and an edit would re-queue
  every finished Phase 1-3 job that used one.
- **The submit script names its targets** rather than building `rule all`, so
  it cannot restart a Phase 3 HMM fit by accident.
- **The SLURM header is copied from Phase 3 verbatim.** A 96 h request was
  rejected twice on `(AssocGrpBillingMinutes)` for reasons that were never
  pinned down; the known-good header is not worth improvising on.
- **Enrichment outputs are written atomically.** A restart-driven workflow
  cannot afford a truncated CSV being mistaken for a finished one.

# 07 Phase 5 -- PDAL-Seq versus control read density in satellites, centromeres and transposable elements

Ask whether PDAL-Seq reads pile up in the repetitive genome, and control for
the obvious alternative explanation by asking the same question of the PCR-free
library.

Adapted from `17_PDALseq_reads_in_CenSat.smk` and its
`Read_enrichment.sh` / `Compile_CenSat_Class_enrichment.R`.

**This phase filters nothing.** No blacklist, no MAPQ threshold, no chromosome
exclusion, whole genome, all reads. Every filter this pipeline has removes
repetitive sequence, and repetitive sequence is the subject here. The coverage
source is `data/phase1/bigwig/<sample>.bigwig`, which Phase 1 built from the
deduplicated `All_reads` BAM.

**All six libraries are analysed separately**, not merged. The two 20 mM and
two 40 mM replicates agreeing on a class is the evidence that the class is
really enriched, and merging would throw that away. The 0 mM library is the
internal negative control -- the PDAL-Seq protocol without permanganate.

## Inputs this phase expects

| Path | What it is |
|---|---|
| `resources/repeats/bTaeGut7v0.4_MT_rDNA.matZ.EDTA2.v0.2.bed` | EDTA2 repeats, class in **column 5** |
| `resources/repeats/bTaeGut7v0.4_MT_rDNA.matZ.Satellites.bed` | satellite families, name in **column 4** |
| `resources/centromere/bTaeGut7v0.4_MT_rDNA.matZ.CEN.bed` | centromeres, one per chromosome (already used by Phase 4) |
| `data/phase1/bigwig/<sample>.bigwig` | per-base read coverage, from Phase 1 |
| `data/phase1/genome_windows/1000_nucleotides/<genome>.bed` | the 1 kbp grid, from Phase 3 |

83 classes are analysed: 26 EDTA2 classes, 54 satellite families, the
centromeres, and two pooled classes (`all_TE`, `all_Satellite`). The class lists
are **hardcoded in `workflow/Snakefile`** from the counts in
`resources/README.md`, so the whole DAG resolves before anything runs.
`split_repeats.py` checks what it finds against those lists and stops with the
offending names if they have drifted -- a class that appeared or was renamed
upstream would otherwise surface as an unexplained missing-output error.

### Three things `resources/README.md` raises, answered

- **`TRF_withMers.bed` is skipped**, as you decided. For the record, its two
  label columns are a Tandem Repeat Finder period *bin* and the exact period of
  the repeat unit: `11-50mer  43-mer` means TRF called a 43 bp unit there. So
  its classes are unit **lengths**, not repeat families, and they cut across
  the satellite families rather than adding a class of their own. Adding it
  later is one more split rule on column 4 plus its own class list; nothing
  else changes.
- **`resources/README.md` says the satellite file has 50 classes.** The listing
  under that sentence has 54 (`TEL` plus 53 `Tgut*`), and 54 is what the
  Snakefile declares. If the split rule fails with a missing class, that is the
  discrepancy talking.
- **Five EDTA2 "classes" are not transposable elements** -- `low_complexity`,
  `repeat_fragment`, `repeat_region`, `target_site_duplication` and
  `rRNA_gene`. They are kept. They are annotated territory that PDAL-Seq reads
  can land in, and dropping a class from a comparison is a decision better made
  after seeing it than before.

## The enrichment calculation

Per library, per class, 100 iterations of:

1. subsample 1,000 annotation fragments at random;
2. `Observed` = summed per-base read coverage over them;
3. shuffle them to random coordinates **on the same chromosome**;
4. `Null` = summed per-base read coverage over the shuffle.

Enrichment is `mean(Observed) / mean(Null)`. Library depth cancels in that
ratio, which is what lets the PCR-free control be compared to a PDAL-Seq
library directly.

### What a "fragment" is

The class BED is first clipped to the 1 kbp window grid with
`bedtools intersect -a windows -b class`. Without `-wa` that reports the
**overlap**, so a fragment is *the part of one annotation that falls in one
1 kbp window*. This is the manuscript's behaviour and it is kept deliberately:

- the 1 kbp cap bounds the length of the sampled intervals, which is what lets
  `bedtools shuffle` place them and what stops one 300 kbp satellite array from
  being a single sampling unit;
- a window covered by many elements of a class yields many fragments, so the
  subsample is weighted by annotation density -- the methods' *"each annotation
  type was first subsampled at random"*, sampling annotations rather than
  windows. Merging the class first, or adding `-u`, would silently turn that
  into uniform sampling over windows.

The one change is that the intersection is computed **once**, not inside every
iteration as the manuscript's script does. Same file, 100 times less work.

### The coverage comes from a bigWig, not the BAM

This is the one substantive deviation from
`17_PDALseq_reads_in_CenSat/Read_enrichment.sh`, and it is a deviation in
*route*, not in statistic.

megadepth reads a BAM front to back -- there is no index seek -- so one
iteration of the manuscript's loop costs a full scan of a 5 GB alignment file.
The loop makes two megadepth calls and this phase has 498 jobs, so the
manuscript's version is ~100,000 whole-BAM scans and would not finish. A bigWig
is indexed, so the same query is a range read.

Phase 1 built `data/phase1/bigwig/<sample>.bigwig` with
`bedtools genomecov -bg -split -ibam` on exactly the `All_reads` BAM this phase
is asking about, so `megadepth <bw> --op sum` and `megadepth <bam> --op sum` are
both "summed per-base read depth over these intervals".

`Read_enrichment.sh` treats empty megadepth output as a hard error rather than
as zero coverage -- an interval with no reads under it still prints a line with
a `0`, so empty stdout means megadepth failed, and 498 files of zeros is not a
failure mode worth allowing. If megadepth's bigWig mode ever stops printing to
stdout, that check is what will say so on the first job.

### Significance: the threshold is the subsample size

The methods pick the test by annotation count: Wilcoxon above 10,000 loci,
empirical below. In the manuscript, 10,000 was *also* the subsample size, and
the subsample size is the reason for the rule -- at or below it, `shuf -n`
returns the whole file every iteration, so `Observed` is a constant and a
Wilcoxon test compares a point mass against a distribution. Phase 5 subsamples
1,000, so **1,000 is the threshold**.

Here that is the normal case, not the exception: most satellite families have a
handful of annotations. Both p-values are in the output either way, and `Test`
says which one `P_value` carries.

## Outputs

| Path | What it is |
|---|---|
| `results/phase5/Repeat_annotation_report.txt` | records, contigs and bp per class -- **read this first** |
| `results/phase5/Repeat_enrichment.csv` | one row per library x class: enrichment, both p-values, fragment count |
| `results/phase5/PDALSeq_versus_control.csv` | the same divided by the PCR-free control |
| `plots/phase5/Repeat_enrichment_heatmap.svg` | log2 enrichment, classes x libraries, PCR-free column first |
| `plots/phase5/PDALSeq_versus_control.svg` | log2(PDAL-Seq / PCR-free) per class, one point per library |

### Reading them

`PDALSeq_versus_control.svg` is the result. `Repeat_enrichment.csv` is what it
is built from.

A class can look enriched in PDAL-Seq for two reasons: permanganate-reactive
ssDNA is there, or the class is a mapping artefact -- multi-mapping reads that
bwa mem scattered into it, or a sequencing bias in a GC-extreme array. The
PCR-free control has the second and not the first, so

```
Ratio = Enrichment(PDAL-Seq library) / Enrichment(PCR-free library)
```

is what survives. **Zero on that axis means the class attracts PDAL-Seq reads
exactly as much as it attracts PCR-free reads**, i.e. everything about it is
mapping and sequencing. Above zero is enrichment over and above that.

Three things to check before believing a row:

1. **The fragment count**, carried on the axis label and in the CSV. A tenfold
   enrichment resting on three fragments is one number, not a mean over 100
   draws -- `Test` will say `empirical` for exactly those rows.
2. **Where the 0 mM library sits.** It went through the PDAL-Seq protocol
   without permanganate. Where it tracks the treated libraries, the signal is
   protocol, not permanganate.
3. **Whether the replicates agree.** Two 20 mM and two 40 mM points on top of
   each other is the evidence; a spread across the panel is not.

Both figures order classes the same way -- by mean log2 ratio within each panel
-- so they can be read side by side. Fill and axis are clipped at ±log2(8);
satellite arrays reach well outside that and the unclipped numbers are in the
CSVs.

## Configuration

Everything is in the Phase 5 block of `workflow/Snakefile`:

| Name | Default | What it does |
|---|---|---|
| `TE_CLASSES`, `SATELLITE_CLASSES` | 26, 54 | the class lists; `split_repeats.py` validates against them |
| `REPEAT_WINDOW` | `HMM_WINDOW` (1000) | the grid the annotations are clipped to |
| `PHASE5_SUBSAMPLES` | `"1000"` | fragments per iteration, and the significance-test threshold |
| `PHASE5_ITERATIONS` | `"100"` | iterations. The methods call 100 a ceiling, not a budget -- do not raise it |

To analyse the merged PDAL-Seq dataset as well as the individual libraries, add
a rule that runs Phase 1's `bam_to_bedgraph` / `bedgraph_to_bigwig` on
`data/phase2/All_reads/<experiment>.bam` and add the experiment names to the
`sample=` list in `compile_repeat_enrichment`. It was left out because five
libraries analysed separately answer the question and the merge would only
average them.

## Run it

```
rsync -Pav ./ submit.hpc.psu.edu:/storage/group/kdm16/default/jus841/js4-2026/js4025_Analyze_Zebra_finch_PDAL-Seq_data

module load anaconda
conda activate js4025_snakemake

# Should report 588 new jobs and nothing from Phase 1, 2, 3 or 4.
snakemake --snakefile workflow/Snakefile -n --cores 20 --use-conda

mkdir -p slurm_logs && sbatch submit_phase5.sh
```

588 jobs: 3 splits, 1 report, 83 window clippings, 498 enrichments, 1 compile,
2 plots. If the count comes back as 590, `chrom_sizes` and
`make_genome_windows` are in there too -- Phase 3 built the 1 kbp grid, so on
Roar they should already exist.

`submit_phase5.sh` uses the same header as Phases 3 and 4 -- `basic`, 1 node,
20 tasks, 4 GB/core, no `--account` line -- with 24 h rather than 48 h. The
phase is ~2-3 h of work at 20-way parallelism and nothing in it re-reads a BAM.
The script names its three targets rather than building `rule all`, so it
cannot restart a Phase 3 HMM fit or 15 h of Phase 4 STAR by accident.

Restarting is safe: `Read_enrichment.sh` stages its output and renames it
atomically, so a job killed at the wall clock leaves nothing that could be
mistaken for a finished result. Resubmit until `snakemake -n` reports nothing
to do.

## Sanity checks after the run

```
# 1. Did every class survive with a sensible number of records? Compare against
#    the counts in resources/README.md.
column -t results/phase5/Repeat_annotation_report.txt | head -40

# 2. Did any class end up resting on a handful of fragments? These are the rows
#    whose enrichment is one draw, not a mean of 100.
awk -F, 'NR>1 && $11 < 1000 {print $8, $11}' results/phase5/Repeat_enrichment.csv | sort -u

# 3. Is the PCR-free control itself enriched anywhere? Those are the mapping
#    artefacts the ratio exists to divide out.
awk -F, 'NR==1 || $2=="WGS"' results/phase5/Repeat_enrichment.csv \
    | sort -t, -k15,15gr | head -15

# 4. Where does the 0 mM negative control sit?
awk -F, 'NR==1 || $2=="0mM"' results/phase5/PDALSeq_versus_control.csv \
    | sort -t, -k12,12gr | head -15

# 5. Did any enrichment job leave a temp directory behind (i.e. was killed)?
ls -d temp/read_enrichment_* 2>/dev/null | wc -l
```

## Decisions recorded for this phase

- **No blacklist, no MAPQ filter, no chromosome exclusion.** The repeats are
  the subject; every filter in this pipeline removes them. This is the same
  reasoning as `phase2_including_repetative`.
- **The coverage source is Phase 1's bigWig, not the BAM.** Same statistic,
  reached by index instead of by ~100,000 whole-BAM scans. This is the only
  substantive change to the manuscript's script.
- **Deduplicated reads (`All_reads`), as the manuscript used.** `samtools
  rmdup -s` collapses single-end reads sharing a start coordinate, which in a
  deep satellite array will remove some real fragments -- but it removes them
  from the PCR-free control identically, and the control is the point.
- **The window intersection is computed once**, not inside every iteration.
- **`bedtools intersect` keeps its default overlap output and its duplicates.**
  Annotation-weighted sampling is what the methods describe; `-u` or a prior
  `bedtools merge` would quietly change it to window-uniform.
- **The Wilcoxon threshold follows the subsample size (1,000), not the literal
  10,000.** In the manuscript those were the same number, and the subsampling
  is what the rule is about.
- **Two pooled classes (`all_TE`, `all_Satellite`) were added.** The manuscript
  had no equivalent; 80 families each resting on a few loci need a line to be
  read against.
- **Class lists are hardcoded and validated.** `split_repeats.py` fails with
  the offending names rather than letting a renamed class surface as a missing
  output.
- **`split_repeats.py` is Python, not awk.** It writes up to 55 files from one
  input, and awk drops output files under parallel jobs.
- **One new environment file** (`phase5_repeats.yml`). No existing yml was
  edited -- snakemake hashes them, and an edit would re-queue every finished
  Phase 1-4 job that used one.
