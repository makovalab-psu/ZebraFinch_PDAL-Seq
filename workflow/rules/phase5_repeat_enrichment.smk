####################################################################
# Phase 5 -- PDAL-Seq versus control read density in satellites,
#            centromeres and transposable elements.
#
# Adapted from the PDAL-Seq manuscript pipeline:
#   rules/17_PDALseq_reads_in_CenSat.smk
#   scripts/17_PDALseq_reads_in_CenSat/Read_enrichment.sh
#   scripts/17_PDALseq_reads_in_CenSat/Compile_CenSat_Class_enrichment.R
#
# This phase is the opposite of Phases 3 and 4 in what it filters. There is no
# blacklist, no MAPQ filter and no chromosome exclusion: the subject is the
# repetitive genome, and every filter this pipeline has removes exactly that.
# The coverage source is the All_reads bigWig Phase 1 already built from the
# deduplicated, unfiltered BAM.
#
# Nothing is remapped and no BAM is read. Phase 1 wrote
# data/phase1/bigwig/<sample>.bigwig from `bedtools genomecov -bg -split -ibam`
# on the All_reads BAM, so the per-base coverage this phase sums is that BAM's
# coverage, reached by index instead of by a full scan. See the header of
# scripts/phase5_repeat_enrichment/Read_enrichment.sh.
#
# The comparison that controls for mapping artefacts is the PCR-free library
# run through the identical procedure: a repeat class that collects reads
# because bwa mem scattered multi-mappers into it collects them in BOTH
# libraries, and dividing removes it. Compile_repeat_enrichment.R does that
# division and writes results/phase5/PDALSeq_versus_control.csv.
####################################################################


####################################################################
# Split the repeat annotations into one BED per class
####################################################################

rule split_te_classes:
    input:
        bed         = ancient(REPEAT_TE_BED),
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        beds = expand(REPEAT_CLASS_DIR + "/{repeat_class}.bed",
                      repeat_class=TE_CLASSES + [TE_ALL]),
        report = "data/phase5/reports/TE.txt",
    params:
        outdir   = REPEAT_CLASS_DIR,
        column   = REPEAT_TE_COLUMN,
        all_name = TE_ALL,
        expected = ",".join(TE_CLASSES),
    log:
        "logs/phase5/split_te_classes.txt"
    benchmark:
        "benchmarks/phase5/split_te_classes/all.tsv"
    conda:
        "../env/phase5_repeats.yml"
    resources:
        mem_mb = 8000
    shell:
        "python3 workflow/scripts/phase5_repeat_enrichment/split_repeats.py "
        "--input {input.bed} --chrom-sizes {input.chrom_sizes} "
        "--outdir {params.outdir} --report {output.report} "
        "--source TE --column {params.column} --all-name {params.all_name} "
        "--expected {params.expected} 2>&1 | tee {log}"


rule split_satellite_classes:
    input:
        bed         = ancient(REPEAT_SATELLITE_BED),
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        beds = expand(REPEAT_CLASS_DIR + "/{repeat_class}.bed",
                      repeat_class=SATELLITE_CLASSES + [SATELLITE_ALL]),
        report = "data/phase5/reports/Satellite.txt",
    params:
        outdir   = REPEAT_CLASS_DIR,
        column   = REPEAT_SATELLITE_COLUMN,
        all_name = SATELLITE_ALL,
        expected = ",".join(SATELLITE_CLASSES),
    log:
        "logs/phase5/split_satellite_classes.txt"
    benchmark:
        "benchmarks/phase5/split_satellite_classes/all.tsv"
    conda:
        "../env/phase5_repeats.yml"
    resources:
        mem_mb = 8000
    shell:
        "python3 workflow/scripts/phase5_repeat_enrichment/split_repeats.py "
        "--input {input.bed} --chrom-sizes {input.chrom_sizes} "
        "--outdir {params.outdir} --report {output.report} "
        "--source Satellite --column {params.column} "
        "--all-name {params.all_name} "
        "--expected {params.expected} 2>&1 | tee {log}"


# The centromere BED has no class column and one record per chromosome, so it
# goes through the same script with --fixed-name. It is cleaned here rather
# than reused from results/phase4/annotations/CEN.bed on purpose: Phase 4's
# copy was clamped against a chrom.sizes it shares with this phase, but making
# Phase 5 depend on a Phase 4 output would put a finished 40 h phase inside
# this phase's DAG.
rule split_centromere:
    input:
        bed         = ancient(CENTROMERE_BED),
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        beds   = expand(REPEAT_CLASS_DIR + "/{repeat_class}.bed",
                        repeat_class=[CENTROMERE_CLASS]),
        report = "data/phase5/reports/Centromere.txt",
    params:
        outdir     = REPEAT_CLASS_DIR,
        fixed_name = CENTROMERE_CLASS,
        expected   = CENTROMERE_CLASS,
    log:
        "logs/phase5/split_centromere.txt"
    benchmark:
        "benchmarks/phase5/split_centromere/all.tsv"
    conda:
        "../env/phase5_repeats.yml"
    resources:
        mem_mb = 2000
    shell:
        "python3 workflow/scripts/phase5_repeat_enrichment/split_repeats.py "
        "--input {input.bed} --chrom-sizes {input.chrom_sizes} "
        "--outdir {params.outdir} --report {output.report} "
        "--source Centromere --fixed-name {params.fixed_name} "
        "--expected {params.expected} 2>&1 | tee {log}"


rule repeat_annotation_report:
    input:
        expand("data/phase5/reports/{source}.txt",
               source=["TE", "Satellite", "Centromere"])
    output:
        "results/phase5/Repeat_annotation_report.txt"
    log:
        "logs/phase5/repeat_annotation_report.txt"
    benchmark:
        "benchmarks/phase5/repeat_annotation_report/all.tsv"
    resources:
        mem_mb = 1000
    shell:
        # Read this before reading any enrichment: a class with three records
        # gives an enrichment computed from three records, however dramatic
        # the number looks.
        "( printf 'Source\\tClass\\tRecords\\tContigs\\tBp\\n'; cat {input} ) "
        "> {output} 2> {log}"


####################################################################
# Clip each class to the 1 kbp window grid
#
# This is the manuscript's `bedtools intersect -a windows_bed -b annotation`,
# hoisted out of the iteration loop: the manuscript recomputed it inside every
# one of its 100 iterations, which produces the same file 100 times.
#
# The intersection is what gets subsampled, so it is worth being clear about
# what the units are. `bedtools intersect` without -wa reports the OVERLAP,
# so a fragment is "the part of one annotation that falls in one 1 kbp
# window" -- at most 1 kbp long, and one fragment per (annotation, window)
# pair. That is the manuscript's behaviour and it is kept:
#
#   * the 1 kbp cap bounds the length of the sampled intervals, which is what
#     makes `bedtools shuffle` able to place them and what keeps one 300 kbp
#     satellite array from being a single sampling unit;
#
#   * a window covered by many elements of a class yields many fragments, so
#     the subsample is weighted by annotation density -- the methods' "each
#     annotation type was first subsampled at random", sampling annotations
#     rather than windows. Merging the class first, or `intersect -u`, would
#     silently change that to uniform sampling over windows.
####################################################################

rule repeat_class_fragments:
    input:
        class_bed   = REPEAT_CLASS_DIR + "/{repeat_class}.bed",
        windows     = "data/phase1/genome_windows/" + REPEAT_WINDOW
                      + "_nucleotides/" + GENOME + ".bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase5/repeat_fragments/{repeat_class}.bed"
    log:
        "logs/phase5/repeat_class_fragments/{repeat_class}.txt"
    benchmark:
        "benchmarks/phase5/repeat_class_fragments/{repeat_class}.tsv"
    conda:
        "../env/phase5_repeats.yml"
    resources:
        # -sorted streams both files instead of loading -b into memory. Both
        # are already in chrom.sizes order -- makewindows emits that order and
        # split_repeats.py sorts by .fai rank -- so -g is a check, not a sort.
        mem_mb = 4000
    shell:
        "bedtools intersect -a {input.windows} -b {input.class_bed} "
        "-sorted -g {input.chrom_sizes} > {output} 2> {log}"


####################################################################
# Enrichment, one job per sample x repeat class
####################################################################

rule repeat_read_enrichment:
    input:
        coverage    = "data/phase1/bigwig/{sample}.bigwig",
        fragments   = "data/phase5/repeat_fragments/{repeat_class}.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase5/repeat_enrichment/{sample}/{repeat_class}.csv"
    params:
        subsamples = PHASE5_SUBSAMPLES,
        iterations = PHASE5_ITERATIONS,
    log:
        "logs/phase5/repeat_read_enrichment/{sample}/{repeat_class}.txt"
    benchmark:
        "benchmarks/phase5/repeat_read_enrichment/{sample}/{repeat_class}.tsv"
    conda:
        "../env/phase5_repeats.yml"
    resources:
        # bedtools holds only the 1,000-record subsample; megadepth holds the
        # bigWig index and one range at a time.
        mem_mb = 2000
    shell:
        "bash workflow/scripts/phase5_repeat_enrichment/Read_enrichment.sh "
        "-c {input.coverage} -a {input.fragments} -g {input.chrom_sizes} "
        "-o {output} -n {params.subsamples} -i {params.iterations} "
        "-d temp 2>&1 | tee {log}"


####################################################################
# Compile and plot
####################################################################

rule compile_repeat_enrichment:
    input:
        enrichment = expand(
            "data/phase5/repeat_enrichment/{sample}/{repeat_class}.csv",
            sample=SAMPLES, repeat_class=REPEAT_CLASSES,
        ),
        sample_table = "Zebrafinch_samples.txt",
    output:
        enrichment = "results/phase5/Repeat_enrichment.csv",
        versus     = "results/phase5/PDALSeq_versus_control.csv",
    params:
        in_dir     = "data/phase5/repeat_enrichment",
        te         = ",".join(TE_CLASSES),
        satellite  = ",".join(SATELLITE_CLASSES),
        other      = ",".join([TE_ALL, SATELLITE_ALL, CENTROMERE_CLASS]),
        subsamples = PHASE5_SUBSAMPLES,
        control    = WGS_SAMPLE,
    log:
        "logs/phase5/compile_repeat_enrichment.txt"
    benchmark:
        "benchmarks/phase5/compile_repeat_enrichment/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        "Rscript workflow/scripts/phase5_repeat_enrichment/Compile_repeat_enrichment.R "
        "{params.in_dir} {input.sample_table} {params.te} {params.satellite} "
        "{params.other} {params.subsamples} {params.control} "
        "{output.enrichment} {output.versus} 2>&1 | tee {log}"


rule plot_repeat_enrichment:
    input:
        enrichment = "results/phase5/Repeat_enrichment.csv",
        versus     = "results/phase5/PDALSeq_versus_control.csv",
    output:
        "plots/phase5/Repeat_enrichment_heatmap.svg"
    log:
        "logs/phase5/plot_repeat_enrichment.txt"
    benchmark:
        "benchmarks/phase5/plot_repeat_enrichment/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase5_repeat_enrichment/Plot_repeat_enrichment.R "
        "{input.enrichment} {input.versus} {output} 2>&1 | tee {log}"


rule plot_repeat_versus_control:
    input:
        "results/phase5/PDALSeq_versus_control.csv"
    output:
        "plots/phase5/PDALSeq_versus_control.svg"
    log:
        "logs/phase5/plot_repeat_versus_control.txt"
    benchmark:
        "benchmarks/phase5/plot_repeat_versus_control/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase5_repeat_enrichment/Plot_versus_control.R "
        "{input} {output} 2>&1 | tee {log}"
