####################################################################
# Phase 2 -- Merge well-correlated datasets and compare PDAL-Seq /
#            control read density to non-B DNA motif density.
#
# Adapted from manuscript steps 03 (blacklist), 04 (non-B comparison)
# and 21 (the K562 per-experiment merge pattern).
#
# Phase 1 established that the PDAL-Seq replicates correlate tightly and
# the WGS control does not resemble them, so 20mM-1/20mM-2 and
# 40mM-1/40mM-2 merge into concentration-level BAMs and all four merge
# into one combined PDAL-Seq dataset.
####################################################################


####################################################################
# Per-experiment BAM merge
####################################################################

def experiment_bams(wildcards):
    """Phase 1 BAMs belonging to this experiment, in sample-sheet order."""
    return expand(
        "data/phase1/{mapq}/{sample}.bam",
        mapq=wildcards.mapq,
        sample=EXPERIMENTS[wildcards.experiment],
    )


rule merge_experiment_bam:
    # Runs for single-library experiments too (WGS, 0mM), where samtools
    # merge is effectively a copy. Uniform paths downstream are worth the
    # duplicated bytes -- every later rule has one input pattern.
    input:
        experiment_bams
    output:
        "data/phase2/{mapq}/{experiment}.bam"
    log:
        "logs/phase2/merge_experiment_bam/{mapq}/{experiment}.txt"
    benchmark:
        "benchmarks/phase2/merge_experiment_bam/{mapq}/{experiment}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        "samtools merge -f -@ {threads} -o {output} {input} 2> {log}"


rule index_experiment_bam:
    input:
        "data/phase2/{mapq}/{experiment}.bam"
    output:
        "data/phase2/{mapq}/{experiment}.bam.bai"
    log:
        "logs/phase2/index_experiment_bam/{mapq}/{experiment}.txt"
    benchmark:
        "benchmarks/phase2/index_experiment_bam/{mapq}/{experiment}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


####################################################################
# Non-B DNA motif annotations
#
# The Zenodo release annotates the FULL diploid assembly, so every
# autosome appears twice (chr10_mat and chr10_pat). Our genome is the
# extracted single haplotype, so roughly half of each file is dropped.
####################################################################

rule clean_nonb:
    wildcard_constraints:
        nonb = "|".join(NONB)
    input:
        bed         = ancient(NONB_PREFIX + ".{nonb}.bed"),
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        bed    = "results/phase2/nonb_annotations/{nonb}.bed",
        report = "data/phase2/nonb_annotations/{nonb}.report.txt",
    log:
        "logs/phase2/clean_nonb/{nonb}.txt"
    benchmark:
        "benchmarks/phase2/clean_nonb/{nonb}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase2_merge_and_nonb/clean_nonb.sh "
        "-i {input.bed} -g {input.chrom_sizes} -o {output.bed} "
        "-r {output.report} -n {wildcards.nonb} -d temp 2>&1 | tee {log}"


rule combine_nonb:
    input:
        beds        = expand("results/phase2/nonb_annotations/{nonb}.bed", nonb=NONB),
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "results/phase2/nonb_annotations/all.bed"
    log:
        "logs/phase2/combine_nonb/all.txt"
    benchmark:
        "benchmarks/phase2/combine_nonb/all.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 8000
    shell:
        # Concatenated, NOT merged -- this matches the manuscript's all.sh.
        # Classes overlap each other, so "all" double counts by design; the
        # base_density metric in the compile step is the overlap-safe one.
        "bash workflow/scripts/phase2_merge_and_nonb/combine_nonb.sh "
        "-g {input.chrom_sizes} -o {output} -d temp {input.beds} 2>&1 | tee {log}"


rule nonb_annotation_report:
    input:
        expand("data/phase2/nonb_annotations/{nonb}.report.txt", nonb=NONB)
    output:
        "results/phase2/nonB_annotation_report.txt"
    log:
        "logs/phase2/nonb_annotation_report.txt"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 1000
    shell:
        "( printf 'Class\\tInput_records\\tRetained_records\\tPercent_retained\\tContigs_retained\\n' ; "
        "cat {input} ) > {output} 2> {log}"


rule nonb_density:
    input:
        windows     = "data/phase1/genome_windows/{size}_nucleotides/" + GENOME + ".bed",
        nonb        = "results/phase2/nonb_annotations/{nonb}.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase2/nonb_density/{size}_nucleotides/{nonb}.bg"
    log:
        "logs/phase2/nonb_density/{size}_nucleotides/{nonb}.txt"
    benchmark:
        "benchmarks/phase2/nonb_density/{size}_nucleotides/{nonb}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 8000
    shell:
        # -sorted streams instead of loading -b into memory, which matters
        # for IR (4.6 M intervals). It requires both files in the -g order,
        # which is .fai order -- clean_nonb.sh guarantees that.
        # Output: chr start end count bases length base_density
        "bedtools coverage -a {input.windows} -b {input.nonb} "
        "-sorted -g {input.chrom_sizes} > {output} 2> {log}"


####################################################################
# Blacklist of unmappable regions (deferred here from Phase 1)
#
# The window coverage this needs was already written by Phase 1, so
# nothing is recomputed.
####################################################################

rule make_raw_blacklist:
    input:
        all_reads = "data/phase1/genome_window_coverage/All_reads/"
                    + BLACKLIST_WINDOW + "_nucleotides/" + WGS_SAMPLE + ".bed",
        high_mapq = "data/phase1/genome_window_coverage/High_MapQ/"
                    + BLACKLIST_WINDOW + "_nucleotides/" + WGS_SAMPLE + ".bed",
    output:
        "data/phase2/blacklist/raw.bed"
    params:
        threshold = BLACKLIST_THRESHOLD
    log:
        "logs/phase2/make_raw_blacklist.txt"
    benchmark:
        "benchmarks/phase2/make_raw_blacklist/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        "Rscript workflow/scripts/phase2_merge_and_nonb/Black_list_windows.R "
        "{input.all_reads} {input.high_mapq} {params.threshold} {output} 2>&1 | tee {log}"


rule merge_blacklist:
    input:
        "data/phase2/blacklist/raw.bed"
    output:
        "data/phase2/blacklist/merged.bed"
    log:
        "logs/phase2/merge_blacklist.txt"
    benchmark:
        "benchmarks/phase2/merge_blacklist/all.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 4000
    shell:
        # bedtools merge wants lexicographic chrom order, which is NOT the
        # .fai order the raw file is written in. Sort into it here and let
        # add_sex_chromosomes put the file back into .fai order.
        "LC_ALL=C sort -k1,1 -k2,2n -S 2G -T temp {input} "
        "| bedtools merge -d 1 -i - > {output} 2> {log}"


rule add_sex_chromosomes_to_blacklist:
    input:
        blacklist   = "data/phase2/blacklist/merged.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "results/phase2/Blacklist.bed"
    params:
        sex = ",".join(SEX_CHROMOSOMES)
    log:
        "logs/phase2/add_sex_chromosomes_to_blacklist.txt"
    benchmark:
        "benchmarks/phase2/add_sex_chromosomes_to_blacklist/all.tsv"
    conda:
        "../env/phase2_nonb.yml"
    shell:
        "bash workflow/scripts/phase2_merge_and_nonb/add_sex_chromosomes_to_blacklist.sh "
        "-b {input.blacklist} -g {input.chrom_sizes} -s {params.sex} "
        "-o {output} -d temp 2>&1 | tee {log}"


rule blacklist_window_flag:
    input:
        windows   = "data/phase1/genome_windows/{size}_nucleotides/" + GENOME + ".bed",
        blacklist = "results/phase2/Blacklist.bed",
    output:
        "data/phase2/blacklist/window_flag/{size}_nucleotides.bed"
    log:
        "logs/phase2/blacklist_window_flag/{size}_nucleotides.txt"
    benchmark:
        "benchmarks/phase2/blacklist_window_flag/{size}_nucleotides.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 4000
    shell:
        # -c appends the number of blacklist intervals hitting each window.
        # Any window with a non-zero count is dropped in the filtered branch.
        "bedtools intersect -a {input.windows} -b {input.blacklist} -c "
        "> {output} 2> {log}"


####################################################################
# Read density in windows, per merged experiment
####################################################################

rule experiment_window_coverage:
    input:
        bam     = "data/phase2/{mapq}/{experiment}.bam",
        windows = "data/phase1/genome_windows/{size}_nucleotides/" + GENOME + ".bed",
    output:
        "data/phase2/experiment_window_coverage/{mapq}/{size}_nucleotides/{experiment}.bed"
    log:
        "logs/phase2/experiment_window_coverage/{mapq}/{size}_nucleotides/{experiment}.txt"
    benchmark:
        "benchmarks/phase2/experiment_window_coverage/{mapq}/{size}_nucleotides/{experiment}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 8000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum "
        "> {output} 2> {log}"


####################################################################
# Compile and plot
####################################################################

rule compile_reads_vs_nonb:
    input:
        coverage = expand(
            "data/phase2/experiment_window_coverage/All_reads/{size}_nucleotides/{experiment}.bed",
            size=WINDOWS, experiment=EXPERIMENT_NAMES,
        ),
        density = expand(
            "data/phase2/nonb_density/{size}_nucleotides/{nonb}.bg",
            size=WINDOWS, nonb=NONB_ALL,
        ),
        flags = expand(
            "data/phase2/blacklist/window_flag/{size}_nucleotides.bed",
            size=WINDOWS,
        ),
        samples = "Zebrafinch_samples.txt",
    output:
        correlation = "results/phase2/nonB_correlation_coefficient.csv",
        compare     = "results/phase2/nonB_coverage_versus_read.csv",
    params:
        experiments    = ",".join(EXPERIMENT_NAMES),
        nonb           = ",".join(NONB_ALL),
        windows        = ",".join(WINDOWS),
        wgs            = WGS_EXPERIMENT,
        exclude        = ",".join(EXCLUDE_CHROMOSOMES),
        compare_window = max(WINDOWS, key=int),
    log:
        "logs/phase2/compile_reads_vs_nonb.txt"
    benchmark:
        "benchmarks/phase2/compile_reads_vs_nonb/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 16000
    shell:
        "Rscript workflow/scripts/phase2_merge_and_nonb/Compile_reads_vs_nonb.R "
        "{params.experiments} {params.nonb} {params.windows} {params.wgs} "
        "{params.exclude} {params.compare_window} "
        "{output.correlation} {output.compare} 2>&1 | tee {log}"


rule plot_nonb_correlation:
    input:
        "results/phase2/nonB_correlation_coefficient.csv"
    output:
        "plots/phase2/nonB_correlation_heatmap.svg"
    log:
        "logs/phase2/plot_nonb_correlation.txt"
    benchmark:
        "benchmarks/phase2/plot_nonb_correlation/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase2_merge_and_nonb/Plot_nonb_correlation.R "
        "{input} {output} 2>&1 | tee {log}"


rule plot_nonb_scatter:
    input:
        compare = "results/phase2/nonB_coverage_versus_read.csv",
        samples = "Zebrafinch_samples.txt",
    output:
        "plots/phase2/nonB_scatter.svg"
    params:
        wgs      = WGS_EXPERIMENT,
        combined = COMBINED_EXPERIMENT,
    log:
        "logs/phase2/plot_nonb_scatter.txt"
    benchmark:
        "benchmarks/phase2/plot_nonb_scatter/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase2_merge_and_nonb/Plot_nonb_scatter.R "
        "{input.compare} {params.wgs} {params.combined} {output} 2>&1 | tee {log}"
