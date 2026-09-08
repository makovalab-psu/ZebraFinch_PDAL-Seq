####################################################################
# Phase 2 including repetitive -- the same read-density versus non-B
# DNA motif density comparison as Phase 2, on EVERY window of the
# genome.
#
# Phase 2 draws its figures from windows that survive two filters: the
# mappability blacklist and the EXCLUDE_CHROMOSOMES contig list. Regions
# land on the blacklist BECAUSE they are repetitive, and bwa mem assigns
# multi-mapping reads to one of the alternatives at random rather than
# discarding them -- so the reads are there, the motifs are there, and
# the filter removes exactly the repeat content worth asking about.
#
# The manuscript agrees: 04_.../compile_reads_vs_nonb.R never opens
# blacklist_files/. Its only filter is !chr %in% c("chrX","chrY"). The
# blacklist first appears at step 09, the HMM, which is where Phase 3
# applies it here.
#
# Nothing expensive is recomputed. The blacklist filter in Phase 2 lives
# entirely inside its R scripts, so the coverage and motif-density files
# underneath already cover the whole genome. This phase adds one small
# bedtools job per window size and four R jobs.
#
# Phase 2 itself is NOT modified. Both results stand side by side.
####################################################################

####################################################################
# How much of each window is blacklisted?
#
# Phase 2's blacklist_window_flag is a boolean: any overlap at all and
# the window is dropped. That is defensible at 10 kbp and misleading at
# 1 Mbp, because merge_blacklist runs `bedtools merge -d 1` and the
# resulting intervals are long. Measured on the manuscript's own human
# blacklist -- 9500 intervals, mean 65 kbp, 14% of autosomal sequence --
# the boolean drops 93.0% of 1 Mbp windows, 37.0% of 100 kbp windows and
# 13.9% of 10 kbp windows. Phase 2's 1 Mbp scatter is therefore fitting
# the cleanest 7% of the genome.
#
# A fraction degrades gracefully where a boolean does not, so this rule
# measures bases rather than counting overlaps.
####################################################################

rule repeat_blacklist_window_fraction:
    input:
        windows     = "data/phase1/genome_windows/{size}_nucleotides/" + GENOME + ".bed",
        blacklist   = "results/phase2/Blacklist.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase2_including_repetative/blacklist_fraction/{size}_nucleotides.bed"
    log:
        "logs/phase2_including_repetative/blacklist_window_fraction/{size}_nucleotides.txt"
    benchmark:
        "benchmarks/phase2_including_repetative/blacklist_window_fraction/{size}_nucleotides.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 4000
    shell:
        # Same invocation as nonb_density. -sorted is safe here: the windows
        # come out of make_windows.sh in .fai order and Blacklist.bed is put
        # back into .fai order by the rank sort at the end of
        # add_sex_chromosomes_to_blacklist.sh.
        # Output: chr start end count bases length fraction
        "bedtools coverage -a {input.windows} -b {input.blacklist} "
        "-sorted -g {input.chrom_sizes} > {output} 2> {log}"


####################################################################
# Compile
####################################################################

rule compile_reads_vs_nonb_all_regions:
    input:
        coverage = expand(
            "data/phase2/experiment_window_coverage/All_reads/{size}_nucleotides/{experiment}.bed",
            size=WINDOWS, experiment=EXPERIMENT_NAMES,
        ),
        density = expand(
            "data/phase2/nonb_density/{size}_nucleotides/{nonb}.bg",
            size=WINDOWS, nonb=NONB_ALL,
        ),
        fraction = expand(
            "data/phase2_including_repetative/blacklist_fraction/{size}_nucleotides.bed",
            size=WINDOWS,
        ),
        samples = "Zebrafinch_samples.txt",
    output:
        correlation = "results/phase2_including_repetative/nonB_correlation_coefficient.csv",
        compare     = "results/phase2_including_repetative/nonB_coverage_versus_read.csv",
        composition = "results/phase2_including_repetative/Region_composition.csv",
    params:
        experiments    = ",".join(EXPERIMENT_NAMES),
        nonb           = ",".join(NONB_ALL),
        windows        = ",".join(WINDOWS),
        wgs            = WGS_EXPERIMENT,
        exclude        = ",".join(EXCLUDE_CHROMOSOMES),
        compare_window = max(WINDOWS, key=int),
        light          = REPEAT_LIGHT_FRACTION,
        heavy          = REPEAT_HEAVY_FRACTION,
    log:
        "logs/phase2_including_repetative/compile_reads_vs_nonb_all_regions.txt"
    benchmark:
        "benchmarks/phase2_including_repetative/compile_reads_vs_nonb_all_regions/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 16000
    shell:
        "Rscript workflow/scripts/phase2_including_repetative/Compile_reads_vs_nonb_all_regions.R "
        "{params.experiments} {params.nonb} {params.windows} {params.wgs} "
        "{params.exclude} {params.compare_window} {params.light} {params.heavy} "
        "{output.correlation} {output.compare} {output.composition} 2>&1 | tee {log}"


####################################################################
# Plot
####################################################################

rule plot_nonb_correlation_all_regions:
    input:
        "results/phase2_including_repetative/nonB_correlation_coefficient.csv"
    output:
        "plots/phase2_including_repetative/nonB_correlation_heatmap.svg"
    log:
        "logs/phase2_including_repetative/plot_nonb_correlation_all_regions.txt"
    benchmark:
        "benchmarks/phase2_including_repetative/plot_nonb_correlation_all_regions/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase2_including_repetative/Plot_nonb_correlation_all_regions.R "
        "{input} {output} 2>&1 | tee {log}"


rule plot_nonb_correlation_delta:
    input:
        "results/phase2_including_repetative/nonB_correlation_coefficient.csv"
    output:
        "plots/phase2_including_repetative/nonB_correlation_delta.svg"
    params:
        light = REPEAT_LIGHT_FRACTION
    log:
        "logs/phase2_including_repetative/plot_nonb_correlation_delta.txt"
    benchmark:
        "benchmarks/phase2_including_repetative/plot_nonb_correlation_delta/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase2_including_repetative/Plot_nonb_correlation_delta.R "
        "{input} {params.light} {output} 2>&1 | tee {log}"


rule plot_nonb_scatter_all_regions:
    input:
        compare = "results/phase2_including_repetative/nonB_coverage_versus_read.csv",
        samples = "Zebrafinch_samples.txt",
    output:
        "plots/phase2_including_repetative/nonB_scatter.svg"
    params:
        wgs      = WGS_EXPERIMENT,
        combined = COMBINED_EXPERIMENT,
        light    = REPEAT_LIGHT_FRACTION,
    log:
        "logs/phase2_including_repetative/plot_nonb_scatter_all_regions.txt"
    benchmark:
        "benchmarks/phase2_including_repetative/plot_nonb_scatter_all_regions/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase2_including_repetative/Plot_nonb_scatter_all_regions.R "
        "{input.compare} {params.wgs} {params.combined} {params.light} {output} 2>&1 | tee {log}"
