####################################################################
# Phase 1 -- Data preprocessing and comparison in genome windows
#
# Adapted from the PDAL-Seq manuscript pipeline:
#   rules/01_preprocess_data.smk                       (trim, map, dedup, MapQ, bigwig)
#   rules/02_preprocessing_statistics.smk              (read counts, bam flags, AS-XS)
#   rules/03_Correlation_and_black_list_unmappable_regions.smk  (windows, megadepth, correlation)
#   rules/21_add_K562_to_analysis.smk                  (self-contained single-genome variant)
#
# Differences from the manuscript, all deliberate:
#   * R1 only. There is no R2 rule anywhere.
#   * One genome, given by GENOME in the Snakefile, so nothing infers a
#     species from a file name the way map_reads.sh and
#     genome_window_coverage.sh do in the manuscript.
#   * bwa mem streams straight into samtools sort instead of writing a
#     temp SAM. A 119 M read SAM is ~40 GB on scratch for no benefit.
#   * chrom.sizes is cut from the existing .fai rather than built with
#     faSize, which drops the ucsc-fasize dependency.
#
# Blacklist generation and per-experiment BAM merging are NOT here. Both
# depend on knowing which datasets correlate, which is what this phase
# answers, so they belong to Phase 2.
####################################################################

# Scripts live in workflow/scripts/phase1_preprocess_and_correlate/ and are
# referenced by literal path in the shell blocks below (a global Python
# name is not reliably visible to Snakemake string formatting).

# BAM stages that get flag counts. Mapped_bam is pre-deduplication.
BAM_STAGES = ["Mapped_bam", "All_reads", "High_MapQ"]


####################################################################
# Genome preparation
####################################################################

rule bwa_index:
    input:
        GENOME_FASTA
    output:
        amb = GENOME_FASTA + ".amb",
        ann = GENOME_FASTA + ".ann",
        pac = GENOME_FASTA + ".pac",
        bwt = GENOME_FASTA + ".bwt",
        sa  = GENOME_FASTA + ".sa",
    log:
        "logs/phase1/bwa_index/" + GENOME + ".txt"
    benchmark:
        "benchmarks/phase1/bwa_index/" + GENOME + ".tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 8000
    shell:
        "bwa index {input} 2>&1 | tee {log}"


rule chrom_sizes:
    input:
        # Built by resources/make_PDAL-Seq_fasta.sh. ancient() so a stray
        # touch on the index does not rebuild the sizes file. The wildcard
        # is carried through so this rule stays correct if a second genome
        # is ever added.
        ancient("resources/genomes/{genome}.fa.fai")
    output:
        "resources/genomes/{genome}.chrom.sizes"
    log:
        "logs/phase1/chrom_sizes/{genome}.txt"
    benchmark:
        "benchmarks/phase1/chrom_sizes/{genome}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 1000
    shell:
        "cut -f 1,2 {input} > {output} 2>&1 | tee {log}"


####################################################################
# Per-sample preprocessing
####################################################################

rule trim_fastq:
    input:
        FASTQ_DIR + "/{sample}.fastq.gz"
    output:
        fastq = temp("data/phase1/Trimmed_fastq/{sample}.fastq.gz"),
        json  = "results/phase1/fastp_reports/{sample}.json",
        html  = "results/phase1/fastp_reports/{sample}.html",
    log:
        "logs/phase1/trim_fastq/{sample}.txt"
    benchmark:
        "benchmarks/phase1/trim_fastq/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 4000
    threads:
        4
    shell:
        # Flags are identical to rule trim_fastq in the manuscript's
        # 01_preprocess_data.smk: overrepresentation analysis, polyG and
        # polyX trimming, 3' quality trimming.
        "fastp --thread {threads} -i {input} -o {output.fastq} "
        "-p -P 50 -y -x -3 "
        "--json {output.json} --html {output.html} 2>&1 | tee {log}"


rule map_reads:
    input:
        fastq  = "data/phase1/Trimmed_fastq/{sample}.fastq.gz",
        genome = GENOME_FASTA,
        amb    = GENOME_FASTA + ".amb",
        ann    = GENOME_FASTA + ".ann",
        pac    = GENOME_FASTA + ".pac",
        bwt    = GENOME_FASTA + ".bwt",
        sa     = GENOME_FASTA + ".sa",
    output:
        temp("data/phase1/Mapped_bam/{sample}.bam")
    log:
        "logs/phase1/map_reads/{sample}.txt"
    benchmark:
        "benchmarks/phase1/map_reads/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 16000
    threads:
        10
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/map_reads.sh -f {input.fastq} -g {input.genome} "
        "-o {output} -t {threads} -d temp 2>&1 | tee {log}"


rule remove_dup:
    input:
        "data/phase1/Mapped_bam/{sample}.bam"
    output:
        "data/phase1/All_reads/{sample}.bam"
    log:
        "logs/phase1/remove_dup/{sample}.txt"
    benchmark:
        "benchmarks/phase1/remove_dup/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 8000
    shell:
        # -s = single-end mode. Applied to the PCR-free control as well,
        # so the Phase 2 denominator is processed exactly like the
        # PDAL-Seq numerator (see README, Phase 1 decisions).
        "samtools rmdup -s {input} {output} 2>&1 | tee {log}"


rule remove_low_mapq:
    input:
        "data/phase1/All_reads/{sample}.bam"
    output:
        "data/phase1/High_MapQ/{sample}.bam"
    log:
        "logs/phase1/remove_low_mapq/{sample}.txt"
    benchmark:
        "benchmarks/phase1/remove_low_mapq/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 8000
    shell:
        "samtools view -b --min-MQ 20 -o {output} {input} 2>&1 | tee {log}"


rule index_bam:
    input:
        "data/phase1/{mapq}/{sample}.bam"
    output:
        "data/phase1/{mapq}/{sample}.bam.bai"
    log:
        "logs/phase1/index_bam/{mapq}/{sample}.txt"
    benchmark:
        "benchmarks/phase1/index_bam/{mapq}/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 4000
    threads:
        2
    shell:
        "samtools index -@ {threads} {input} 2>&1 | tee {log}"


####################################################################
# Preprocessing statistics
####################################################################

rule count_raw_reads:
    input:
        FASTQ_DIR + "/{sample}.fastq.gz"
    output:
        "data/phase1/statistics/Raw_fastq/{sample}.txt"
    log:
        "logs/phase1/count_raw_reads/{sample}.txt"
    benchmark:
        "benchmarks/phase1/count_raw_reads/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 1000
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/read_count_fastq.gz.sh {input} {output} 2>&1 | tee {log}"


rule count_trimmed_reads:
    input:
        "data/phase1/Trimmed_fastq/{sample}.fastq.gz"
    output:
        "data/phase1/statistics/Trimmed_fastq/{sample}.txt"
    log:
        "logs/phase1/count_trimmed_reads/{sample}.txt"
    benchmark:
        "benchmarks/phase1/count_trimmed_reads/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 1000
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/read_count_fastq.gz.sh {input} {output} 2>&1 | tee {log}"


rule bam_flags:
    input:
        "data/phase1/{stage}/{sample}.bam"
    output:
        "data/phase1/statistics/bam_flags/{stage}/{sample}.txt"
    wildcard_constraints:
        stage = "Mapped_bam|All_reads|High_MapQ"
    log:
        "logs/phase1/bam_flags/{stage}/{sample}.txt"
    benchmark:
        "benchmarks/phase1/bam_flags/{stage}/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/bam_flags.sh {input} {output} temp 2> {log}"


rule mapping_quality:
    input:
        "data/phase1/All_reads/{sample}.bam"
    output:
        "data/phase1/statistics/mapping_quality/{sample}.txt"
    log:
        "logs/phase1/mapping_quality/{sample}.txt"
    benchmark:
        "benchmarks/phase1/mapping_quality/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 4000
    shell:
        "samtools view {input} | cut -f 5 | sort -T temp/ | uniq -c > {output} 2> {log}"


rule estimate_unique_maps:
    input:
        "data/phase1/All_reads/{sample}.bam"
    output:
        "data/phase1/statistics/unique_maps/{sample}.txt"
    log:
        "logs/phase1/estimate_unique_maps/{sample}.txt"
    benchmark:
        "benchmarks/phase1/estimate_unique_maps/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/Estimate_unique_maps.sh {input} temp > {output} 2> {log}"


rule compile_preprocessing_statistics:
    input:
        raw     = expand("data/phase1/statistics/Raw_fastq/{sample}.txt", sample=SAMPLES),
        trimmed = expand("data/phase1/statistics/Trimmed_fastq/{sample}.txt", sample=SAMPLES),
        flags   = expand("data/phase1/statistics/bam_flags/{stage}/{sample}.txt",
                         stage=BAM_STAGES, sample=SAMPLES),
        mapq    = expand("data/phase1/statistics/mapping_quality/{sample}.txt", sample=SAMPLES),
        unique  = expand("data/phase1/statistics/unique_maps/{sample}.txt", sample=SAMPLES),
        samples = "Zebrafinch_samples.txt",
    output:
        "results/phase1/Preprocessing_statistics.csv"
    log:
        "logs/phase1/compile_preprocessing_statistics.txt"
    benchmark:
        "benchmarks/phase1/compile_preprocessing_statistics/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase1_preprocess_and_correlate/Compile_preprocessing_statistics.R {input.samples} {output} 2>&1 | tee {log}"


####################################################################
# Window quantification and correlation
####################################################################

rule make_genome_windows:
    input:
        "resources/genomes/{genome}.chrom.sizes"
    output:
        "data/phase1/genome_windows/{size}_nucleotides/{genome}.bed"
    log:
        "logs/phase1/make_genome_windows/{size}_nucleotides/{genome}.txt"
    benchmark:
        "benchmarks/phase1/make_genome_windows/{size}_nucleotides/{genome}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 2000
    shell:
        "bash workflow/scripts/phase1_preprocess_and_correlate/make_windows.sh {input} {output} {wildcards.size} 2>&1 | tee {log}"


rule genome_window_coverage:
    input:
        bam     = "data/phase1/{mapq}/{sample}.bam",
        windows = "data/phase1/genome_windows/{size}_nucleotides/" + GENOME + ".bed",
    output:
        "data/phase1/genome_window_coverage/{mapq}/{size}_nucleotides/{sample}.bed"
    log:
        "logs/phase1/genome_window_coverage/{mapq}/{size}_nucleotides/{sample}.txt"
    benchmark:
        "benchmarks/phase1/genome_window_coverage/{mapq}/{size}_nucleotides/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 8000
    shell:
        # --op sum gives the summed per-base coverage in each window, which
        # is what the manuscript's compile_correlation.R calls "reads".
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum "
        "> {output} 2> {log}"


rule compile_correlation:
    input:
        coverage = expand(
            "data/phase1/genome_window_coverage/{mapq}/{size}_nucleotides/{sample}.bed",
            mapq=MAPQ, size=WINDOWS, sample=SAMPLES
        ),
        samples = "Zebrafinch_samples.txt",
    output:
        correlation = "results/phase1/Window_coverage_correlation.csv",
        matrix      = expand("results/phase1/Window_coverage_{size}_CPM.csv", size=[PRIMARY_WINDOW]),
    params:
        primary_window = PRIMARY_WINDOW,
        windows        = ",".join(WINDOWS),
        mapq           = ",".join(MAPQ),
    log:
        "logs/phase1/compile_correlation.txt"
    benchmark:
        "benchmarks/phase1/compile_correlation/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        "Rscript workflow/scripts/phase1_preprocess_and_correlate/Compile_correlation.R {input.samples} {params.windows} "
        "{params.mapq} {params.primary_window} {output.correlation} {output.matrix} "
        "2>&1 | tee {log}"


rule plot_correlation:
    input:
        "results/phase1/Window_coverage_correlation.csv"
    output:
        "plots/phase1/Correlation_heatmap.svg"
    log:
        "logs/phase1/plot_correlation.txt"
    benchmark:
        "benchmarks/phase1/plot_correlation/all.tsv"
    conda:
        # Rplot.yml, not Rutility.yml: ggsave() to .svg needs svglite, and
        # adding it to Rutility.yml would rehash that env and rerun the two
        # completed compile_* jobs.
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase1_preprocess_and_correlate/Plot_correlation.R {input} {output} 2>&1 | tee {log}"


####################################################################
# Browser tracks
####################################################################

rule bam_to_bedgraph:
    input:
        bam = "data/phase1/All_reads/{sample}.bam"
    output:
        temp("data/phase1/bedgraph/{sample}.bg")
    log:
        "logs/phase1/bam_to_bedgraph/{sample}.txt"
    benchmark:
        "benchmarks/phase1/bam_to_bedgraph/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 8000
    shell:
        "bedtools genomecov -bg -split -ibam {input.bam} > {output} 2> {log}"


rule bedgraph_to_bigwig:
    input:
        bg          = "data/phase1/bedgraph/{sample}.bg",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase1/bigwig/{sample}.bigwig"
    log:
        "logs/phase1/bedgraph_to_bigwig/{sample}.txt"
    benchmark:
        "benchmarks/phase1/bedgraph_to_bigwig/{sample}.tsv"
    conda:
        "../env/phase1_preprocess.yml"
    resources:
        mem_mb = 12000
    threads:
        2
    shell:
        # bedtools emits the bedgraph in BAM header (i.e. .fai) order, but
        # bedGraphToBigWig requires ASCII chromosome order, and the .fai is
        # chr1, chr1A, chr2, ... not chr1, chr10, chr11. The sort is not
        # optional here.
        "bash workflow/scripts/phase1_preprocess_and_correlate/bg_to_bigwig.sh -b {input.bg} -c {input.chrom_sizes} "
        "-o {output} -d temp 2>&1 | tee {log}"
