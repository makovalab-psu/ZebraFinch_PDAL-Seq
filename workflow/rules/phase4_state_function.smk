####################################################################
# Phase 4 -- Intersect HMM states with non-B DNA motifs, 5mC rates,
#            RNA-Seq and functional genomic annotations.
#
# Adapted from the PDAL-Seq manuscript pipeline:
#   rules/07_PhyloHGMP_model.smk    (state annotation, Enrichment.sh,
#                                    BG_enrichment.sh, the compile pattern)
#   rules/10_Methylation_analysis.smk (5mC rates in states)
#   rules/18_ATACseq_analysis.smk     (read mapping and filtering protocol)
#   scripts/Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R (the plot layout)
#
# Nothing is re-fit here. Phase 3 already fit every model in PHASE4_STATES, so
# this phase reads data/phase3/GHMM/{k}/GHMM_states.csv.gz and Parameters.csv
# and asks what the states MEAN. Every rule carries {states} as a wildcard, so
# adding or removing a k is a one-line edit to PHASE4_STATES.
#
# Two protocol deviations, both discussed in README section 06:
#   * RNA-Seq is mapped with STAR, not bwa mem. bwa mem is not splice aware.
#   * There is no duplicate-removal step. The ATAC-Seq protocol runs
#     `samtools rmdup -s`, which is right for genomic fragments; on RNA-Seq a
#     highly expressed gene produces many independent fragments from the same
#     start position, so rmdup would strip the most reads from exactly the
#     regions this analysis is trying to detect.
#
# The enrichment null shuffles across the whole genome, as the manuscript does,
# rather than being restricted to the segmented genome.
####################################################################


####################################################################
# RNA-Seq: STAR index
####################################################################

rule star_index:
    input:
        fasta = GENOME_FASTA
    output:
        # Naming the files STAR always writes, rather than directory(), so a
        # partial index from a killed job cannot satisfy the DAG.
        genome  = STAR_INDEX + "/Genome",
        sa      = STAR_INDEX + "/SA",
        saindex = STAR_INDEX + "/SAindex",
        sizes   = STAR_INDEX + "/chrNameLength.txt",
    params:
        index_dir = STAR_INDEX,
        # min(14, log2(genomeLength)/2 - 1). This assembly is 1.14 Gbp, which
        # puts the formula at 14.05, so 14 -- STAR's default -- is correct. It
        # is written out because the default is tuned for the human genome and
        # silently wrong for a small one.
        sa_index_nbases = 14,
    log:
        "logs/phase4/star_index/" + GENOME + ".txt"
    benchmark:
        "benchmarks/phase4/star_index/" + GENOME + ".tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        # Peak is the suffix array sort, roughly 10 bytes per base plus
        # workspace. 40 GB is ~2x what a 1.14 Gbp genome needs and still fits
        # 10 cores on the 4 GB/core basic partition.
        mem_mb = 40000
    threads:
        20
    shell:
        "mkdir -p {params.index_dir} && "
        "STAR --runMode genomeGenerate "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--genomeSAindexNbases {params.sa_index_nbases} "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.index_dir}/ 2>&1 | tee {log}"


####################################################################
# RNA-Seq: trim, map, quantify
####################################################################

rule rna_trim_fastq:
    input:
        R1 = RNA_DIR + "/{sra}_1.fastq.gz",
        R2 = RNA_DIR + "/{sra}_2.fastq.gz",
    output:
        R1   = temp("data/phase4/RNA_seq/trimmed/{sra}_1.fastq.gz"),
        R2   = temp("data/phase4/RNA_seq/trimmed/{sra}_2.fastq.gz"),
        json = "results/phase4/fastp_reports/{sra}.json",
        html = "results/phase4/fastp_reports/{sra}.html",
    params:
        read_limit = RNA_READ_LIMIT_FLAG
    log:
        "logs/phase4/rna_trim_fastq/{sra}.txt"
    benchmark:
        "benchmarks/phase4/rna_trim_fastq/{sra}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 8000
    threads:
        16
    shell:
        # Flags are identical to ATACseq_trim_fastq in the manuscript's
        # 18_ATACseq_analysis.smk. fastp ignores --thread above 16.
        "fastp --thread {threads} --in1 {input.R1} --in2 {input.R2} "
        "--out1 {output.R1} --out2 {output.R2} "
        "-p -P 50 -y -x -3 {params.read_limit} "
        "--json {output.json} --html {output.html} 2>&1 | tee {log}"


rule rna_map:
    input:
        R1      = "data/phase4/RNA_seq/trimmed/{sra}_1.fastq.gz",
        R2      = "data/phase4/RNA_seq/trimmed/{sra}_2.fastq.gz",
        genome  = STAR_INDEX + "/Genome",
        sa      = STAR_INDEX + "/SA",
        saindex = STAR_INDEX + "/SAindex",
        annotation = rna_annotation_input,
    output:
        # NOT temp(). With no duplicate-removal step this is the analysis BAM,
        # and regenerating it costs another ~8-10 h of STAR.
        bam       = "data/phase4/RNA_seq/sorted/{sra}.bam",
        log_final = "results/phase4/STAR_logs/{sra}.Log.final.out",
        junctions = "results/phase4/STAR_logs/{sra}.SJ.out.tab",
    params:
        index_dir  = STAR_INDEX,
        prefix     = "results/phase4/STAR_logs/{sra}.",
        annotation = RNA_ANNOTATION_FLAG,
    log:
        "logs/phase4/rna_map/{sra}.txt"
    benchmark:
        "benchmarks/phase4/rna_map/{sra}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        # STAR holds the ~12-14 GB index for the whole run and samtools sort
        # runs concurrently in the same pipe at a fixed 6 GB (see
        # map_rna_reads.sh). 28 GB covers both with headroom; it is also 7
        # cores' worth on basic, so two of these cannot run at once on a
        # 20-core allocation, which is what we want.
        mem_mb = 28000
    threads:
        20
    shell:
        "bash workflow/scripts/phase4_state_function/map_rna_reads.sh "
        "-1 {input.R1} -2 {input.R2} -x {params.index_dir} "
        "-o {output.bam} -p {params.prefix} -t {threads} -d temp "
        "{params.annotation} 2>&1 | tee {log}"


rule rna_flagstat:
    input:
        "data/phase4/RNA_seq/sorted/{sra}.bam"
    output:
        "data/phase4/RNA_seq/flagstat/{sra}.txt"
    log:
        "logs/phase4/rna_flagstat/{sra}.txt"
    benchmark:
        "benchmarks/phase4/rna_flagstat/{sra}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 4000
    threads:
        4
    shell:
        "samtools flagstat -@ {threads} {input} > {output} 2> {log}"


rule rna_signal:
    input:
        bam     = "data/phase4/RNA_seq/sorted/{sra}.bam",
        windows = "data/phase1/genome_windows/" + HMM_WINDOW
                  + "_nucleotides/" + GENOME + ".bed",
    output:
        "data/phase4/RNA_signal/{sra}.bg"
    log:
        "logs/phase4/rna_signal/{sra}.txt"
    benchmark:
        "benchmarks/phase4/rna_signal/{sra}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 8000
    shell:
        # The same 1 kbp grid the HMM was fit to, and the same --op sum
        # quantity Phase 1 and Phase 2 use, so "RNA-Seq signal in a state" is
        # measured on the identical windows as "PDAL-Seq signal in a state".
        # See the header of BG_enrichment.sh for why this is not a raw
        # `bedtools genomecov -bg`.
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum "
        "> {output} 2> {log}"


rule rna_mapping_report:
    input:
        star     = expand("results/phase4/STAR_logs/{sra}.Log.final.out", sra=RNA_SRA),
        flagstat = expand("data/phase4/RNA_seq/flagstat/{sra}.txt", sra=RNA_SRA),
    output:
        "results/phase4/RNA_seq_mapping_statistics.txt"
    params:
        sra = " ".join(RNA_SRA)
    log:
        "logs/phase4/rna_mapping_report.txt"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 1000
    shell:
        "( for S in {params.sra}; do "
        "    echo \"### ${{S}} -- STAR\"; "
        "    cat results/phase4/STAR_logs/${{S}}.Log.final.out; "
        "    echo \"### ${{S}} -- samtools flagstat\"; "
        "    cat data/phase4/RNA_seq/flagstat/${{S}}.txt; "
        "    echo; "
        "  done ) > {output} 2> {log}"


####################################################################
# Functional annotations
####################################################################

rule clean_annotation:
    input:
        # An input function, so it cannot be wrapped in ancient() the way the
        # literal paths elsewhere are -- ancient() takes a path, not a callable.
        bed         = annotation_source,
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        bed    = "results/phase4/annotations/{ann}.bed",
        report = "data/phase4/annotations/{ann}.report.txt",
    log:
        "logs/phase4/clean_annotation/{ann}.txt"
    benchmark:
        "benchmarks/phase4/clean_annotation/{ann}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase4_state_function/clean_annotation.sh "
        "-i {input.bed} -g {input.chrom_sizes} -o {output.bed} "
        "-r {output.report} -n {wildcards.ann} -d temp 2>&1 | tee {log}"


rule annotation_report:
    input:
        expand("data/phase4/annotations/{ann}.report.txt", ann=ANNOTATIONS)
    output:
        "results/phase4/Annotation_report.txt"
    log:
        "logs/phase4/annotation_report.txt"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 1000
    shell:
        "( printf 'Annotation\\tInput_records\\tRetained_records\\tPercent_retained\\tContigs_retained\\n' ; "
        "cat {input} ) > {output} 2> {log}"


####################################################################
# HMM state segments
####################################################################

rule segmented_genome:
    input:
        "data/phase3/HMM_input.csv"
    output:
        "results/phase4/Segmented_genome.bed"
    log:
        "logs/phase4/segmented_genome.txt"
    benchmark:
        "benchmarks/phase4/segmented_genome/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Segmented_genome.R "
        "{input} {output} 2>&1 | tee {log}"


# One rule per state count rather than one rule with a {states} wildcard.
# Snakemake evaluates `output:` at parse time, before wildcards exist, so
# `expand(..., state=range(int(wildcards.states)))` cannot appear there. The
# state counts are literals in PHASE4_STATES, so the list is known at parse
# time and the loop bakes it in -- the same reason the skill generates
# per-genome rules when a chromosome set is wildcard derived.
#
# The alternative, a .done marker with the BEDs referenced through params,
# would hide 77 real files from the DAG for no gain here.
for PHASE4_K in PHASE4_STATES:

    rule:
        name:
            "annotate_states_{}".format(PHASE4_K)
        input:
            PHASE3_HMM_DIR + "/{}/GHMM_states.csv.gz".format(PHASE4_K)
        output:
            expand(STATE_DIR + "/{states}/state_{state}.bed",
                   states=PHASE4_K, state=range(int(PHASE4_K)))
        params:
            n_states = PHASE4_K,
            out_dir  = STATE_DIR + "/" + PHASE4_K,
        log:
            "logs/phase4/annotate_states/{}.txt".format(PHASE4_K)
        benchmark:
            "benchmarks/phase4/annotate_states/{}.tsv".format(PHASE4_K)
        conda:
            "../env/Rutility.yml"
        resources:
            mem_mb = 8000
        shell:
            "Rscript workflow/scripts/phase4_state_function/Annotate_states.R "
            "{input} {params.n_states} {params.out_dir} 2>&1 | tee {log}"


####################################################################
# Enrichment, one job per state x annotation
####################################################################

rule nonb_enrichment_in_state:
    input:
        state       = STATE_DIR + "/{states}/state_{state}.bed",
        nonb        = "results/phase2/nonb_annotations/{nonb}.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase4/nonB_enrichment/{states}/state_{state}/{nonb}.csv"
    params:
        subsamples = ENRICHMENT_SUBSAMPLES,
        iterations = ENRICHMENT_ITERATIONS,
    log:
        "logs/phase4/nonb_enrichment_in_state/{states}/state_{state}/{nonb}.txt"
    benchmark:
        "benchmarks/phase4/nonb_enrichment_in_state/{states}/state_{state}/{nonb}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        # bedtools holds only the 10,000-record subsample in memory; the state
        # BED streams past it.
        mem_mb = 2000
    shell:
        "bash workflow/scripts/phase4_state_function/Enrichment.sh "
        "-s {input.state} -a {input.nonb} -g {input.chrom_sizes} "
        "-o {output} -n {params.subsamples} -i {params.iterations} "
        "-d temp 2>&1 | tee {log}"


rule functional_enrichment_in_state:
    input:
        state       = STATE_DIR + "/{states}/state_{state}.bed",
        annotation  = "results/phase4/annotations/{ann}.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase4/functional_enrichment/{states}/state_{state}/{ann}.csv"
    params:
        subsamples = ENRICHMENT_SUBSAMPLES,
        iterations = ENRICHMENT_ITERATIONS,
    log:
        "logs/phase4/functional_enrichment_in_state/{states}/state_{state}/{ann}.txt"
    benchmark:
        "benchmarks/phase4/functional_enrichment_in_state/{states}/state_{state}/{ann}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 2000
    shell:
        "bash workflow/scripts/phase4_state_function/Enrichment.sh "
        "-s {input.state} -a {input.annotation} -g {input.chrom_sizes} "
        "-o {output} -n {params.subsamples} -i {params.iterations} "
        "-d temp 2>&1 | tee {log}"


rule rna_enrichment_in_state:
    input:
        state       = STATE_DIR + "/{states}/state_{state}.bed",
        signal      = "data/phase4/RNA_signal/{sra}.bg",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        "data/phase4/RNA_enrichment/{states}/state_{state}/{sra}.csv"
    params:
        subsamples = RNA_ENRICHMENT_SUBSAMPLES,
        iterations = ENRICHMENT_ITERATIONS,
    log:
        "logs/phase4/rna_enrichment_in_state/{states}/state_{state}/{sra}.txt"
    benchmark:
        "benchmarks/phase4/rna_enrichment_in_state/{states}/state_{state}/{sra}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase4_state_function/BG_enrichment.sh "
        "-s {input.state} -a {input.signal} -g {input.chrom_sizes} "
        "-o {output} -n {params.subsamples} -i {params.iterations} "
        "-d temp 2>&1 | tee {log}"


rule methylation_in_state:
    input:
        state = STATE_DIR + "/{states}/state_{state}.bed",
        m5C   = ancient(METHYLATION_BED),
    output:
        "data/phase4/methylation/{states}/state_{state}.txt"
    log:
        "logs/phase4/methylation_in_state/{states}/state_{state}.txt"
    benchmark:
        "benchmarks/phase4/methylation_in_state/{states}/state_{state}.tsv"
    conda:
        "../env/phase4_functional.yml"
    resources:
        mem_mb = 4000
    shell:
        "bash workflow/scripts/phase4_state_function/methylation_in_state.sh "
        "-s {input.state} -m {input.m5C} -o {output} 2>&1 | tee {log}"


####################################################################
# Compile
####################################################################

def state_range(wildcards):
    """0 .. k-1 for the model this job is compiling."""
    return [str(i) for i in range(int(wildcards.states))]


rule compile_nonb_enrichment:
    input:
        lambda wildcards: expand(
            "data/phase4/nonB_enrichment/{states}/state_{state}/{nonb}.csv",
            states=wildcards.states, state=state_range(wildcards), nonb=NONB_ALL,
        )
    output:
        "results/phase4/{states}/nonB_enrichment.csv"
    params:
        in_dir = "data/phase4/nonB_enrichment/{states}",
        nonb   = ",".join(NONB_ALL),
    log:
        "logs/phase4/compile_nonb_enrichment/{states}.txt"
    benchmark:
        "benchmarks/phase4/compile_nonb_enrichment/{states}.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Compile_enrichment.R "
        "{params.in_dir} {wildcards.states} {params.nonb} nonB {output} "
        "2>&1 | tee {log}"


rule compile_functional_enrichment:
    input:
        lambda wildcards: expand(
            "data/phase4/functional_enrichment/{states}/state_{state}/{ann}.csv",
            states=wildcards.states, state=state_range(wildcards), ann=ANNOTATIONS,
        )
    output:
        "results/phase4/{states}/Functional_enrichment.csv"
    params:
        in_dir = "data/phase4/functional_enrichment/{states}",
        ann    = ",".join(ANNOTATIONS),
    log:
        "logs/phase4/compile_functional_enrichment/{states}.txt"
    benchmark:
        "benchmarks/phase4/compile_functional_enrichment/{states}.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Compile_enrichment.R "
        "{params.in_dir} {wildcards.states} {params.ann} ann {output} "
        "2>&1 | tee {log}"


rule compile_rna_enrichment:
    input:
        lambda wildcards: expand(
            "data/phase4/RNA_enrichment/{states}/state_{state}/{sra}.csv",
            states=wildcards.states, state=state_range(wildcards), sra=RNA_SRA,
        )
    output:
        "results/phase4/{states}/RNA_enrichment.csv"
    params:
        in_dir = "data/phase4/RNA_enrichment/{states}",
        sra    = ",".join(RNA_SRA),
    log:
        "logs/phase4/compile_rna_enrichment/{states}.txt"
    benchmark:
        "benchmarks/phase4/compile_rna_enrichment/{states}.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Compile_enrichment.R "
        "{params.in_dir} {wildcards.states} {params.sra} SRA {output} "
        "2>&1 | tee {log}"


rule compile_methylation:
    input:
        lambda wildcards: expand(
            "data/phase4/methylation/{states}/state_{state}.txt",
            states=wildcards.states, state=state_range(wildcards),
        )
    output:
        "results/phase4/{states}/Methylation_by_state.csv"
    params:
        in_dir = "data/phase4/methylation/{states}"
    log:
        "logs/phase4/compile_methylation/{states}.txt"
    benchmark:
        "benchmarks/phase4/compile_methylation/{states}.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Compile_methylation.R "
        "{params.in_dir} {wildcards.states} {output} 2>&1 | tee {log}"


rule phase4_state_summary:
    input:
        parameters  = PHASE3_HMM_DIR + "/{states}/Parameters.csv",
        states      = lambda wildcards: expand(
            STATE_DIR + "/{states}/state_{state}.bed",
            states=wildcards.states, state=state_range(wildcards),
        ),
        segmented   = "results/phase4/Segmented_genome.bed",
        chrom_sizes = "resources/genomes/" + GENOME + ".chrom.sizes",
    output:
        summary  = "results/phase4/{states}/State_summary.csv",
        by_chrom = "results/phase4/{states}/State_by_chromosome.csv",
    params:
        state_dir = STATE_DIR + "/{states}"
    log:
        "logs/phase4/phase4_state_summary/{states}.txt"
    benchmark:
        "benchmarks/phase4/phase4_state_summary/{states}.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Compile_state_summary.R "
        "{input.parameters} {params.state_dir} {wildcards.states} "
        "{input.segmented} {input.chrom_sizes} "
        "{output.summary} {output.by_chrom} 2>&1 | tee {log}"


####################################################################
# Plot
####################################################################

rule plot_state_function:
    input:
        summary     = "results/phase4/{states}/State_summary.csv",
        nonb        = "results/phase4/{states}/nonB_enrichment.csv",
        functional  = "results/phase4/{states}/Functional_enrichment.csv",
        rna         = "results/phase4/{states}/RNA_enrichment.csv",
        methylation = "results/phase4/{states}/Methylation_by_state.csv",
    output:
        "plots/phase4/{states}/HMM_functional_summary.svg"
    log:
        "logs/phase4/plot_state_function/{states}.txt"
    benchmark:
        "benchmarks/phase4/plot_state_function/{states}.tsv"
    conda:
        "../env/Rphase4.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Plot_state_function.R "
        "{input.summary} {input.nonb} {input.functional} {input.rna} "
        "{input.methylation} {wildcards.states} {output} 2>&1 | tee {log}"


rule plot_state_distribution:
    input:
        "results/phase4/{states}/State_by_chromosome.csv"
    output:
        "plots/phase4/{states}/State_distribution_by_chromosome.svg"
    log:
        "logs/phase4/plot_state_distribution/{states}.txt"
    benchmark:
        "benchmarks/phase4/plot_state_distribution/{states}.tsv"
    conda:
        "../env/Rphase4.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase4_state_function/Plot_state_distribution.R "
        "{input} {wildcards.states} {output} 2>&1 | tee {log}"
