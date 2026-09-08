####################################################################
# Phase 3 -- Determine the optimum number of HMM states.
#
# Adapted from manuscript step 09 (09_Human_HMM.smk) and the model
# scripts it borrows from step 07 (07_PhyloHGMP_model).
#
# The manuscript had four human cell lines, so its Gaussian HMM emission
# was four dimensional. We have one, so ours is univariate. Everything
# below is written over HMM_FEATURES rather than a fixed column, so
# adding the WGS control or the 0 mM library as a second dimension is a
# Snakefile edit and no rule changes.
#
# This phase stops at the diagnostics. Which k to use is a judgement made
# by eye from the BIC scree, the k-means scree and the Gaussian fit
# plots; the functional interpretation of the states is Phase 4.
####################################################################

import os

HMM_DIR = "data/phase3/GHMM"


####################################################################
# Window signal, blacklist filtered
####################################################################

rule hmm_filter_windows:
    # One rule for what the manuscript splits between
    # Remove_human_black_list and Remove_human_sex_chromosomes: our Phase 2
    # blacklist already carries chrW_mat and chrZ_pat, so the awk only has to
    # drop chrMT and the three rDNA morphs.
    #
    # awk, not `grep -v`. grep exits 1 when it filters nothing, and a non-zero
    # exit anywhere in a Snakemake shell block fails the job -- the same trap
    # that add_sex_chromosomes_to_blacklist.sh works around.
    input:
        coverage = "data/phase2/experiment_window_coverage/All_reads/"
                   + HMM_WINDOW + "_nucleotides/{experiment}.bed",
        blacklist = "results/phase2/Blacklist.bed",
    output:
        "data/phase3/filtered_coverage/{experiment}.bg"
    params:
        exclude = ",".join(EXCLUDE_CHROMOSOMES)
    log:
        "logs/phase3/hmm_filter_windows/{experiment}.txt"
    benchmark:
        "benchmarks/phase3/hmm_filter_windows/{experiment}.tsv"
    conda:
        "../env/phase2_nonb.yml"
    resources:
        mem_mb = 4000
    shell:
        # -a order is preserved, so the output stays in .fai order.
        "bedtools intersect -v -a {input.coverage} -b {input.blacklist} 2> {log} "
        "| awk -v drop='{params.exclude}' "
        "'BEGIN {{ n = split(drop, a, \",\"); for (i = 1; i <= n; i++) skip[a[i]] = 1 }} "
        "!($1 in skip)' > {output}"


rule normalize_hmm_coverage:
    input:
        expand("data/phase3/filtered_coverage/{experiment}.bg",
               experiment=HMM_FEATURES)
    output:
        "data/phase3/HMM_input.csv"
    params:
        features  = ",".join(HMM_FEATURES),
        transform = HMM_TRANSFORM,
    log:
        "logs/phase3/normalize_hmm_coverage.txt"
    benchmark:
        "benchmarks/phase3/normalize_hmm_coverage/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 8000
    shell:
        # HMM_TRANSFORM is a params value, and snakemake counts params as a
        # rerun trigger by default, so flipping it to log1p re-runs this rule
        # and everything after it without touching Phase 1 or Phase 2.
        "Rscript workflow/scripts/phase3_hmm_states/Normalize_hmm_coverage.R "
        "{params.features} {params.transform} {output} 2>&1 | tee {log}"


####################################################################
# k-means scree -- the model-free read on how many clusters there are
####################################################################

rule kmeans_scree:
    input:
        "data/phase3/HMM_input.csv"
    output:
        "results/phase3/Kmeans_scree.csv"
    params:
        features = ",".join(HMM_FEATURES),
        states   = ",".join(HMM_STATES),
    log:
        "logs/phase3/kmeans_scree.txt"
    benchmark:
        "benchmarks/phase3/kmeans_scree/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 12000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Kmeans_scree.R "
        "{input} {params.features} {params.states} {output} 2>&1 | tee {log}"


rule plot_kmeans_scree:
    input:
        "results/phase3/Kmeans_scree.csv"
    output:
        "plots/phase3/Kmeans_scree.svg"
    log:
        "logs/phase3/plot_kmeans_scree.txt"
    benchmark:
        "benchmarks/phase3/plot_kmeans_scree/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Plot_Kmeans_scree.R "
        "{input} {output} 2>&1 | tee {log}"


####################################################################
# Gaussian HMM, one fit per state count
####################################################################

rule fit_ghmm:
    input:
        "data/phase3/HMM_input.csv"
    output:
        states = HMM_DIR + "/{states}/GHMM_states.csv.gz",
        pkl    = HMM_DIR + "/{states}/GHMM.pkl",
    params:
        features = ",".join(HMM_FEATURES)
    log:
        "logs/phase3/fit_ghmm/{states}.txt"
    benchmark:
        "benchmarks/phase3/fit_ghmm/{states}.tsv"
    conda:
        "../env/hmm.yml"
    resources:
        # The forward/backward lattices dominate: ~10^6 windows x k x 8 bytes,
        # several copies live at once. Measured on 10^6 windows x 1 feature:
        # 1.19 GB peak RSS at k = 10, 2.62 GB at k = 30. This leaves roughly a
        # 4x margin (10.5 GB at k = 30) without over-reserving, which matters
        # because mem_mb is a global cap and every GB claimed here is a fit
        # that cannot start alongside.
        mem_mb = lambda wildcards: 3000 + 250 * int(wildcards.states)
    shell:
        "python workflow/scripts/phase3_hmm_states/GHMM.py "
        "{input} {params.features} {wildcards.states} "
        "{output.states} {output.pkl} 2>&1 | tee {log}"


rule extract_ghmm_parameters:
    input:
        pkl    = HMM_DIR + "/{states}/GHMM.pkl",
        states = HMM_DIR + "/{states}/GHMM_states.csv.gz",
    output:
        parameters = HMM_DIR + "/{states}/Parameters.csv",
        histogram  = HMM_DIR + "/{states}/Histogram.csv",
    params:
        features = ",".join(HMM_FEATURES)
    log:
        "logs/phase3/extract_ghmm_parameters/{states}.txt"
    benchmark:
        "benchmarks/phase3/extract_ghmm_parameters/{states}.tsv"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb = 8000
    shell:
        "python workflow/scripts/phase3_hmm_states/extract_GHMM_parameters.py "
        "{input.pkl} {input.states} {params.features} "
        "{output.parameters} {output.histogram} 2>&1 | tee {log}"


####################################################################
# BIC scree
####################################################################

rule hmm_bic:
    input:
        # Parameters.csv is what supplies occupied_states, so it is a real
        # dependency and not just a convenience.
        pkl = expand(HMM_DIR + "/{states}/GHMM.pkl", states=HMM_STATES),
        parameters = expand(HMM_DIR + "/{states}/Parameters.csv", states=HMM_STATES),
        data = "data/phase3/HMM_input.csv",
    output:
        "results/phase3/HMM_BIC.csv"
    params:
        # Derived from the inputs rather than hardcoded, so the path is
        # correct wherever the job runs: dirname twice off any pkl gives the
        # directory that holds the per-state model directories.
        model_dir = lambda wildcards, input: os.path.dirname(
            os.path.dirname(input.pkl[0])
        ),
        states    = ",".join(HMM_STATES),
        features  = ",".join(HMM_FEATURES),
    log:
        "logs/phase3/hmm_bic.txt"
    benchmark:
        "benchmarks/phase3/hmm_bic/all.tsv"
    conda:
        "../env/hmm.yml"
    resources:
        # Scores one model at a time, so the peak is a single k = 30 forward
        # pass over ~10^6 windows -- the same 2.6 GB the fit itself needs.
        mem_mb = 8000
    shell:
        "python workflow/scripts/phase3_hmm_states/calc_BIC.py "
        "{params.model_dir} {params.states} {input.data} {params.features} "
        "{output} 2>&1 | tee {log}"


rule plot_hmm_bic:
    input:
        "results/phase3/HMM_BIC.csv"
    output:
        "plots/phase3/HMM_BIC_scree.svg"
    log:
        "logs/phase3/plot_hmm_bic.txt"
    benchmark:
        "benchmarks/phase3/plot_hmm_bic/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Plot_BIC_scree.R "
        "{input} {output} 2>&1 | tee {log}"


####################################################################
# Is the model fitting the data?
####################################################################

rule plot_gaussian_fit:
    input:
        parameters = HMM_DIR + "/{states}/Parameters.csv",
        histogram  = HMM_DIR + "/{states}/Histogram.csv",
    output:
        "plots/phase3/gaussian_fit/{states}_states.svg"
    log:
        "logs/phase3/plot_gaussian_fit/{states}.txt"
    benchmark:
        "benchmarks/phase3/plot_gaussian_fit/{states}.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        # Cheap because extract_GHMM_parameters.py already binned the signal;
        # this reads a few hundred rows, not ~10^6 windows.
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Plot_gaussian_fit.R "
        "{input.parameters} {input.histogram} {output} 2>&1 | tee {log}"


rule compile_state_summary:
    input:
        expand(HMM_DIR + "/{states}/Parameters.csv", states=HMM_STATES)
    output:
        "results/phase3/HMM_state_summary.csv"
    log:
        "logs/phase3/compile_state_summary.txt"
    benchmark:
        "benchmarks/phase3/compile_state_summary/all.tsv"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Compile_state_summary.R "
        "{output} {input} 2>&1 | tee {log}"


rule plot_model_vs_empirical_means:
    input:
        "results/phase3/HMM_state_summary.csv"
    output:
        "plots/phase3/HMM_model_vs_empirical_means.svg"
    log:
        "logs/phase3/plot_model_vs_empirical_means.txt"
    benchmark:
        "benchmarks/phase3/plot_model_vs_empirical_means/all.tsv"
    conda:
        "../env/Rplot.yml"
    resources:
        mem_mb = 4000
    shell:
        "Rscript workflow/scripts/phase3_hmm_states/Plot_model_vs_empirical_means.R "
        "{input} {output} 2>&1 | tee {log}"
