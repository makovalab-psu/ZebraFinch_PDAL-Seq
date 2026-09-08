# All-pairs Pearson and Spearman correlation of window coverage between
# every pair of Zebra finch sequencing samples, at each window size and
# for both MapQ branches. Also writes the wide coverage matrix that
# Phase 2 consumes.
#
# Adapted from
#   workflow/scripts/03_Correlation_and_black_list_unmappable_regions/compile_correlation.R
#   workflow/scripts/21_add_K562_to_analysis/compile_K562_correlation.R
# in the PDAL-Seq manuscript pipeline. Changes: one genome instead of a
# per-species split, an added MapQ-branch dimension, an added Comparison
# label, and a check that the window rows actually line up between files.
#
# Note on units: megadepth --op sum reports the summed per-base coverage
# in each window, not a read count. The manuscript's script calls this
# column "reads"; that name is kept for continuity, but the normalized
# matrix is labelled CPM (coverage per million) because that is what it
# is: window coverage divided by the genome-wide coverage total, times 1e6.
#
# Usage:
#   Rscript Compile_correlation.R <samples.txt> <windows_csv> <mapq_csv> \
#           <primary_window> <out_correlation.csv> <out_matrix.csv>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
	stop("Usage: Compile_correlation.R <samples.txt> <windows_csv> <mapq_csv> <primary_window> <out_correlation.csv> <out_matrix.csv>")
}

samples_file    <- args[1]
Windows         <- strsplit(args[2], ",")[[1]]
Mapq_branches   <- strsplit(args[3], ",")[[1]]
primary_window  <- args[4]
out_correlation <- args[5]
out_matrix      <- args[6]

# The wide matrix is written from the All_reads branch, matching the
# manuscript convention that all downstream analyses run on deduplicated
# BAMs with no MapQ filter.
matrix_branch <- "All_reads"

coverage_dir <- "data/phase1/genome_window_coverage"

df.index <- read.delim(samples_file, stringsAsFactors = FALSE)
samples <- df.index$Sequencing_sample

get.path <- function(mapq, window, sample) {
	file.path(coverage_dir, mapq, paste0(window, "_nucleotides"), paste0(sample, ".bed"))
}

read.bg <- function(x) {
	df <- read.delim(x, sep = "", header = FALSE)
	colnames(df) <- c("chr", "start", "end", "reads")
	df
}

# A short label for plots and for reading the correlation table by eye.
short.label <- function(sample) {
	row <- df.index[df.index$Sequencing_sample == sample, ]
	if (row$Library == "PCRfree") {
		"WGS"
	} else {
		paste0(row$Treatment, "-", row$Replicate)
	}
}

labels <- setNames(vapply(samples, short.label, character(1)), samples)
library_of   <- setNames(df.index$Library, df.index$Sequencing_sample)
treatment_of <- setNames(df.index$Treatment, df.index$Sequencing_sample)
experiment_of <- setNames(df.index$Experiment, df.index$Sequencing_sample)

comparison.label <- function(s1, s2) {
	libs <- sort(c(library_of[[s1]], library_of[[s2]]))
	if (identical(libs, c("PCRfree", "PCRfree"))) {
		"WGS_vs_WGS"
	} else if (identical(libs, c("PDALSeq", "PDALSeq"))) {
		"PDALSeq_vs_PDALSeq"
	} else {
		"PDALSeq_vs_WGS"
	}
}

# ---------------------------------------------------------------------
# Correlation over every mapq branch x window size x sample pair
# ---------------------------------------------------------------------

results <- list()
k <- 0

for (mapq in Mapq_branches) {
	for (window in Windows) {

		message("loading ", mapq, " ", window, " nucleotide windows")

		# Load all samples for this branch/window once. At 10 kbp there are
		# ~114,000 windows over the 1.14 Gbp assembly, so the whole set is
		# a few MB and holding it is cheaper than re-reading per pair.
		cov <- list()
		coords <- NULL
		for (s in samples) {
			df.s <- read.bg(get.path(mapq, window, s))
			if (is.null(coords)) {
				coords <- df.s[, c("chr", "start", "end")]
			} else if (!identical(coords$chr, df.s$chr) || !identical(coords$start, df.s$start)) {
				stop("window rows do not line up between samples for ", mapq, " ", window,
				     " (offending sample: ", s, ")")
			}
			cov[[s]] <- df.s$reads
		}

		for (s1 in samples) {
			for (s2 in samples) {
				k <- k + 1
				t.pearson  <- suppressWarnings(cor.test(cov[[s1]], cov[[s2]], method = "pearson"))
				t.spearman <- suppressWarnings(cor.test(cov[[s1]], cov[[s2]], method = "spearman"))

				results[[k]] <- data.frame(
					MapQ_branch     = mapq,
					Window          = window,
					Sample_1        = s1,
					Sample_2        = s2,
					Label_1         = labels[[s1]],
					Label_2         = labels[[s2]],
					Library_1       = library_of[[s1]],
					Library_2       = library_of[[s2]],
					Treatment_1     = treatment_of[[s1]],
					Treatment_2     = treatment_of[[s2]],
					Comparison      = comparison.label(s1, s2),
					Same_experiment = experiment_of[[s1]] == experiment_of[[s2]],
					Pearson         = unname(t.pearson$estimate),
					Spearman        = unname(t.spearman$estimate),
					stringsAsFactors = FALSE
				)
			}
		}

		# The wide matrix Phase 2 consumes, written while this branch and
		# window are already in memory.
		if (mapq == matrix_branch && window == primary_window) {
			df.matrix <- coords
			total <- vapply(samples, function(s) sum(cov[[s]]), numeric(1))
			for (s in samples) {
				df.matrix[[s]] <- 1e6 * cov[[s]] / total[[s]]
			}
			dir.create(dirname(out_matrix), recursive = TRUE, showWarnings = FALSE)
			write.csv(df.matrix, out_matrix, row.names = FALSE)
			message("wrote ", out_matrix, " (", nrow(df.matrix), " windows x ",
			        length(samples), " samples, CPM)")
		}
	}
}

df <- bind_rows(results)

dir.create(dirname(out_correlation), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out_correlation, row.names = FALSE)

message("wrote ", out_correlation, " (", nrow(df), " pairs)")
