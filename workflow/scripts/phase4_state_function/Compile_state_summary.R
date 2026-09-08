# Size, signal and chromosome distribution of every state in one fitted HMM.
#
# Covers what three manuscript pieces did separately:
#   07_PhyloHGMP_model/BP_in_state.sh              -- bp per state
#   07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R -- mean signal per state
#   the per-chromosome loop inside Figure_3_Ape_PDAL-Seq_MVGHMM.R
#
# The signal columns are not recomputed here: Phase 3's
# extract_GHMM_parameters.py already wrote the Gaussian and empirical means to
# Parameters.csv for exactly these models, and re-deriving them from the state
# BEDs would be a second implementation of the same number that could disagree
# with the fit plots.
#
# "Percent" in the summary is of the SEGMENTED genome (what any state could
# have covered), not of the assembly. "Percent" in the per-chromosome table is
# of the chromosome, and each chromosome carries an "unsegmented" row for the
# blacklisted remainder, so the bars in Plot_state_distribution.R reach 100%.
#
# Usage:
#   Rscript Compile_state_summary.R <Parameters.csv> <state_dir> <n_states> \
#           <segmented.bed> <chrom.sizes> <out_summary.csv> <out_by_chrom.csv>

suppressPackageStartupMessages(library(tidyverse))

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
	stop("Usage: Compile_state_summary.R <Parameters.csv> <state_dir> <n_states> <segmented.bed> <chrom.sizes> <out_summary.csv> <out_by_chrom.csv>")
}

parameters_csv <- args[1]
state_dir      <- args[2]
n_states       <- as.integer(args[3])
segmented_bed  <- args[4]
chrom_sizes    <- args[5]
out_summary    <- args[6]
out_by_chrom   <- args[7]

read_bed <- function(path) {
	# An unoccupied state leaves an empty file; read.delim would abort on it.
	if (!file.exists(path)) {
		stop("missing state BED: ", path)
	}
	if (file.size(path) == 0) {
		return(tibble::tibble(chr = character(), start = numeric(), end = numeric()))
	}
	readr::read_tsv(
		path,
		col_names = c("chr", "start", "end"),
		col_types = readr::cols(
			chr = readr::col_character(),
			start = readr::col_double(),
			end = readr::col_double()
		)
	)
}

df.sizes <- readr::read_tsv(
	chrom_sizes,
	col_names = c("chr", "length"),
	col_types = readr::cols(chr = readr::col_character(), length = readr::col_double())
)

df.segmented <- read_bed(segmented_bed)
segmented_bp <- sum(df.segmented$end - df.segmented$start + 1)

df.parameters <- readr::read_csv(parameters_csv, col_types = readr::cols())
df.parameters <- df.parameters %>% dplyr::mutate(State = as.character(state))

l.summary  <- list()
l.by_chrom <- list()

for (state in seq_len(n_states) - 1L) {

	df.state <- read_bed(file.path(state_dir, paste0("state_", state, ".bed")))

	bp <- if (nrow(df.state) == 0) 0 else sum(df.state$end - df.state$start + 1)

	l.summary[[length(l.summary) + 1]] <- tibble::tibble(
		n_states       = n_states,
		State          = as.character(state),
		Segments       = nrow(df.state),
		BP             = bp,
		Percent        = 100 * bp / segmented_bp,
		Segmented_BP   = segmented_bp
	)

	if (nrow(df.state) > 0) {
		l.by_chrom[[length(l.by_chrom) + 1]] <- df.state %>%
			dplyr::group_by(chr) %>%
			dplyr::summarize(BP = sum(end - start + 1), .groups = "drop") %>%
			dplyr::mutate(n_states = n_states, State = as.character(state))
	}
}

df.summary <- dplyr::bind_rows(l.summary)

# Carry the Phase 3 signal columns across. left_join, not cbind: Parameters.csv
# has one row per state per feature, and a model with more than one feature
# should widen this table rather than silently pair the wrong rows.
df.summary <- df.summary %>%
	dplyr::left_join(
		df.parameters %>% dplyr::select(
			State, feature, gaussian_mean, gaussian_var,
			empirical_mean, empirical_sd, n_windows, weight
		),
		by = "State"
	)

if (any(is.na(df.summary$empirical_mean) & df.summary$BP > 0)) {
	stop("a state with sequence in it has no row in ", parameters_csv,
	     " -- do the two files describe the same model?")
}

df.by_chrom <- dplyr::bind_rows(l.by_chrom)

# One "unsegmented" row per chromosome: blacklisted windows plus any contig the
# HMM never saw. Without it the stacked bars stop short of 100% with no
# explanation of where the rest went.
df.unsegmented <- df.by_chrom %>%
	dplyr::group_by(chr) %>%
	dplyr::summarize(Segmented = sum(BP), .groups = "drop") %>%
	dplyr::right_join(df.sizes, by = "chr") %>%
	tidyr::replace_na(list(Segmented = 0)) %>%
	dplyr::mutate(
		n_states = n_states,
		State    = "unsegmented",
		BP       = length - Segmented
	) %>%
	dplyr::select(chr, BP, n_states, State)

df.by_chrom <- dplyr::bind_rows(df.by_chrom, df.unsegmented) %>%
	dplyr::left_join(df.sizes, by = "chr") %>%
	dplyr::mutate(Percent = 100 * BP / length) %>%
	dplyr::rename(Chromosome_length = length) %>%
	dplyr::arrange(match(chr, df.sizes$chr), State)

dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df.summary, out_summary)
readr::write_csv(df.by_chrom, out_by_chrom)

message("k = ", n_states, ": ", format(segmented_bp, big.mark = ","),
        " bp segmented across ", nrow(df.summary), " states")
message("wrote ", out_summary, " and ", out_by_chrom)
