# Compile the Phase 5 enrichment sweep (sample x repeat class) into two tables.
#
# Table 1, Repeat_enrichment.csv, is one row per sample x repeat class:
# enrichment is mean(Observed) / mean(Null) over the iterations, exactly as
# the manuscript's Compile_CenSat_Class_enrichment.R computes it.
#
# Table 2, PDALSeq_versus_control.csv, is the comparison the phase exists for.
# A repeat class can look enriched in PDAL-Seq for two reasons: because
# permanganate-reactive ssDNA is there, or because the class is a mapping
# artefact -- multi-mapping reads that bwa mem scattered into it, or a
# sequencing bias in a GC-extreme array. The PCR-free control has the second
# and not the first, so
#
#     Ratio = Enrichment(PDAL-Seq sample) / Enrichment(PCR-free control)
#
# is what survives after the artefact is divided out. Read table 2, not table
# 1, for the biology; table 1 is what it is built from.
#
# Significance follows the manuscript's rule, with the threshold set to the
# subsample size rather than hardcoded at 10,000. In the manuscript those two
# numbers were the same (both 10,000), and the reason for the rule is a
# property of the subsampling rather than of the number itself: at or below
# the subsample size, `shuf -n` returns the whole file every iteration, so
# Observed is a constant and a Wilcoxon test would be comparing a point mass
# against a distribution. Phase 5 subsamples 1,000, so 1,000 is the threshold.
#
# That case is the norm here, not the exception -- most satellite families
# have a handful of annotations -- so the empirical p-value is doing most of
# the work. Both p-values are reported for every row either way.
#
# Usage:
#   Rscript Compile_repeat_enrichment.R <in_dir> <sample_table> <te_classes> \
#           <satellite_classes> <other_classes> <subsample_size> \
#           <control_sample> <out_enrichment.csv> <out_versus_control.csv>
#
#   <in_dir> holds <sample>/<repeat>.csv

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
	stop("Usage: Compile_repeat_enrichment.R <in_dir> <sample_table> <te_classes> <satellite_classes> <other_classes> <subsample_size> <control_sample> <out_enrichment.csv> <out_versus_control.csv>")
}

in_dir          <- args[1]
sample_table    <- args[2]
te_classes      <- strsplit(args[3], ",")[[1]]
sat_classes     <- strsplit(args[4], ",")[[1]]
other_classes   <- strsplit(args[5], ",")[[1]]
subsample_size  <- as.numeric(args[6])
control_sample  <- args[7]
out_enrichment  <- args[8]
out_versus      <- args[9]

repeat_classes <- c(te_classes, sat_classes, other_classes)

# Source is what the class was cut out of, and it is also how the plots are
# faceted. "all_TE" and "all_Satellite" are every record of their annotation
# pooled; they are the one-line summary of their panel and are kept apart from
# it so they cannot be read as just another family.
source_of <- function(x) {
	dplyr::case_when(
		x == "all_TE"        ~ "All TEs",
		x == "all_Satellite" ~ "All satellites",
		x == "CEN"           ~ "Centromere",
		x %in% te_classes    ~ "Transposable element",
		x %in% sat_classes   ~ "Satellite",
		TRUE                 ~ "Other"
	)
}

df.samples <- readr::read_tsv(sample_table, col_types = readr::cols(.default = "c"))

# Labels for the figures: the PCR-free library is "WGS", and a PDAL-Seq
# library is its permanganate concentration, with the replicate appended only
# when the concentration has more than one. Same convention as Phase 1 and 2.
df.samples <- df.samples %>%
	dplyr::group_by(Treatment) %>%
	dplyr::mutate(n_replicates = dplyr::n()) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(Label = dplyr::if_else(
		Library == "PCRfree",
		"WGS",
		dplyr::if_else(n_replicates > 1,
		               paste0(Treatment, "-", Replicate),
		               Treatment)
	))

samples <- df.samples$Sequencing_sample

if (!(control_sample %in% samples)) {
	stop("control sample ", control_sample, " is not in ", sample_table)
}

grid <- expand.grid(
	Sample = samples,
	Repeat = repeat_classes,
	stringsAsFactors = FALSE
)

Mean.observed <- numeric(nrow(grid))
Mean.null     <- numeric(nrow(grid))
Wilcox.p      <- numeric(nrow(grid))
Empirical.p   <- numeric(nrow(grid))
N.annotation  <- numeric(nrow(grid))
Iterations    <- integer(nrow(grid))

for (i in seq_len(nrow(grid))) {

	path <- file.path(in_dir, grid$Sample[i], paste0(grid$Repeat[i], ".csv"))

	if (!file.exists(path)) {
		stop("missing enrichment file: ", path)
	}

	df.i <- read.csv(path)

	Mean.observed[i] <- mean(df.i$Observed)
	Mean.null[i]     <- mean(df.i$Null)
	N.annotation[i]  <- df.i$N_annotation[1]
	Iterations[i]    <- nrow(df.i)

	# wilcox.test throws when both samples are constant, which happens
	# whenever a class is small enough that the subsample is the whole file
	# AND the null happens not to vary. NA is the honest answer there.
	Wilcox.p[i] <- tryCatch(
		wilcox.test(df.i$Observed, df.i$Null)$p.value,
		error = function(e) NA_real_,
		warning = function(w) suppressWarnings(
			wilcox.test(df.i$Observed, df.i$Null)$p.value
		)
	)

	# How often did the null land at least as far from its own mean as the
	# observed mean did? +1 in both places is the standard finite-sampling
	# correction, so p is never reported as exactly 0.
	delta_observed <- abs(Mean.observed[i] - Mean.null[i])
	Empirical.p[i] <- (sum(abs(df.i$Null - Mean.null[i]) >= delta_observed) + 1) /
	                  (nrow(df.i) + 1)
}

df <- tibble::tibble(
	Sample        = grid$Sample,
	Repeat        = grid$Repeat,
	Mean.observed = Mean.observed,
	Mean.null     = Mean.null,
	N_annotation  = N.annotation,
	Iterations    = Iterations,
	Wilcox.p      = Wilcox.p,
	Empirical.p   = Empirical.p
)

df$Enrichment     <- df$Mean.observed / df$Mean.null
df$log2Enrichment <- log2(df$Enrichment)

df$Test    <- ifelse(df$N_annotation > subsample_size, "wilcoxon", "empirical")
df$P_value <- ifelse(df$Test == "wilcoxon", df$Wilcox.p, df$Empirical.p)

df$Source <- source_of(df$Repeat)

df <- df %>%
	dplyr::left_join(
		df.samples %>% dplyr::select(Sequencing_sample, Experiment, Library,
		                             Treatment, Replicate, Label),
		by = c("Sample" = "Sequencing_sample")
	) %>%
	dplyr::select(Sample, Label, Experiment, Library, Treatment, Replicate,
	              Source, Repeat, dplyr::everything())

dir.create(dirname(out_enrichment), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df, out_enrichment)

message("wrote ", nrow(df), " rows (", length(samples), " samples x ",
        length(repeat_classes), " repeat classes) to ", out_enrichment)

n_empirical <- sum(df$Test == "empirical")
if (n_empirical > 0) {
	message("  ", n_empirical, " rows have <= ", subsample_size,
	        " annotations and use the empirical p-value")
}
if (any(!is.finite(df$Enrichment))) {
	message("  ", sum(!is.finite(df$Enrichment)),
	        " rows have a non-finite enrichment (a null summed to zero)")
}

####################################################################
# PDAL-Seq versus the PCR-free control
####################################################################

df.control <- df %>%
	dplyr::filter(Sample == control_sample) %>%
	dplyr::select(Repeat,
	              Control.enrichment = Enrichment,
	              Control.P_value    = P_value,
	              Control.N          = N_annotation)

df.versus <- df %>%
	dplyr::filter(Sample != control_sample) %>%
	dplyr::left_join(df.control, by = "Repeat") %>%
	dplyr::mutate(
		Ratio     = Enrichment / Control.enrichment,
		log2Ratio = log2(Ratio)
	) %>%
	dplyr::select(Sample, Label, Experiment, Treatment, Replicate,
	              Source, Repeat, N_annotation,
	              PDALSeq.enrichment = Enrichment,
	              Control.enrichment,
	              Ratio, log2Ratio,
	              PDALSeq.P_value = P_value,
	              Control.P_value,
	              Test)

readr::write_csv(df.versus, out_versus)

message("wrote ", nrow(df.versus), " rows to ", out_versus,
        " (control = ", control_sample, ")")
