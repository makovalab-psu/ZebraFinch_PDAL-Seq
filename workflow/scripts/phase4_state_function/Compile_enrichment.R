# Compile one enrichment sweep (state x annotation) into a single table.
#
# Replaces three near-identical manuscript scripts with one, because they
# differ only in the directory they read and the name of the category column:
#
#   07_PhyloHGMP_model/compile_nonB_enrichment_in_state.R
#   07_PhyloHGMP_model/compile_functional_element_enrichment_in_state.R
#   07_PhyloHGMP_model/compile_RNA_expression_in_state.R
#
# Enrichment is mean(Observed) / mean(Null) over the 100 iterations, exactly as
# the manuscript computes it.
#
# Significance is where this goes further than the manuscript's compile
# scripts, which call wilcox.test unconditionally. The methods say:
#
#   "Significance was determined using a Wilcoxon test in R, if the number of
#    annotations in the original query was >10,000. Alternatively, for
#    annotations with a smaller number of loci, significance was determined by
#    determining how many times the null distribution (shuffling) produced a
#    value that differed from the null mean at least as much as the value
#    computed on the original annotations."
#
# Below 10,000 loci every iteration draws the SAME annotations -- `shuf -n` on
# a file with fewer records returns all of them -- so Observed is constant and
# a Wilcoxon test is comparing a point mass against a distribution. That is
# what the empirical test is for. The centromere annotation has 41 records, so
# this is not hypothetical here.
#
# Both p-values are reported for every row. P_value carries whichever one the
# methods select, and Test says which that was, so nothing has to be
# recomputed to check the other.
#
# Usage:
#   Rscript Compile_enrichment.R <in_dir> <n_states> <categories> \
#           <category_column> <out.csv>
#
#   <in_dir>           holds state_<i>/<category>.csv
#   <categories>       comma separated
#   <category_column>  what to call the category column in the output

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
	stop("Usage: Compile_enrichment.R <in_dir> <n_states> <categories> <category_column> <out.csv>")
}

in_dir       <- args[1]
n_states     <- as.integer(args[2])
categories   <- strsplit(args[3], ",")[[1]]
category_col <- args[4]
out_csv      <- args[5]

WILCOXON_MINIMUM <- 10000

grid <- expand.grid(
	State = as.character(seq_len(n_states) - 1L),
	Category = categories,
	stringsAsFactors = FALSE
)

Mean.observed <- numeric(nrow(grid))
Mean.null     <- numeric(nrow(grid))
Wilcox.p      <- numeric(nrow(grid))
Empirical.p   <- numeric(nrow(grid))
N.annotation  <- numeric(nrow(grid))
Iterations    <- integer(nrow(grid))

for (i in seq_len(nrow(grid))) {

	path <- file.path(in_dir,
	                  paste0("state_", grid$State[i]),
	                  paste0(grid$Category[i], ".csv"))

	if (!file.exists(path)) {
		stop("missing enrichment file: ", path)
	}

	df.i <- read.csv(path)

	Mean.observed[i] <- mean(df.i$Observed)
	Mean.null[i]     <- mean(df.i$Null)
	N.annotation[i]  <- df.i$N_annotation[1]
	Iterations[i]    <- nrow(df.i)

	# wilcox.test throws when both samples are constant, which is exactly the
	# unoccupied-state case (all zeros). NA is the honest answer there.
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
	n_states      = n_states,
	State         = grid$State,
	Category      = grid$Category,
	Mean.observed = Mean.observed,
	Mean.null     = Mean.null,
	N_annotation  = N.annotation,
	Iterations    = Iterations,
	Wilcox.p      = Wilcox.p,
	Empirical.p   = Empirical.p
)

df$Enrichment     <- df$Mean.observed / df$Mean.null
df$log2Enrichment <- log2(df$Enrichment)

df$Test    <- ifelse(df$N_annotation > WILCOXON_MINIMUM, "wilcoxon", "empirical")
df$P_value <- ifelse(df$Test == "wilcoxon", df$Wilcox.p, df$Empirical.p)

names(df)[names(df) == "Category"] <- category_col

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df, out_csv)

message("wrote ", nrow(df), " rows (", n_states, " states x ",
        length(categories), " categories) to ", out_csv)

n_empirical <- sum(df$Test == "empirical")
if (n_empirical > 0) {
	message("  ", n_empirical, " rows fell below ", WILCOXON_MINIMUM,
	        " annotations and use the empirical p-value")
}
if (any(!is.finite(df$Enrichment))) {
	message("  ", sum(!is.finite(df$Enrichment)),
	        " rows have a non-finite enrichment (a state or a null was empty)")
}
