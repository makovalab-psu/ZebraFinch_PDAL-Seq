# k-means scree data: total within-cluster sum of squares against k.
#
# Adapted from the scree half of
#   workflow/scripts/09_Human_HMM/Human_Kmeans_scree.R
# in the PDAL-Seq manuscript pipeline. Differences:
#
#   * the sample-vs-sample correlation heatmap that script also builds is
#     dropped -- it is meaningless with one feature, and Phase 1 already
#     answered the question it was asking;
#   * iter.max is raised from R's default of 10. At k = 30 on ~10^6 points the
#     default stops before convergence and warns, which silently inflates wss
#     for exactly the models at the right-hand end of the scree;
#   * plotting is a separate rule, so a change to the figure does not re-run
#     twelve k-means fits.
#
# Usage:
#   Rscript Kmeans_scree.R <input.csv> <features> <states> <out.csv>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
	stop("Usage: Kmeans_scree.R <input.csv> <features> <states> <out.csv>")
}

in_csv   <- args[1]
Features <- strsplit(args[2], ",")[[1]]
States   <- as.integer(strsplit(args[3], ",")[[1]])
out_csv  <- args[4]

df <- readr::read_csv(in_csv, show_col_types = FALSE)

missing <- setdiff(Features, colnames(df))
if (length(missing) > 0) {
	stop("columns missing from ", in_csv, ": ", paste(missing, collapse = ", "))
}

mat <- as.matrix(df[, Features, drop = FALSE])

message(nrow(mat), " windows x ", ncol(mat), " features")

set.seed(123)

results <- list()

for (i in seq_along(States)) {

	k <- States[i]
	message("k = ", k)

	km <- kmeans(mat, centers = k, nstart = 20, iter.max = 50)

	# ifault is only returned by the default Hartigan-Wong algorithm: 2 means
	# iter.max was hit, 4 means the quick-transfer stage was. Either way the
	# wss for this k is an upper bound rather than the converged value.
	ifault <- if (is.null(km$ifault)) 0L else as.integer(km$ifault)
	converged <- !(ifault %in% c(2L, 4L))
	if (!converged) {
		warning("k-means did not converge at k = ", k, " (ifault ", ifault, ")")
	}

	results[[i]] <- data.frame(
		n_clusters      = k,
		wss             = km$tot.withinss,
		betweenss       = km$betweenss,
		totss           = km$totss,
		variance_explained = km$betweenss / km$totss,
		iterations      = km$iter,
		converged       = converged
	)

	message("  wss = ", signif(km$tot.withinss, 6),
	        ", variance explained = ", signif(km$betweenss / km$totss, 4))
}

df.out <- dplyr::bind_rows(results)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df.out, out_csv)

message("wrote ", out_csv)
