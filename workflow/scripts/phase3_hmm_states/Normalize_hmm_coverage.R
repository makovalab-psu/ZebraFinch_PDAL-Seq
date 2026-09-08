# Build the feature matrix the Gaussian HMM is fit to.
#
# Adapted from
#   workflow/scripts/09_Human_HMM/Normalize_human_window_coverage.R
# in the PDAL-Seq manuscript pipeline. Two deliberate differences:
#
#   * The manuscript computes `SD = mean(x)` and then divides by it, so its
#     "z-score" is really x/mean - 1. We use sd(). With a single feature this
#     is only a scale factor -- it shifts every model's log likelihood by the
#     same N*log(scale), so the argmin of BIC is unchanged -- but the fit
#     plots are read on a z-scale, so it is fixed here.
#
#   * The manuscript drops the coordinates and re-attaches them later by row
#     position (Annotate_state.R). We carry chr/start/end through the whole
#     chain instead, so nothing downstream depends on two files having been
#     written in the same order.
#
# Usage:
#   Rscript Normalize_hmm_coverage.R <features> <transform> <out.csv>
#
#   <features>  comma separated experiment names, one filtered .bg each
#   <transform> "none" or "log1p", applied to the per-bp per-billion signal
#               before z-scoring

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Normalize_hmm_coverage.R <features> <transform> <out.csv>")
}

Features  <- strsplit(args[1], ",")[[1]]
Transform <- args[2]
out_csv   <- args[3]

if (!Transform %in% c("none", "log1p")) {
	stop("unknown transform '", Transform, "' -- expected none or log1p")
}

cov_dir <- "data/phase3/filtered_coverage"

coords <- NULL
values <- list()

for (feature in Features) {

	path <- file.path(cov_dir, paste0(feature, ".bg"))
	message("reading ", path)

	df <- readr::read_tsv(
		path,
		col_names = c("chr", "start", "end", "coverage"),
		col_types = readr::cols(
			chr = readr::col_character(),
			start = readr::col_double(),
			end = readr::col_double(),
			coverage = readr::col_double()
		)
	)

	if (nrow(df) == 0) {
		stop("no windows left in ", path,
		     " -- the blacklist filter removed everything")
	}

	# Every feature must describe the same windows in the same order, or the
	# columns of the matrix below would refer to different parts of the genome.
	if (is.null(coords)) {
		coords <- df[, c("chr", "start", "end")]
	} else if (!identical(coords$chr, df$chr) || !identical(coords$start, df$start)) {
		stop("window rows do not line up between ", Features[1], " and ", feature)
	}

	# Manuscript quantity: summed per-base coverage per bp of window, per
	# billion of the genome-wide total. Windows are a fixed width, so the
	# per-bp division only matters for the final partial window of a contig.
	total_billion <- sum(df$coverage, na.rm = TRUE) / 1e9
	if (!is.finite(total_billion) || total_billion <= 0) {
		stop("total coverage for ", feature, " is ", total_billion,
		     " -- the coverage file is empty or unreadable")
	}

	length_bp <- df$end - df$start + 1
	signal <- df$coverage / length_bp / total_billion

	if (Transform == "log1p") {
		signal <- log1p(signal)
	}

	MEAN <- mean(signal, na.rm = TRUE)
	SD   <- sd(signal, na.rm = TRUE)
	if (!is.finite(SD) || SD == 0) {
		stop("signal for ", feature, " has zero variance -- nothing to segment")
	}

	values[[feature]] <- (signal - MEAN) / SD

	message("  ", nrow(df), " windows, transform=", Transform,
	        ", mean=", signif(MEAN, 4), ", sd=", signif(SD, 4))
}

df.out <- dplyr::bind_cols(coords, as.data.frame(values, check.names = FALSE))

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df.out, out_csv)

message("wrote ", nrow(df.out), " windows x ", length(Features), " features to ", out_csv)
