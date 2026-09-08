# Correlate merged-experiment read density against non-B DNA motif density
# in genome windows.
#
# Adapted from
#   workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/compile_reads_vs_nonb.R
# in the PDAL-Seq manuscript pipeline. Changes:
#
#   * one genome instead of a seven-species loop
#   * TWO motif metrics, not one. The manuscript correlates against the
#     motif COUNT from bedtools coverage. Six of the seven Zebra finch
#     classes self-overlap (see the merge tests in resources/README.md) and
#     the classes overlap each other, so "all" in particular double counts.
#     base_density -- bases of the window covered by >= 1 motif -- is
#     immune to that, so both are reported.
#   * TWO signals. CPM per experiment (the manuscript's quantity) and the
#     PDAL-Seq/WGS enrichment ratio, which is what "PDAL-Seq/Control read
#     density" asks for directly.
#   * a blacklist branch. The manuscript's step 04 does not apply its
#     blacklist here; both are reported so the choice stays open.
#
# Usage:
#   Rscript Compile_reads_vs_nonb.R <experiments> <nonb> <windows> <wgs> \
#           <exclude_chroms> <compare_window> <out_correlation.csv> <out_compare.csv>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
	stop("Usage: Compile_reads_vs_nonb.R <experiments> <nonb> <windows> <wgs> <exclude> <compare_window> <out_correlation.csv> <out_compare.csv>")
}

Experiments    <- strsplit(args[1], ",")[[1]]
Nonb           <- strsplit(args[2], ",")[[1]]
Windows        <- strsplit(args[3], ",")[[1]]
wgs            <- args[4]
Exclude        <- strsplit(args[5], ",")[[1]]
compare_window <- args[6]
out_corr       <- args[7]
out_compare    <- args[8]

cov_dir  <- "data/phase2/experiment_window_coverage/All_reads"
dens_dir <- "data/phase2/nonb_density"
flag_dir <- "data/phase2/blacklist/window_flag"

# Short, readable labels. "Tguttata_CFS414_PDALSeq" (all treated libraries
# merged) has no concentration suffix, so it becomes "combined".
label.of <- function(e) {
	if (e == wgs) return("WGS")
	suffix <- sub("^.*PDALSeq_?", "", e)
	if (suffix == "") "combined" else suffix
}

read.cov <- function(x) {
	df <- read.delim(x, sep = "", header = FALSE)
	colnames(df) <- c("chr", "start", "end", "coverage")
	df
}

read.dens <- function(x) {
	df <- read.delim(x, sep = "", header = FALSE)
	colnames(df) <- c("chr", "start", "end", "count", "bases", "length", "base_density")
	df
}

results <- list()
compare <- list()
k <- 0

for (window in Windows) {

	message("loading ", window, " nucleotide windows")

	df.flag <- read.delim(file.path(flag_dir, paste0(window, "_nucleotides.bed")),
	                      sep = "", header = FALSE)
	colnames(df.flag) <- c("chr", "start", "end", "blacklist_n")

	coords <- df.flag[, c("chr", "start", "end")]

	check.coords <- function(df, what) {
		if (!identical(coords$chr, df$chr) || !identical(coords$start, df$start)) {
			stop("window rows do not line up for ", what, " at ", window, " nucleotides")
		}
	}

	# Coverage per experiment, normalised genome wide BEFORE any window is
	# excluded, so CPM means the same thing in every branch below.
	cpm <- list()
	for (e in Experiments) {
		df.e <- read.cov(file.path(cov_dir, paste0(window, "_nucleotides"), paste0(e, ".bed")))
		check.coords(df.e, e)
		cpm[[e]] <- 1e6 * df.e$coverage / sum(df.e$coverage)
	}

	dens <- list()
	for (nb in Nonb) {
		df.n <- read.dens(file.path(dens_dir, paste0(window, "_nucleotides"), paste0(nb, ".bg")))
		check.coords(df.n, nb)
		dens[[nb]] <- df.n
	}

	keep_chr <- !(coords$chr %in% Exclude)
	not_blacklisted <- df.flag$blacklist_n == 0

	message("  ", sum(keep_chr), " of ", nrow(coords), " windows on included chromosomes; ",
	        sum(keep_chr & not_blacklisted), " of those outside the blacklist")

	for (e in Experiments) {

		signals <- list(CPM = cpm[[e]])
		if (e != wgs) {
			signals[["Ratio_vs_WGS"]] <- cpm[[e]] / cpm[[wgs]]
		}

		for (nb in Nonb) {
			for (metric in c("count", "base_density")) {

				motif <- dens[[nb]][[metric]]

				for (filtered in c(FALSE, TRUE)) {

					mask <- keep_chr
					if (filtered) mask <- mask & not_blacklisted

					for (sig in names(signals)) {

						x <- signals[[sig]][mask]
						y <- motif[mask]

						# Ratio_vs_WGS is not finite where the control has
						# no coverage in the window. Those windows are
						# dropped from the ratio branch only, and the count
						# that survived is recorded.
						ok <- is.finite(x) & is.finite(y)
						x <- x[ok]
						y <- y[ok]

						k <- k + 1
						if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
							results[[k]] <- data.frame(
								Experiment = e, Label = label.of(e),
								Library = if (e == wgs) "WGS" else "PDALSeq",
								nonB = nb, Window = window, Metric = metric,
								Signal = sig, Blacklist_filtered = filtered,
								N_windows = length(x),
								Pearson = NA_real_, Spearman = NA_real_,
								stringsAsFactors = FALSE
							)
							next
						}

						t.pearson  <- suppressWarnings(cor.test(x, y, method = "pearson"))
						t.spearman <- suppressWarnings(cor.test(x, y, method = "spearman"))

						results[[k]] <- data.frame(
							Experiment = e, Label = label.of(e),
							Library = if (e == wgs) "WGS" else "PDALSeq",
							nonB = nb, Window = window, Metric = metric,
							Signal = sig, Blacklist_filtered = filtered,
							N_windows = length(x),
							Pearson  = unname(t.pearson$estimate),
							Spearman = unname(t.spearman$estimate),
							stringsAsFactors = FALSE
						)
					}
				}
			}
		}
	}

	# The per-window table the scatter plot reads, written from the widest
	# windows so the panels are not a solid block of ink.
	if (window == compare_window) {
		for (e in Experiments) {
			for (nb in Nonb) {
				df.c <- coords[keep_chr, ]
				df.c$Window             <- window
				df.c$Experiment         <- e
				df.c$Label              <- label.of(e)
				df.c$CPM                <- cpm[[e]][keep_chr]
				df.c$Ratio_vs_WGS       <- if (e == wgs) NA_real_ else (cpm[[e]] / cpm[[wgs]])[keep_chr]
				df.c$nonB               <- nb
				df.c$motif_count        <- dens[[nb]]$count[keep_chr]
				df.c$motif_base_density <- dens[[nb]]$base_density[keep_chr]
				df.c$Blacklisted        <- !not_blacklisted[keep_chr]
				compare[[length(compare) + 1]] <- df.c
			}
		}
	}
}

df.corr <- bind_rows(results)
dir.create(dirname(out_corr), recursive = TRUE, showWarnings = FALSE)
write.csv(df.corr, out_corr, row.names = FALSE)
message("wrote ", out_corr, " (", nrow(df.corr), " rows)")

df.compare <- bind_rows(compare)
write.csv(df.compare, out_compare, row.names = FALSE)
message("wrote ", out_compare, " (", nrow(df.compare), " rows at ",
        compare_window, " nucleotide windows)")
