# Correlate merged-experiment read density against non-B DNA motif density
# in genome windows, over EVERY window of the genome.
#
# This is the "including repetitive" version of
#   workflow/scripts/phase2_merge_and_nonb/Compile_reads_vs_nonb.R
# and it differs from it in exactly one respect: nothing is filtered out
# before the correlation. Regions are blacklisted BECAUSE they are
# repetitive, bwa mem assigns their multi-mapping reads to one of the
# alternatives at random rather than discarding them, and those repeats
# are what the comparison is about.
#
# The manuscript agrees. 04_.../compile_reads_vs_nonb.R never opens
# blacklist_files/; its only filter is !chr %in% c("chrX", "chrY"). The
# blacklist first appears at step 09, the HMM.
#
# Rather than swap one hard-coded filter for another, both filters are
# emitted as COLUMNS, so every reading of "all regions" is answerable
# from one table:
#
#   Region_set       All_windows      every window -- the headline
#                    Blacklist_free   blacklisted fraction == 0, i.e.
#                                     exactly Phase 2's filtered branch
#                    Blacklist_light  fraction <= <light>
#                    Blacklist_heavy  fraction >= <heavy>
#
#   Chromosome_set   All_contigs      everything, sex chromosomes and
#                                     chrMT and the rDNA morphs included
#                    Phase2_contigs   EXCLUDE_CHROMOSOMES removed
#
# All_windows x Phase2_contigs reproduces the manuscript's step 04
# convention exactly, so the Zebra finch numbers stay comparable with the
# published human ones.
#
# A FRACTION, not Phase 2's boolean. merge_blacklist runs `bedtools merge
# -d 1`, so blacklist intervals are long (9500 intervals of mean 65 kbp in
# the manuscript's human blacklist, 14% of autosomal sequence). A boolean
# "any overlap" test therefore drops 93% of 1 Mbp windows while dropping
# only 14% of the sequence. The fraction degrades gracefully where the
# boolean does not.
#
# CPM is normalised genome wide over ALL windows before any subsetting,
# exactly as Phase 2 does, so CPM means the same number in both phases.
#
# Usage:
#   Rscript Compile_reads_vs_nonb_all_regions.R <experiments> <nonb> <windows> \
#           <wgs> <exclude_chroms> <compare_window> <light> <heavy> \
#           <out_correlation.csv> <out_compare.csv> <out_composition.csv>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
	stop("Usage: Compile_reads_vs_nonb_all_regions.R <experiments> <nonb> <windows> <wgs> <exclude> <compare_window> <light> <heavy> <out_correlation.csv> <out_compare.csv> <out_composition.csv>")
}

Experiments    <- strsplit(args[1], ",")[[1]]
Nonb           <- strsplit(args[2], ",")[[1]]
Windows        <- strsplit(args[3], ",")[[1]]
wgs            <- args[4]
Exclude        <- strsplit(args[5], ",")[[1]]
compare_window <- args[6]
light          <- as.numeric(args[7])
heavy          <- as.numeric(args[8])
out_corr       <- args[9]
out_compare    <- args[10]
out_comp       <- args[11]

if (!is.finite(light) || !is.finite(heavy)) {
	stop("light and heavy must be numeric fractions, got '", args[7], "' and '", args[8], "'")
}

cov_dir  <- "data/phase2/experiment_window_coverage/All_reads"
dens_dir <- "data/phase2/nonb_density"
frac_dir <- "data/phase2_including_repetative/blacklist_fraction"

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

# bedtools coverage -a windows -b blacklist
read.frac <- function(x) {
	df <- read.delim(x, sep = "", header = FALSE)
	colnames(df) <- c("chr", "start", "end", "bl_count", "bl_bases", "length", "bl_fraction")
	df
}

results <- list()
compare <- list()
composition <- list()
k <- 0

for (window in Windows) {

	message("loading ", window, " nucleotide windows")

	df.frac <- read.frac(file.path(frac_dir, paste0(window, "_nucleotides.bed")))
	coords  <- df.frac[, c("chr", "start", "end")]
	n_win   <- nrow(coords)

	check.coords <- function(df, what) {
		if (!identical(coords$chr, df$chr) || !identical(coords$start, df$start)) {
			stop("window rows do not line up for ", what, " at ", window, " nucleotides")
		}
	}

	# Coverage per experiment, normalised genome wide over every window.
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

	# The two filter axes, precomputed once per window size rather than
	# rebuilt inside the innermost loop.
	region_masks <- list(
		All_windows     = rep(TRUE, n_win),
		Blacklist_free  = df.frac$bl_fraction == 0,
		Blacklist_light = df.frac$bl_fraction <= light,
		Blacklist_heavy = df.frac$bl_fraction >= heavy
	)

	excluded_contig <- coords$chr %in% Exclude
	chrom_masks <- list(
		All_contigs    = rep(TRUE, n_win),
		Phase2_contigs = !excluded_contig
	)

	window_bp <- df.frac$end - df.frac$start

	for (cs in names(chrom_masks)) {
		denom_n  <- sum(chrom_masks[[cs]])
		denom_bp <- sum(window_bp[chrom_masks[[cs]]])
		for (rs in names(region_masks)) {
			m <- chrom_masks[[cs]] & region_masks[[rs]]
			composition[[length(composition) + 1]] <- data.frame(
				Window = window, Chromosome_set = cs, Region_set = rs,
				N_windows = sum(m),
				Mbp = sum(window_bp[m]) / 1e6,
				Pct_windows = if (denom_n  > 0) 100 * sum(m) / denom_n else NA_real_,
				Pct_Mbp     = if (denom_bp > 0) 100 * sum(window_bp[m]) / denom_bp else NA_real_,
				stringsAsFactors = FALSE
			)
		}
	}

	message("  ", n_win, " windows; ",
	        sum(region_masks$Blacklist_free), " with no blacklisted base (",
	        sprintf("%.1f", 100 * sum(region_masks$Blacklist_free) / n_win), "%), ",
	        sum(region_masks$Blacklist_heavy), " at least ",
	        sprintf("%.0f", 100 * heavy), "% blacklisted")

	for (e in Experiments) {

		signals <- list(CPM = cpm[[e]])
		if (e != wgs) {
			signals[["Ratio_vs_WGS"]] <- cpm[[e]] / cpm[[wgs]]
		}

		for (nb in Nonb) {
			for (metric in c("count", "base_density")) {

				motif <- dens[[nb]][[metric]]

				for (cs in names(chrom_masks)) {
					for (rs in names(region_masks)) {

						mask <- chrom_masks[[cs]] & region_masks[[rs]]

						for (sig in names(signals)) {

							x <- signals[[sig]][mask]
							y <- motif[mask]

							# Ratio_vs_WGS is not finite where the control
							# has no coverage in the window. Those windows
							# are dropped from the ratio branch only, and
							# the count that survived is recorded.
							ok <- is.finite(x) & is.finite(y)
							x <- x[ok]
							y <- y[ok]

							k <- k + 1
							if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
								results[[k]] <- data.frame(
									Experiment = e, Label = label.of(e),
									Library = if (e == wgs) "WGS" else "PDALSeq",
									nonB = nb, Window = window, Metric = metric,
									Signal = sig, Chromosome_set = cs, Region_set = rs,
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
								Signal = sig, Chromosome_set = cs, Region_set = rs,
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
	}

	# The per-window table the scatter plot reads, written from the widest
	# windows so the panels are not a solid block of ink. Unlike Phase 2
	# this keeps EVERY window and carries the blacklisted fraction and the
	# contig class along, so the plot can show what the filter would have
	# removed instead of silently removing it.
	if (window == compare_window) {
		for (e in Experiments) {
			for (nb in Nonb) {
				df.c <- coords
				df.c$Window                 <- window
				df.c$Experiment             <- e
				df.c$Label                  <- label.of(e)
				df.c$CPM                    <- cpm[[e]]
				df.c$Ratio_vs_WGS           <- if (e == wgs) NA_real_ else cpm[[e]] / cpm[[wgs]]
				df.c$nonB                   <- nb
				df.c$motif_count            <- dens[[nb]]$count
				df.c$motif_base_density     <- dens[[nb]]$base_density
				df.c$blacklist_fraction     <- df.frac$bl_fraction
				df.c$Phase2_excluded_contig <- excluded_contig
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
        compare_window, " nucleotide windows, every window kept)")

df.comp <- bind_rows(composition)
write.csv(df.comp, out_comp, row.names = FALSE)
message("wrote ", out_comp, " (", nrow(df.comp), " rows)")

# The headline number: how much of the genome Phase 2 was discarding.
for (window in Windows) {
	sub <- df.comp[df.comp$Window == window & df.comp$Chromosome_set == "Phase2_contigs", ]
	all_n  <- sub$N_windows[sub$Region_set == "All_windows"]
	free_n <- sub$N_windows[sub$Region_set == "Blacklist_free"]
	message("at ", window, " nt, Phase 2 kept ", free_n, " of ", all_n, " windows (",
	        sprintf("%.1f", 100 * free_n / all_n), "%)")
}
