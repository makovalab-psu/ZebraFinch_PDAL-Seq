# Bin the per-state CpG methylation histograms into rate ranges.
#
# Adapted from workflow/scripts/10_Methylation_analysis/
# compile_methylation_in_state.R in the PDAL-Seq manuscript pipeline, which
# bins into five ranges and reports the percent of CpGs in each.
#
# One scale difference: the manuscript's methylation bedGraph held fractions
# (0-1) and its bins were 0-0.2, 0.2-0.4, ... The Zebra finch file
# (bTaeGut7v0.4_MT_rDNA.matZ.PBmethylation.v0.1.bed) holds percentages, so the
# bins here are 0-20, 20-40, ... The Range labels are written on the 0-100
# scale so nothing silently looks like the manuscript's output when it is not.
#
# Usage:
#   Rscript Compile_methylation.R <in_dir> <n_states> <out.csv>
#
#   <in_dir> holds state_<i>.txt, each `count rate` as written by
#            methylation_in_state.sh

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Compile_methylation.R <in_dir> <n_states> <out.csv>")
}

in_dir   <- args[1]
n_states <- as.integer(args[2])
out_csv  <- args[3]

Start <- c(0, 20, 40, 60, 80)
End   <- c(20, 40, 60, 80, 100)

l.df <- list()

for (state in seq_len(n_states) - 1L) {

	path <- file.path(in_dir, paste0("state_", state, ".txt"))
	if (!file.exists(path)) {
		stop("missing methylation histogram: ", path)
	}

	if (file.size(path) == 0) {
		# Unoccupied state: no CpGs, so every bin is zero and the percent is
		# undefined rather than zero.
		message("state ", state, " has no CpGs (unoccupied)")
		l.df[[length(l.df) + 1]] <- tibble::tibble(
			n_states = n_states,
			State = as.character(state),
			CpG_5mC_rate_start = Start,
			CpG_5mC_rate_end = End,
			Count = 0,
			Total = 0,
			Percent = NA_real_
		)
		next
	}

	df.i <- read.delim(path, header = FALSE, sep = "",
	                   col.names = c("Count", "CpG_rate"))

	Count <- numeric(length(Start))
	for (j in seq_along(Start)) {
		# The top bin is closed on both sides so a fully methylated CpG at
		# exactly 100 is not dropped.
		in_bin <- if (End[j] == 100) {
			df.i$CpG_rate >= Start[j] & df.i$CpG_rate <= End[j]
		} else {
			df.i$CpG_rate >= Start[j] & df.i$CpG_rate < End[j]
		}
		Count[j] <- sum(df.i$Count[in_bin])
	}

	Total <- sum(Count)
	if (Total != sum(df.i$Count)) {
		stop("state ", state, ": ", sum(df.i$Count) - Total,
		     " CpGs fell outside 0-100 -- is this file on the 0-1 scale?")
	}

	l.df[[length(l.df) + 1]] <- tibble::tibble(
		n_states = n_states,
		State = as.character(state),
		CpG_5mC_rate_start = Start,
		CpG_5mC_rate_end = End,
		Count = Count,
		Total = Total,
		Percent = 100 * Count / Total
	)
}

df.out <- dplyr::bind_rows(l.df)
df.out$Range <- paste(df.out$CpG_5mC_rate_start, df.out$CpG_5mC_rate_end, sep = "_")

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df.out, out_csv)

message("wrote ", nrow(df.out), " rows to ", out_csv)
