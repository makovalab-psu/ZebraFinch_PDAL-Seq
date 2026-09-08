# Split one fitted HMM into a BED file per hidden state.
#
# Adapted from workflow/scripts/07_PhyloHGMP_model/Annotate_state.R in the
# PDAL-Seq manuscript pipeline. Two differences:
#
#   * The manuscript reads the alignment BED and the state assignment as two
#     separate files and pairs them BY ROW POSITION. Phase 3's GHMM.py carries
#     chr/start/end through into GHMM_states.csv.gz, so the coordinates come
#     from the same row as the state and nothing depends on two files having
#     been written in the same order.
#
#   * Adjacent 1 kbp windows in the same state are merged into a segment. The
#     manuscript's states were alignment blocks and were already segment-like;
#     ours are a fixed grid, so without this a single state is ~10^5 separate
#     one-window records. Merging is what makes "the HMM state segments" in
#     the methods mean what it says, and it also removes the double counting
#     that `sum += $3 - $2 + 1` produces at every window boundary.
#
# Window coordinates follow the project convention set by make_windows.sh:
# 1-based inclusive, so window i+1 starts at end_i + 1 and a window's length
# is end - start + 1. The contiguity test below depends on that.
#
# A variational HMM is free to leave a state unoccupied. That is information,
# not an error, so an empty state still gets an (empty) BED file -- every
# downstream expand() over 0..k-1 stays satisfied and the enrichment script
# short circuits on it.
#
# Usage:
#   Rscript Annotate_states.R <GHMM_states.csv.gz> <n_states> <out_dir>

suppressPackageStartupMessages(library(tidyverse))

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Annotate_states.R <GHMM_states.csv.gz> <n_states> <out_dir>")
}

in_csv   <- args[1]
n_states <- as.integer(args[2])
out_dir  <- args[3]

if (is.na(n_states) || n_states < 1) {
	stop("n_states must be a positive integer, got '", args[2], "'")
}

df <- readr::read_csv(
	in_csv,
	col_types = readr::cols(
		chr = readr::col_character(),
		start = readr::col_double(),
		end = readr::col_double(),
		hidden_state = readr::col_integer(),
		.default = readr::col_double()
	)
)

for (column in c("chr", "start", "end", "hidden_state")) {
	if (!column %in% names(df)) {
		stop("column '", column, "' is missing from ", in_csv)
	}
}

n <- nrow(df)
if (n == 0) {
	stop("no windows in ", in_csv)
}

observed_states <- sort(unique(df$hidden_state))
if (max(observed_states) >= n_states) {
	stop("state ", max(observed_states), " appears in ", in_csv,
	     " but only ", n_states, " states were requested")
}

# A new segment starts wherever the contig changes, the state changes, or the
# window grid is broken (a blacklisted window was dropped between the two).
run_starts <- c(
	TRUE,
	df$chr[-1] != df$chr[-n] |
	df$hidden_state[-1] != df$hidden_state[-n] |
	df$start[-1] != df$end[-n] + 1
)

idx_start <- which(run_starts)
idx_end   <- c(idx_start[-1] - 1, n)

segments <- tibble::tibble(
	chr          = df$chr[idx_start],
	start        = df$start[idx_start],
	end          = df$end[idx_end],
	hidden_state = df$hidden_state[idx_start]
)

message(n, " windows collapsed into ", nrow(segments), " segments")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (state in seq_len(n_states) - 1L) {

	df.state <- segments %>% dplyr::filter(hidden_state == state)
	path <- file.path(out_dir, paste0("state_", state, ".bed"))

	# write.table on a zero-row frame makes an empty file, which is exactly
	# what an unoccupied state should leave behind.
	write.table(
		df.state[, c("chr", "start", "end")],
		path,
		sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
	)

	bp <- if (nrow(df.state) == 0) 0 else sum(df.state$end - df.state$start + 1)
	message(sprintf(
		"  state %-3d %8d segments  %10.0f bp  %5.2f%%",
		state, nrow(df.state), bp,
		100 * bp / sum(segments$end - segments$start + 1)
	))
}

if (length(observed_states) < n_states) {
	message("NOTE: ", n_states - length(observed_states),
	        " of ", n_states, " states are unoccupied and got empty BED files")
}
