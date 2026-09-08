# The union of everything the HMM segmented: genome minus the Phase 2
# blacklist, minus the sex chromosomes, chrMT and the three rDNA morphs.
#
# Written once, from the HMM input matrix rather than from any one model, so
# it is identical for every k. Used as the denominator when Phase 4 reports
# what fraction of the segmentable genome each state covers.
#
# It is NOT used to constrain `bedtools shuffle`. The enrichment null shuffles
# across the whole genome, as the manuscript's Enrichment.sh does. That means
# the null includes placements in regions no state covers, which dilutes it;
# the enrichment ratios are therefore relative to a whole-genome expectation,
# not to a within-segmentation expectation. See README section 06.
#
# Usage:
#   Rscript Segmented_genome.R <HMM_input.csv> <out.bed>

suppressPackageStartupMessages(library(tidyverse))

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Segmented_genome.R <HMM_input.csv> <out.bed>")
}

df <- readr::read_csv(
	args[1],
	col_types = readr::cols(
		chr = readr::col_character(),
		start = readr::col_double(),
		end = readr::col_double(),
		.default = readr::col_double()
	)
)

n <- nrow(df)
if (n == 0) {
	stop("no windows in ", args[1])
}

# Same contiguity rule as Annotate_states.R: 1-based inclusive windows, so a
# gap is any row whose start is not the previous row's end + 1.
run_starts <- c(TRUE, df$chr[-1] != df$chr[-n] | df$start[-1] != df$end[-n] + 1)
idx_start  <- which(run_starts)
idx_end    <- c(idx_start[-1] - 1, n)

segments <- tibble::tibble(
	chr   = df$chr[idx_start],
	start = df$start[idx_start],
	end   = df$end[idx_end]
)

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.table(segments, args[2], sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message(nrow(segments), " segments covering ",
        format(sum(segments$end - segments$start + 1), big.mark = ","),
        " bp on ", length(unique(segments$chr)), " contigs")
