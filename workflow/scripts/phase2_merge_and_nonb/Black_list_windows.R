# Blacklist unmappable 10 kbp windows from the WGS control.
#
# Adapted from
#   workflow/scripts/03_Correlation_and_black_list_unmappable_regions/black_list_windows.R
# in the PDAL-Seq manuscript pipeline. Changes: one genome instead of a
# seven-species loop, paths passed as arguments instead of hard coded, and
# the input is the window coverage Phase 1 already wrote rather than a
# freshly computed one.
#
# The rule: a window whose All_reads coverage exceeds its High_MapQ
# coverage by >= Threshold is losing that fraction of its signal to the
# MAPQ >= 20 filter, i.e. it is multi-mappable. The manuscript uses 1.1,
# so >10% lost. Ratios that are infinite (no high-MapQ coverage at all),
# undefined (no coverage at all), or extreme are all pinned to 2, which
# puts them above any sane threshold and so inside the blacklist.
#
# Note the manuscript calls the megadepth column "count"; --op sum makes it
# summed per-base coverage. The ratio of two such sums is unaffected.
#
# Usage:
#   Rscript Black_list_windows.R <All_reads.bed> <High_MapQ.bed> <threshold> <out.bed>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
	stop("Usage: Black_list_windows.R <All_reads.bed> <High_MapQ.bed> <threshold> <out.bed>")
}

all_file  <- args[1]
high_file <- args[2]
threshold <- as.numeric(args[3])
out_file  <- args[4]

read.bed <- function(x) {
	df <- read.delim(x, sep = "", header = FALSE)
	colnames(df) <- c("chr", "start", "end", "count")
	df
}

df.all  <- read.bed(all_file)
df.high <- read.bed(high_file)

if (!identical(df.all$chr, df.high$chr) || !identical(df.all$start, df.high$start)) {
	stop("window rows do not line up between the All_reads and High_MapQ coverage files")
}

df <- df.all[, c("chr", "start", "end")]
df$All       <- df.all$count
df$High_MapQ <- df.high$count
df$Ratio     <- df$All / df$High_MapQ

df$Ratio[which(is.infinite(df$Ratio))] <- 2
df$Ratio[which(is.na(df$Ratio))]       <- 2
df$Ratio[which(df$Ratio > 2)]          <- 2

df.bl <- df[which(df$Ratio >= threshold), ]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(df.bl[, c("chr", "start", "end")], out_file,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

message("blacklisted ", nrow(df.bl), " of ", nrow(df), " windows (",
        sprintf("%.1f", 100 * nrow(df.bl) / nrow(df)), "%) at threshold ", threshold)
message("median All/High_MapQ ratio: ", sprintf("%.3f", median(df$Ratio)))
