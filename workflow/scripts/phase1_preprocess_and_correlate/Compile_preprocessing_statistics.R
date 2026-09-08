# Compile one row of preprocessing statistics per sequencing sample.
#
# Adapted from
# workflow/scripts/02_preprocessing_statistics/Compile_preprocessing_statistics.R
# in the PDAL-Seq manuscript pipeline. Restricted to R1 samples, and the
# file paths point at data/phase1/statistics/ instead of the manuscript's
# data/02_preprocessing_statistics/.
#
# Usage: Rscript Compile_preprocessing_statistics.R <samples.txt> <out.csv>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Compile_preprocessing_statistics.R <samples.txt> <out.csv>")
}

samples_file <- args[1]
out_file <- args[2]

stat_dir <- "data/phase1/statistics"

df.index <- read.delim(samples_file, stringsAsFactors = FALSE)
samples <- df.index$Sequencing_sample

# ---------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------

# A single integer written by read_count_fastq.gz.sh.
read.count <- function(path) {
	as.numeric(readLines(path)[1])
}

# `uniq -c` output: leading whitespace, count, value. sep = "" splits on
# any run of whitespace, which is how the manuscript reads these files.
read.uniq.c <- function(path, value_name) {
	df <- read.delim(path, sep = "", header = FALSE,
	                 col.names = c("count", value_name),
	                 stringsAsFactors = FALSE)
	df$count <- as.numeric(df$count)
	df[[value_name]] <- as.numeric(df[[value_name]])
	df
}

# Estimate_unique_maps.sh prints a header line before the uniq -c block.
read.unique.maps <- function(path) {
	lines <- readLines(path)
	lines <- lines[-1]
	lines <- lines[trimws(lines) != ""]
	parts <- strsplit(trimws(lines), "[[:space:]]+")
	data.frame(
		count = as.numeric(sapply(parts, `[`, 1)),
		score_difference = as.numeric(sapply(parts, `[`, 2)),
		stringsAsFactors = FALSE
	)
}

# ---------------------------------------------------------------------
# Flag arithmetic
#
# These are single-end R1 libraries, so every trimmed read contributes
# exactly one primary record, plus optional secondary (256) and
# supplementary (2048) records.
# ---------------------------------------------------------------------

flag.summary <- function(df) {
	unmapped      <- bitwAnd(df$flag, 4)    > 0
	secondary     <- bitwAnd(df$flag, 256)  > 0
	supplementary <- bitwAnd(df$flag, 2048) > 0
	primary       <- !secondary & !supplementary

	list(
		total         = sum(df$count),
		primary       = sum(df$count[primary]),
		mapped        = sum(df$count[primary & !unmapped]),
		unmapped      = sum(df$count[primary & unmapped]),
		secondary     = sum(df$count[secondary]),
		supplementary = sum(df$count[supplementary])
	)
}

# ---------------------------------------------------------------------
# Per-sample compilation
# ---------------------------------------------------------------------

rows <- list()

for (i in seq_along(samples)) {
	s <- samples[i]
	message("compiling ", s)

	raw     <- read.count(file.path(stat_dir, "Raw_fastq", paste0(s, ".txt")))
	trimmed <- read.count(file.path(stat_dir, "Trimmed_fastq", paste0(s, ".txt")))

	f.mapped  <- flag.summary(read.uniq.c(file.path(stat_dir, "bam_flags", "Mapped_bam", paste0(s, ".txt")), "flag"))
	f.dedup   <- flag.summary(read.uniq.c(file.path(stat_dir, "bam_flags", "All_reads", paste0(s, ".txt")), "flag"))
	f.highmq  <- flag.summary(read.uniq.c(file.path(stat_dir, "bam_flags", "High_MapQ", paste0(s, ".txt")), "flag"))

	df.mapq <- read.uniq.c(file.path(stat_dir, "mapping_quality", paste0(s, ".txt")), "mapq")
	# Median MAPQ over the deduplicated alignments, from the count histogram.
	mapq.expanded.median <- {
		o <- order(df.mapq$mapq)
		v <- df.mapq$mapq[o]
		n <- df.mapq$count[o]
		cum <- cumsum(n)
		v[which(cum >= sum(n) / 2)[1]]
	}

	df.unique <- read.unique.maps(file.path(stat_dir, "unique_maps", paste0(s, ".txt")))
	unique.total <- sum(df.unique$count)
	unique.best  <- sum(df.unique$count[df.unique$score_difference > 0])
	unique.tied  <- sum(df.unique$count[df.unique$score_difference == 0])

	rows[[i]] <- data.frame(
		Sequencing_sample                 = s,
		Raw_reads                         = raw,
		Trimmed_reads                     = trimmed,
		Percent_surviving_trimming        = 100 * trimmed / raw,
		Total_alignments                  = f.mapped$total,
		Primary_alignments                = f.mapped$primary,
		Mapped_reads                      = f.mapped$mapped,
		Unmapped_reads                    = f.mapped$unmapped,
		Secondary_alignments              = f.mapped$secondary,
		Supplementary_alignments          = f.mapped$supplementary,
		Percent_mapped                    = 100 * f.mapped$mapped / trimmed,
		Alignments_after_deduplication    = f.dedup$total,
		Duplicates_removed                = f.mapped$total - f.dedup$total,
		Percent_duplicates_removed        = 100 * (f.mapped$total - f.dedup$total) / f.mapped$total,
		Alignments_MapQ20                 = f.highmq$total,
		Percent_MapQ20_of_deduplicated    = 100 * f.highmq$total / f.dedup$total,
		Median_MapQ                       = mapq.expanded.median,
		Uniquely_placed_alignments        = unique.best,
		Percent_uniquely_placed           = 100 * unique.best / unique.total,
		Alignments_with_equal_alternative  = unique.tied,
		Percent_with_equal_alternative    = 100 * unique.tied / unique.total,
		stringsAsFactors = FALSE
	)
}

df <- bind_rows(rows)
df <- left_join(df.index, df, by = "Sequencing_sample")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out_file, row.names = FALSE)

message("wrote ", out_file, " (", nrow(df), " rows)")
