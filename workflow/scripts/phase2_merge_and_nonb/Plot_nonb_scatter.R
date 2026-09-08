# Read density against non-B DNA motif density, window by window.
#
# The heatmap gives one number per class; this shows the relationship those
# numbers summarise, so a correlation driven by a handful of outlier
# windows is visible as such.
#
# Only two series are drawn -- the combined PDAL-Seq dataset and the WGS
# control -- because the question is whether PDAL-Seq tracks motif density
# in a way the control does not. The concentration-level experiments are in
# the heatmap and the CSV. Two categorical hues, validated for CVD
# separation (protan dE 22.8) against a white surface.
#
# Usage: Rscript Plot_nonb_scatter.R <nonB_coverage_versus_read.csv> <wgs> <combined> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
	stop("Usage: Plot_nonb_scatter.R <compare.csv> <wgs_experiment> <combined_experiment> <out.svg>")
}

in_file  <- args[1]
wgs      <- args[2]
combined <- args[3]
out_file <- args[4]

df <- read.csv(in_file, stringsAsFactors = FALSE)

df <- df %>% filter(Experiment %in% c(wgs, combined), !Blacklisted)

if (nrow(df) == 0) stop("no rows left to plot after filtering")

nonb_levels <- c(setdiff(sort(unique(df$nonB)), "all"), "all")
df$nonB <- factor(df$nonB, levels = nonb_levels)

# Window size comes from the table rather than being hardcoded, so the
# title cannot drift from what compile_reads_vs_nonb actually wrote.
w <- as.numeric(unique(df$Window))
if (length(w) != 1) stop("compare table mixes window sizes: ", paste(w, collapse = ", "))
window_label <- if (w >= 1e6) paste0(w / 1e6, " Mbp windows") else paste0(w / 1e3, " kbp windows")

df$Series <- factor(
	ifelse(df$Experiment == wgs, "PCR-free WGS control", "PDAL-Seq (20 + 40 mM)"),
	levels = c("PDAL-Seq (20 + 40 mM)", "PCR-free WGS control")
)

series_colours <- c(
	"PDAL-Seq (20 + 40 mM)" = "#D95F02",
	"PCR-free WGS control"  = "#2166AC"
)

p <- ggplot(df, aes(x = motif_base_density, y = CPM, colour = Series, fill = Series)) +
	geom_point(size = 0.5, alpha = 0.25, stroke = 0) +
	geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.9, se = TRUE, alpha = 0.25) +
	facet_wrap(~ nonB, scales = "free_x", nrow = 2) +
	scale_colour_manual(values = series_colours, name = NULL) +
	scale_fill_manual(values = series_colours, name = NULL) +
	labs(
		title = paste("Read density versus non-B DNA motif density in", window_label),
		subtitle = "CFS414 Zebra finch, bTaeGut7v0.4; blacklisted and sex-chromosome windows excluded",
		x = "Fraction of window covered by motif",
		y = "Coverage per million"
	) +
	theme_minimal(base_size = 9) +
	theme(
		panel.grid.minor = element_blank(),
		# free_x scales end their tick labels at the panel edge, so adjacent
		# panels run their numbers together (0.125 / 0.000 -> "0.1250.000")
		# without extra gutter.
		panel.spacing.x = grid::unit(1.1, "lines"),
		panel.spacing.y = grid::unit(0.9, "lines"),
		strip.text = element_text(face = "bold"),
		plot.title = element_text(face = "bold"),
		legend.position = "top",
		plot.background  = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 11, height = 6, units = "in", bg = "white")

message("wrote ", out_file)
