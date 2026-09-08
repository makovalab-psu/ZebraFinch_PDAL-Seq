# Correlation-matrix heatmaps for Phase 1.
#
# One tile per sample pair, faceted by window size (columns) and MapQ
# branch (rows). Spearman is plotted because the coverage distributions
# are heavy tailed; Pearson is in the CSV alongside it.
#
# Form: a matrix of one magnitude, so a sequential ramp, monotonic in
# lightness, with every cell directly labelled (36 cells per facet - the
# labels are the point of a correlation matrix, not decoration). Viridis
# is used because it is perceptually uniform, CVD-safe, and is what the
# rest of the lab's figures use.
#
# Usage: Rscript Plot_correlation.R <Window_coverage_correlation.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_correlation.R <correlation.csv> <out.svg>")
}

in_file <- args[1]
out_file <- args[2]

df <- read.csv(in_file, stringsAsFactors = FALSE)

# Order samples WGS first, then by permanganate concentration, so the
# replicate blocks sit next to each other on the diagonal.
label_levels <- unique(c("WGS", sort(setdiff(df$Label_1, "WGS"))))

df$Label_1 <- factor(df$Label_1, levels = label_levels)
df$Label_2 <- factor(df$Label_2, levels = rev(label_levels))

# Facet labels: window size in human units, MapQ branch spelled out.
window_label <- function(x) {
	n <- as.numeric(x)
	ifelse(n >= 1e6, paste0(n / 1e6, " Mbp windows"), paste0(n / 1e3, " kbp windows"))
}
df$Window_label <- factor(window_label(df$Window),
                          levels = window_label(sort(unique(as.numeric(df$Window)))))
df$MapQ_label <- factor(
	ifelse(df$MapQ_branch == "All_reads", "All reads", "MapQ >= 20"),
	levels = c("All reads", "MapQ >= 20")
)

# Cell labels flip to white over the dark end of the ramp. Without this
# the low-correlation cells (the PDAL-Seq vs WGS block, which is the one
# worth reading) are dark grey text on dark purple.
fill_limits <- c(min(df$Spearman), 1)
df$Label_dark <- (df$Spearman - fill_limits[1]) / diff(fill_limits) < 0.55

p <- ggplot(df, aes(x = Label_1, y = Label_2, fill = Spearman)) +
	geom_tile(colour = "white", linewidth = 0.6) +
	geom_text(aes(label = sprintf("%.2f", Spearman), colour = Label_dark), size = 2.6) +
	facet_grid(MapQ_label ~ Window_label) +
	scale_fill_viridis(name = expression(Spearman~rho), option = "viridis",
	                   limits = fill_limits) +
	scale_colour_manual(values = c("TRUE" = "grey95", "FALSE" = "grey10"), guide = "none") +
	labs(
		title = "Pairwise correlation of window coverage, CFS414 Zebra finch",
		subtitle = "PDAL-Seq replicates and the PCR-free WGS control, bTaeGut7v0.4",
		x = NULL, y = NULL
	) +
	coord_equal() +
	theme_minimal(base_size = 9) +
	theme(
		panel.grid = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1),
		strip.text = element_text(face = "bold"),
		plot.title = element_text(face = "bold"),
		# The viewer supplies no background for a transparent plot, which
		# leaves the dark title text unreadable on a dark backdrop.
		plot.background = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 10, height = 7, units = "in", bg = "white")

message("wrote ", out_file)
