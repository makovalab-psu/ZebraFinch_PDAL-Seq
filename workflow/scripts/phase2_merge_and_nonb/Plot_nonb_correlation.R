# Heatmap of read density versus non-B DNA motif density.
#
# Form: a matrix of one statistic, but unlike the Phase 1 replicate matrix
# this one carries POLARITY -- a class can be depleted as easily as
# enriched -- so the ramp is diverging (two hues, neutral grey at zero)
# with symmetric limits, not sequential. A sequential ramp here would hide
# the sign, which is the whole result.
#
# Every cell is labelled: 40 cells per panel and the numbers are the point.
#
# Shows the CPM signal on blacklist-filtered windows. The ratio signal and
# the unfiltered branch are in the CSV.
#
# Usage: Rscript Plot_nonb_correlation.R <nonB_correlation_coefficient.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_nonb_correlation.R <correlation.csv> <out.svg>")
}

in_file  <- args[1]
out_file <- args[2]

df <- read.csv(in_file, stringsAsFactors = FALSE)

df <- df %>% filter(Signal == "CPM", Blacklist_filtered == TRUE, !is.na(Spearman))

if (nrow(df) == 0) stop("no rows left to plot after filtering")

# Motif classes in release order, with the combined set last.
nonb_levels <- c(setdiff(sort(unique(df$nonB)), "all"), "all")
df$nonB <- factor(df$nonB, levels = nonb_levels)

# Experiments: control first, then increasing permanganate, then combined.
label_order <- c("WGS", "0mM", "20mM", "40mM", "combined")
df$Label <- factor(df$Label, levels = rev(intersect(label_order, unique(df$Label))))

window_label <- function(x) {
	n <- as.numeric(x)
	ifelse(n >= 1e6, paste0(n / 1e6, " Mbp windows"), paste0(n / 1e3, " kbp windows"))
}
df$Window_label <- factor(window_label(df$Window),
                          levels = window_label(sort(unique(as.numeric(df$Window)))))

df$Metric_label <- factor(
	ifelse(df$Metric == "count", "Motif count", "Motif base density"),
	levels = c("Motif count", "Motif base density")
)

# Symmetric limits so zero sits exactly on the neutral midpoint.
lim <- max(abs(df$Spearman))
lim <- ceiling(lim * 20) / 20

# Text colour is chosen from the ACTUAL luminance of each cell's fill, not
# from a fraction of the ramp. Phase 1 could use "55% along the ramp"
# because viridis is monotonic in lightness; a diverging ramp is lightest
# in the MIDDLE, so the same rule puts white text on pale mid-tone cells.
# Rebuild the fill here and pick whichever of white/black has the better
# WCAG contrast against it.
fill_rgb <- colorRamp(c("#2166AC", "#F2F2F0", "#B2182B"), space = "Lab")(
	(df$Spearman + lim) / (2 * lim)
)

relative_luminance <- function(m) {
	cs <- m / 255
	cs <- ifelse(cs <= 0.03928, cs / 12.92, ((cs + 0.055) / 1.055) ^ 2.4)
	0.2126 * cs[, 1] + 0.7152 * cs[, 2] + 0.0722 * cs[, 3]
}

L <- relative_luminance(fill_rgb)
df$Label_light <- (1.05 / (L + 0.05)) > ((L + 0.05) / 0.05)

p <- ggplot(df, aes(x = nonB, y = Label, fill = Spearman)) +
	geom_tile(colour = "white", linewidth = 0.6) +
	geom_text(aes(label = sprintf("%.2f", Spearman), colour = Label_light), size = 2.5) +
	facet_grid(Metric_label ~ Window_label) +
	scale_fill_gradient2(
		name = expression(Spearman~rho),
		low = "#2166AC", mid = "#F2F2F0", high = "#B2182B",
		midpoint = 0, limits = c(-lim, lim)
	) +
	scale_colour_manual(values = c("TRUE" = "grey97", "FALSE" = "grey10"), guide = "none") +
	labs(
		title = "Read density versus non-B DNA motif density, CFS414 Zebra finch",
		subtitle = "Spearman rho per genome window, blacklisted and sex-chromosome windows excluded",
		x = "Non-B DNA motif class", y = NULL
	) +
	theme_minimal(base_size = 9) +
	theme(
		panel.grid = element_blank(),
		strip.text = element_text(face = "bold"),
		plot.title = element_text(face = "bold"),
		# A transparent background leaves the dark title unreadable in a
		# dark viewer; paint it explicitly, as in Phase 1.
		plot.background  = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 11, height = 6, units = "in", bg = "white")

message("wrote ", out_file)
