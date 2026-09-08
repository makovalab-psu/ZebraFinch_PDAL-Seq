# Heatmap of read density versus non-B DNA motif density, on EVERY window
# of the genome.
#
# The same figure as plots/phase2/nonB_correlation_heatmap.svg, drawn from
# Region_set == "All_windows" and Chromosome_set == "All_contigs" instead
# of Phase 2's blacklist-filtered, sex-chromosome-excluded subset. Read
# the two side by side: where a cell moves, repeats were carrying the
# signal.
#
# Form, ramp and label-contrast logic are deliberately identical to the
# Phase 2 heatmap so the two are visually comparable. The diverging ramp
# carries POLARITY -- a class can be depleted as easily as enriched -- and
# is symmetric about zero. Text colour is chosen from the ACTUAL WCAG
# luminance of each cell's fill, because a diverging ramp is lightest in
# the MIDDLE and "55% along the ramp" would put white text on pale cells.
#
# Cells the compile step could not compute (fewer than three windows, or
# no variance) are drawn grey and labelled n/a rather than dropped, so a
# missing combination is visible as missing.
#
# Usage: Rscript Plot_nonb_correlation_all_regions.R <nonB_correlation_coefficient.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_nonb_correlation_all_regions.R <correlation.csv> <out.svg>")
}

in_file  <- args[1]
out_file <- args[2]

df <- read.csv(in_file, stringsAsFactors = FALSE)

df <- df %>% filter(
	Signal == "CPM",
	Region_set == "All_windows",
	Chromosome_set == "All_contigs"
)

if (nrow(df) == 0) stop("no rows left to plot after filtering")

n_na <- sum(is.na(df$Spearman))
if (n_na > 0) message("note: ", n_na, " of ", nrow(df), " cells have no correlation")
if (all(is.na(df$Spearman))) stop("every cell is NA -- nothing to plot")

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
lim <- max(abs(df$Spearman), na.rm = TRUE)
lim <- ceiling(lim * 20) / 20
if (!is.finite(lim) || lim == 0) lim <- 0.05

# Rebuild each cell's fill and pick whichever of white/black has the better
# WCAG contrast against it.
fill_rgb <- colorRamp(c("#2166AC", "#F2F2F0", "#B2182B"), space = "Lab")(
	pmin(pmax((df$Spearman + lim) / (2 * lim), 0), 1)
)

relative_luminance <- function(m) {
	cs <- m / 255
	cs <- ifelse(cs <= 0.03928, cs / 12.92, ((cs + 0.055) / 1.055) ^ 2.4)
	0.2126 * cs[, 1] + 0.7152 * cs[, 2] + 0.0722 * cs[, 3]
}

L <- relative_luminance(fill_rgb)
df$Label_light <- (1.05 / (L + 0.05)) > ((L + 0.05) / 0.05)
# NA cells get the grey na.value fill, which is dark text either way.
df$Label_light[is.na(df$Spearman)] <- FALSE
df$Cell_label <- ifelse(is.na(df$Spearman), "n/a", sprintf("%.2f", df$Spearman))

p <- ggplot(df, aes(x = nonB, y = Label, fill = Spearman)) +
	geom_tile(colour = "white", linewidth = 0.6) +
	geom_text(aes(label = Cell_label, colour = Label_light), size = 2.5) +
	facet_grid(Metric_label ~ Window_label) +
	scale_fill_gradient2(
		name = expression(Spearman~rho),
		low = "#2166AC", mid = "#F2F2F0", high = "#B2182B",
		midpoint = 0, limits = c(-lim, lim), na.value = "#D9D9D9"
	) +
	scale_colour_manual(values = c("TRUE" = "grey97", "FALSE" = "grey10"), guide = "none") +
	labs(
		title = "Read density versus non-B DNA motif density, CFS414 Zebra finch",
		subtitle = "Spearman rho per genome window, EVERY window included -- blacklisted, sex-chromosome, chrMT and rDNA windows all kept",
		x = "Non-B DNA motif class", y = NULL
	) +
	theme_minimal(base_size = 9) +
	theme(
		panel.grid = element_blank(),
		strip.text = element_text(face = "bold"),
		plot.title = element_text(face = "bold"),
		plot.subtitle = element_text(size = 7.5),
		# A transparent background leaves the dark title unreadable in a
		# dark viewer; paint it explicitly, as in Phases 1 and 2.
		plot.background  = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 11, height = 6, units = "in", bg = "white")

message("wrote ", out_file)
