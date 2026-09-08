# Read density against non-B DNA motif density, window by window, with
# every window drawn.
#
# Phase 2's scatter removes blacklisted windows before plotting, which at
# 1 Mbp removes most of the genome. This one keeps them and ENCODES the
# thing that was being removed:
#
#   colour  blacklisted fraction of the window (viridis)
#   shape   whether the contig is one Phase 2 excluded outright
#           (chrW_mat, chrZ_pat, chrMT, the three rDNA morphs)
#   lines   solid = lm over all windows in the panel
#           dashed = lm over windows at most <light> blacklisted
#
# If the dashed line sits on the solid one, the repeats are not driving
# the correlation. If it swings away, they are, and the heatmap drawn
# from the filtered branch is answering a different question from the one
# drawn from all windows.
#
# Rows are the two series that matter -- the combined PDAL-Seq dataset and
# the WGS control -- because the question is whether PDAL-Seq tracks motif
# density in a way the control does not. The concentration-level
# experiments are in the CSV.
#
# facet_grid(..., scales = "free") gives free x per column (motif classes
# are on wildly different scales) and free y per row.
#
# Usage: Rscript Plot_nonb_scatter_all_regions.R <compare.csv> <wgs> <combined> <light> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
	stop("Usage: Plot_nonb_scatter_all_regions.R <compare.csv> <wgs_experiment> <combined_experiment> <light> <out.svg>")
}

in_file  <- args[1]
wgs      <- args[2]
combined <- args[3]
light    <- as.numeric(args[4])
out_file <- args[5]

if (!is.finite(light)) stop("light must be a numeric fraction, got '", args[4], "'")

df <- read.csv(in_file, stringsAsFactors = FALSE)

df <- df %>% filter(Experiment %in% c(wgs, combined))

if (nrow(df) == 0) stop("no rows left to plot after filtering")

nonb_levels <- c(setdiff(sort(unique(df$nonB)), "all"), "all")
df$nonB <- factor(df$nonB, levels = nonb_levels)

# Window size comes from the table rather than being hardcoded, so the
# title cannot drift from what the compile step actually wrote.
w <- as.numeric(unique(df$Window))
if (length(w) != 1) stop("compare table mixes window sizes: ", paste(w, collapse = ", "))
window_label <- if (w >= 1e6) paste0(w / 1e6, " Mbp windows") else paste0(w / 1e3, " kbp windows")

df$Series <- factor(
	ifelse(df$Experiment == wgs, "PCR-free WGS control", "PDAL-Seq (20 + 40 mM)"),
	levels = c("PDAL-Seq (20 + 40 mM)", "PCR-free WGS control")
)

df$Contig <- factor(
	ifelse(df$Phase2_excluded_contig, "Excluded by Phase 2", "Autosome"),
	levels = c("Autosome", "Excluded by Phase 2")
)

# Draw the handful of Phase-2-excluded windows last so they land on top of
# the autosomal cloud rather than under it.
df <- df[order(df$Contig), ]

# The filtered comparison line, restricted to panels that have enough
# windows to fit one. Without this, a panel with fewer than three points
# makes geom_smooth warn and drop the group anyway, but silently.
df.light <- df %>%
	filter(blacklist_fraction <= light) %>%
	group_by(Series, nonB) %>%
	filter(dplyr::n() >= 3, sd(motif_base_density) > 0) %>%
	ungroup()

n_panels <- length(unique(df$Series)) * length(unique(df$nonB))
n_light_panels <- nrow(distinct(df.light, Series, nonB))
message("comparison line drawn in ", n_light_panels, " of ", n_panels,
        " panels (windows at most ", light, " blacklisted)")

line_colour <- "#B2182B"

# Name the actual threshold on the figure rather than calling it "outside
# the blacklist", which would be Blacklist_free and is not what is drawn.
lab_all   <- "All windows"
lab_light <- sprintf("At most %g%% blacklisted", 100 * light)

p <- ggplot(df, aes(x = motif_base_density, y = CPM)) +
	geom_point(aes(colour = blacklist_fraction, shape = Contig, size = Contig),
	           alpha = 0.45, stroke = 0) +
	geom_smooth(aes(linetype = lab_all),
	            method = "lm", formula = y ~ x, se = FALSE,
	            colour = line_colour, linewidth = 0.8) +
	{ if (nrow(df.light) > 0)
		geom_smooth(data = df.light, aes(linetype = lab_light),
		            method = "lm", formula = y ~ x, se = FALSE,
		            colour = line_colour, linewidth = 0.8)
	  else NULL } +
	facet_grid(Series ~ nonB, scales = "free") +
	scale_colour_viridis_c(name = "Fraction of window blacklisted", limits = c(0, 1)) +
	scale_shape_manual(values = c("Autosome" = 16, "Excluded by Phase 2" = 17), name = NULL) +
	scale_size_manual(values = c("Autosome" = 0.45, "Excluded by Phase 2" = 1.5), name = NULL) +
	scale_linetype_manual(values = setNames(c("solid", "22"), c(lab_all, lab_light)),
	                      breaks = c(lab_all, lab_light), name = NULL) +
	guides(
		# barwidth / barheight / title.vjust were soft deprecated in ggplot2
		# 3.5; Rplot.yml does not pin ggplot2, so stay on the stable path and
		# size the key from the theme instead.
		# title.position = "top" rather than the default flush-left, which
		# collides with the bar and with its first tick label.
		colour = guide_colourbar(order = 1, title.position = "top"),
		shape  = guide_legend(order = 2, override.aes = list(size = 2, alpha = 1,
		                                                     colour = "grey30")),
		size   = "none",
		linetype = guide_legend(order = 3,
		                        override.aes = list(colour = line_colour, linewidth = 0.8))
	) +
	labs(
		title = paste("Read density versus non-B DNA motif density in", window_label),
		subtitle = paste0(
			"CFS414 Zebra finch, bTaeGut7v0.4; EVERY window drawn -- blacklisted, ",
			"sex-chromosome, chrMT and rDNA windows all kept"
		),
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
		strip.text.y = element_text(size = 7),
		plot.title = element_text(face = "bold"),
		plot.subtitle = element_text(size = 7.5),
		legend.position = "top",
		legend.box = "horizontal",
		legend.key.height = grid::unit(0.45, "lines"),
		legend.key.width  = grid::unit(1.6, "lines"),
		legend.title = element_text(size = 8, hjust = 0),
		plot.background  = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 13, height = 6, units = "in", bg = "white")

message("wrote ", out_file)
