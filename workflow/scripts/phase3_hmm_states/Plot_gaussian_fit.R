# Do the fitted Gaussians actually describe the signal?
#
# Adapted from workflow/scripts/07_PhyloHGMP_model/plot_to_test_Guassian.R in
# the PDAL-Seq manuscript pipeline. Three differences:
#
#   * it reads a pre-binned histogram written by extract_GHMM_parameters.py
#     instead of ~10^6 raw windows, so twelve of these can run concurrently;
#   * it is built from ggplot2 facets rather than nested gridExtra::grid.arrange
#     calls, because gridExtra is not in Rplot.yml and adding it there would
#     rehash the environment and re-run every finished Phase 1/2 figure;
#   * it adds the panel that actually answers the question. The manuscript
#     plots each state's windows against that state's Gaussian, which can look
#     fine state by state while the mixture still misses the overall
#     distribution. The first panel here is every window against
#     sum_k w_k * dnorm(x; mu_k, sd_k).
#
# Densities are on the same footing as the curves: per-state bars integrate to
# 1 over that state, the "all windows" bars integrate to 1 over the genome.
#
# Usage:
#   Rscript Plot_gaussian_fit.R <Parameters.csv> <Histogram.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Plot_gaussian_fit.R <Parameters.csv> <Histogram.csv> <out.svg>")
}

df.par  <- readr::read_csv(args[1], show_col_types = FALSE)
df.hist <- readr::read_csv(args[2], col_types = readr::cols(state = readr::col_character(),
                                                            .default = readr::col_guess()))

Features <- unique(df.par$feature)
N.states <- unique(df.par$n_states)
if (length(N.states) != 1) {
	stop("Parameters.csv mixes models with ", paste(N.states, collapse = ", "), " states")
}

GLOBAL <- "all windows"
CURVE_POINTS <- 400

# Panels are ordered by fitted mean so the figure reads left to right from the
# lowest-signal state to the highest, with the overall fit first.
df.order <- df.par %>%
	dplyr::filter(feature == Features[1]) %>%
	dplyr::arrange(gaussian_mean)

# One label per (state, feature). Built as a lookup rather than computed in
# two places, so the histogram, the curves and the factor levels can never
# disagree about what a panel is called.
df.label <- dplyr::bind_rows(
	df.par %>%
		dplyr::transmute(
			state = as.character(state),
			feature = feature,
			label = sprintf("state %s  (%.1f%% of windows)", state, 100 * weight)
		),
	data.frame(state = "global", feature = Features, label = GLOBAL,
	           stringsAsFactors = FALSE)
)

if (length(Features) > 1) {
	df.label$label <- paste0(df.label$label, "\n", df.label$feature)
}

panel.label <- function(state, feature) {
	df.label$label[match(paste(state, feature),
	                     paste(df.label$state, df.label$feature))]
}

df.hist$panel <- panel.label(df.hist$state, df.hist$feature)

if (anyNA(df.hist$panel)) {
	stop("Histogram.csv holds a state that Parameters.csv does not describe")
}

panel_levels <- unlist(lapply(Features, function(f) {
	panel.label(c("global", as.character(df.order$state)), f)
}))
panel_levels <- unique(panel_levels[panel_levels %in% df.hist$panel])
df.hist$panel <- factor(df.hist$panel, levels = panel_levels)

####Fitted curves#####################################################

curves <- list()
k <- 0

for (feature in Features) {

	par.f <- df.par %>% dplyr::filter(feature == !!feature)

	for (state in c("global", as.character(par.f$state))) {

		hist.i <- df.hist %>%
			dplyr::filter(state == !!state, feature == !!feature)

		if (nrow(hist.i) == 0) {
			# An empty state has no windows and therefore no histogram; the
			# variational HMM is allowed to leave states unused.
			next
		}

		half <- hist.i$bin_width[1] / 2
		x <- seq(min(hist.i$bin_mid) - half, max(hist.i$bin_mid) + half,
		         length.out = CURVE_POINTS)

		if (state == "global") {
			# The mixture, weighted by how many windows each state holds.
			y <- rep(0, length(x))
			for (i in seq_len(nrow(par.f))) {
				if (!is.finite(par.f$gaussian_var[i]) || par.f$gaussian_var[i] <= 0) next
				y <- y + par.f$weight[i] *
					dnorm(x, mean = par.f$gaussian_mean[i],
					      sd = sqrt(par.f$gaussian_var[i]))
			}
		} else {
			i <- which(as.character(par.f$state) == state)
			if (!is.finite(par.f$gaussian_var[i]) || par.f$gaussian_var[i] <= 0) next
			y <- dnorm(x, mean = par.f$gaussian_mean[i],
			           sd = sqrt(par.f$gaussian_var[i]))
		}

		k <- k + 1
		curves[[k]] <- data.frame(
			panel   = panel.label(state, feature),
			feature = feature,
			x       = x,
			y       = y
		)
	}
}

df.curve <- dplyr::bind_rows(curves)
df.curve$panel <- factor(df.curve$panel, levels = levels(df.hist$panel))

####Plot##############################################################

n_panels <- nlevels(df.hist$panel)
n_col    <- min(4, n_panels)
n_row    <- ceiling(n_panels / n_col)

# geom_col ignores `width` supplied as an aesthetic, so it goes in as a
# parameter. The bins come off one shared grid per feature, so there is a
# single width unless several features are being plotted at once; with more
# than one, fall back to ggplot's default (0.9 * the bin spacing).
bar_widths <- unique(df.hist$bin_width)
bar_width  <- if (length(bar_widths) == 1) bar_widths else NULL

P <- ggplot() +
	geom_col(data = df.hist,
	         mapping = aes(x = bin_mid, y = density),
	         width = bar_width, fill = "#BDBDBD") +
	geom_line(data = df.curve,
	          mapping = aes(x = x, y = y),
	          colour = "#B2182B", linewidth = 0.7) +
	facet_wrap(~ panel, ncol = n_col, scales = "free") +
	theme_classic() +
	theme(strip.background = element_blank(),
	      strip.text = element_text(hjust = 0, size = 8),
	      panel.spacing.x = grid::unit(1.1, "lines"),
	      axis.text = element_text(size = 7)) +
	xlab("normalized PDAL-Seq signal (z)") +
	ylab("density") +
	ggtitle(paste0(N.states, "-state model: fitted Gaussians against the signal"))

dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[3], P,
       width = 2.6 * n_col + 0.6,
       height = 2.1 * n_row + 0.8,
       units = "in", limitsize = FALSE)

message("wrote ", args[3], " (", n_panels, " panels)")
