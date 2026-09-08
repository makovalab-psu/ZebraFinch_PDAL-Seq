# Does each state's fitted Gaussian mean match the windows it was given?
#
# workflow/scripts/07_PhyloHGMP_model/Compare_GHMM_model_to_means.R makes this
# comparison one model at a time, as two heatmaps and a scatter. Here it is
# one figure across every model, because the useful reading is comparative:
# a model whose points leave the y = x line is not describing its own states,
# and that is a reason to reject a k regardless of what BIC says about it.
#
# Point colour is the fraction of the genome in the state, so a state that
# sits off the line but holds 0.1% of windows is visibly not the same problem
# as one that holds 20%.
#
# Usage:
#   Rscript Plot_model_vs_empirical_means.R <HMM_state_summary.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_model_vs_empirical_means.R <HMM_state_summary.csv> <out.svg>")
}

df <- readr::read_csv(args[1], show_col_types = FALSE)

# Unused states have no windows and therefore no empirical mean to compare to.
df <- df %>% dplyr::filter(n_windows > 0)

if (nrow(df) == 0) {
	stop("no occupied states in ", args[1])
}

df$panel <- factor(paste0("k = ", df$n_states),
                   levels = paste0("k = ", sort(unique(df$n_states))))

if (length(unique(df$feature)) > 1) {
	df$panel <- interaction(df$panel, df$feature, sep = "\n", lex.order = TRUE)
}

worst <- df %>%
	dplyr::mutate(gap = abs(gaussian_mean - empirical_mean)) %>%
	dplyr::slice_max(gap, n = 1)
message("largest model/data gap: k = ", worst$n_states[1],
        ", state ", worst$state[1],
        ", fitted ", signif(worst$gaussian_mean[1], 4),
        " vs observed ", signif(worst$empirical_mean[1], 4),
        " (", signif(100 * worst$weight[1], 3), "% of windows)")

n_col <- min(4, nlevels(df$panel))
n_row <- ceiling(nlevels(df$panel) / n_col)

# Equal x and y limits rather than coord_fixed(): coord_fixed inside
# facet_wrap squashes the panels to different heights and silently drops the
# axis from most of them. Square panels come from the width/height below
# instead, and the y = x line reads the same either way.
lim <- range(c(df$gaussian_mean, df$empirical_mean), na.rm = TRUE)
lim <- lim + c(-1, 1) * 0.06 * diff(lim)

P <- ggplot(df, aes(x = gaussian_mean, y = empirical_mean)) +
	geom_abline(intercept = 0, slope = 1, colour = "grey70", linetype = "22") +
	geom_point(aes(fill = 100 * weight), shape = 21, colour = "white",
	           stroke = 0.3, size = 2.6) +
	scale_fill_viridis(name = "% of windows", option = "viridis",
	                   trans = "sqrt") +
	scale_x_continuous(limits = lim) +
	scale_y_continuous(limits = lim) +
	facet_wrap(~ panel, ncol = n_col) +
	theme_classic() +
	theme(strip.background = element_blank(),
	      strip.text = element_text(hjust = 0, face = "bold"),
	      legend.position = "bottom",
	      legend.key.width = grid::unit(1.4, "lines")) +
	xlab("fitted Gaussian mean (z)") +
	ylab("mean of the windows in the state (z)")

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[2], P,
       width = 2.2 * n_col + 1.0,
       height = 2.2 * n_row + 1.4,
       units = "in", limitsize = FALSE)

message("wrote ", args[2])
