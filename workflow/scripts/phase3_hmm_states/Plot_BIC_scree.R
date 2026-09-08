# BIC scree plot: how many hidden states does the PDAL-Seq signal support?
#
# Adapted from workflow/scripts/07_PhyloHGMP_model/plot_BIC.R in the PDAL-Seq
# manuscript pipeline, which plots BIC alone. Two panels are added, both free:
#
#   * log likelihood. Every model is fit to the same windows, so a model with
#     more states cannot legitimately explain them worse. Where this curve
#     DIPS, that k's fit fell into a bad optimum, and its BIC is an artefact
#     of the fit rather than a statement about the data. Without this panel a
#     jagged BIC curve is indistinguishable from a real one.
#
#   * occupied states -- how many states the variational HMM actually assigns
#     windows to. Where that curve leaves the diagonal, extra states are being
#     fit but not used, which is independent evidence about the optimum.
#
# The dips are not a bug to be tuned away: hmmlearn initialises the emissions
# from k-means, and k-means splits wide clusters and merges narrow ones, so
# some k land badly however many times they are restarted. The remedy is to
# read the three panels together and cross-check the chosen k against
# plots/phase3/gaussian_fit/.
#
# Usage:
#   Rscript Plot_BIC_scree.R <HMM_BIC.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_BIC_scree.R <HMM_BIC.csv> <out.svg>")
}

df <- readr::read_csv(args[1], show_col_types = FALSE)

if (nrow(df) == 0) {
	stop("no models in ", args[1])
}

best <- df$States[which.min(df$BIC)]
message("lowest BIC at ", best, " states")

PANELS <- c("BIC",
            "log likelihood (must not dip)",
            "States the model actually uses")

df.long <- dplyr::bind_rows(
	data.frame(States = df$States, value = df$BIC, panel = PANELS[1]),
	data.frame(States = df$States, value = df$log_l, panel = PANELS[2]),
	data.frame(States = df$States, value = df$occupied_states, panel = PANELS[3])
)
df.long$panel <- factor(df.long$panel, levels = PANELS)

df.best <- data.frame(States = best,
                      value = min(df$BIC),
                      panel = factor(PANELS[1], levels = PANELS))

# A model with more states cannot legitimately fit the same windows worse, so
# any k below the running maximum of the log likelihood is a failed fit.
running_max <- cummax(df$log_l)
suspect <- df$log_l < running_max - .Machine$double.eps^0.5
df.suspect <- data.frame(States = df$States[suspect],
                         value  = df$log_l[suspect],
                         panel  = factor(PANELS[2], levels = PANELS))
if (nrow(df.suspect) > 0) {
	message("fits that landed in a worse optimum than a smaller k: ",
	        paste(df.suspect$States, collapse = ", "),
	        " -- treat their BIC with suspicion")
}

# Push the label into the empty half of the panel: the minimum sits at the
# bottom of a trough, so whichever side of the x range it falls on, the space
# on the other side is clear. Fixed nudges clip once k = 30 wins.
label_hjust <- if (best <= mean(range(df$States))) -0.08 else 1.08

# y = x reference for the occupancy panel only.
df.diagonal <- data.frame(
	States = df$States,
	value  = df$States,
	panel  = factor(PANELS[3], levels = PANELS)
)

P <- ggplot(df.long, aes(x = States, y = value)) +
	geom_line(data = df.diagonal, colour = "grey70", linetype = "22") +
	geom_line(colour = "#2166AC") +
	geom_point(colour = "#2166AC", size = 1.8) +
	geom_point(data = df.best, colour = "#B2182B", size = 3.4, shape = 21,
	           fill = NA, stroke = 1.1) +
	geom_text(data = df.best,
	          aes(label = paste0("lowest BIC: ", States, " states")),
	          colour = "#B2182B", vjust = -0.9, hjust = label_hjust,
	          size = 3.1) +
	geom_point(data = df.suspect, colour = "#B2182B", size = 3.4, shape = 4,
	           stroke = 1.1) +
	facet_wrap(~ panel, ncol = 1, scales = "free_y") +
	scale_x_continuous(breaks = df$States) +
	scale_y_continuous(expand = expansion(mult = 0.12)) +
	theme_classic() +
	theme(strip.background = element_blank(),
	      strip.text = element_text(hjust = 0, face = "bold"),
	      panel.grid.major.y = element_line(colour = "grey93")) +
	xlab("number of hidden states") +
	ylab(NULL)

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[2], P, width = 6, height = 8, units = "in")

message("wrote ", args[2])
