# The Figure 3 composite, one per fitted model.
#
# Layout and conventions follow
# workflow/scripts/Revised_figures/Figure_3_Ape_PDAL-Seq_MVGHMM.R in the
# PDAL-Seq manuscript: a row of panels sharing one state axis, log2 enrichment
# on a blue-white-red scale clipped at +/- log2(3), each tile labelled with the
# raw enrichment above its significance mark.
#
# Two things are deliberately NOT copied from that script:
#
#   * State.levels and State.labels are hardcoded there ("1-origin",
#     "2-recomb", ...) because the ape states had already been interpreted.
#     Ours have not, and there are seven models to compare, so states are
#     ordered by their mean PDAL-Seq signal and labelled s<i>. The ordering is
#     recomputed per model, so s3 in the k = 9 panel is not s3 in the k = 12
#     panel -- read the state id, not the row position.
#
#   * The manuscript's leftmost panel is a tile of 21 experiments. We have one,
#     so it is a bar of the mean z-scored signal with the empirical SD on it.
#
# Significance marks come from P_value, which Compile_enrichment.R has already
# set to the Wilcoxon or the empirical p according to the annotation count.
#
# Usage:
#   Rscript Plot_state_function.R <summary.csv> <nonB.csv> <functional.csv> \
#           <RNA.csv> <methylation.csv> <n_states> <out.svg>

suppressPackageStartupMessages({
	library(tidyverse)
	library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
	stop("Usage: Plot_state_function.R <summary.csv> <nonB.csv> <functional.csv> <RNA.csv> <methylation.csv> <n_states> <out.svg>")
}

summary_csv      <- args[1]
nonb_csv         <- args[2]
functional_csv   <- args[3]
rna_csv          <- args[4]
methylation_csv  <- args[5]
n_states         <- as.integer(args[6])
out_svg          <- args[7]

CLIP <- log2(3)

df.summary <- readr::read_csv(summary_csv, col_types = readr::cols(State = readr::col_character()))
df.nonB    <- readr::read_csv(nonb_csv, col_types = readr::cols(State = readr::col_character()))
df.ann     <- readr::read_csv(functional_csv, col_types = readr::cols(State = readr::col_character()))
df.RNA     <- readr::read_csv(rna_csv, col_types = readr::cols(State = readr::col_character()))
df.5mC     <- readr::read_csv(methylation_csv, col_types = readr::cols(State = readr::col_character()))

# Lowest mean signal first, so the highest-signal state sits at the top of the
# y axis (ggplot draws the first level at the bottom). An unoccupied state has
# NA for its mean; dplyr::arrange sorts NA last, so those land above the
# highest-signal state rather than being dropped.
df.summary <- df.summary %>% dplyr::arrange(empirical_mean)
State.levels <- df.summary$State
State.labels <- paste0("s", State.levels)

as_state <- function(x) factor(x, levels = State.levels, labels = State.labels)

df.summary$State <- as_state(df.summary$State)
df.nonB$State    <- as_state(df.nonB$State)
df.ann$State     <- as_state(df.ann$State)
df.RNA$State     <- as_state(df.RNA$State)
df.5mC$State     <- as_state(df.5mC$State)

BASE <- theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
	      axis.title.x = element_text(size = 8, color = "black"),
	      axis.title.y = element_blank(),
	      axis.text.y = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 7, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      legend.key.size = unit(0.1, "in"),
	      plot.title = element_text(size = 8, color = "black"),
	      plot.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "pt"))

significance <- function(p) {
	mark <- rep("ns", length(p))
	mark[which(p <= 0.05)] <- "*"
	mark[which(p <= 0.01)] <- "**"
	mark[is.na(p)] <- ""
	mark
}

heatmap_panel <- function(df, category, legend_title) {
	# Clip for colour only; the printed number stays the real enrichment, which
	# is how the manuscript figure reads a 12-fold enrichment off a scale that
	# stops at 3.
	df$fill <- pmin(pmax(df$log2Enrichment, -CLIP), CLIP)
	df$fill[!is.finite(df$log2Enrichment)] <- NA
	df$LABEL <- paste(
		ifelse(is.finite(df$Enrichment), format(round(df$Enrichment, 1), nsmall = 1), "-"),
		significance(df$P_value),
		sep = "\n"
	)
	ggplot(df, aes(x = .data[[category]], y = State, fill = fill, label = LABEL)) +
		geom_tile() +
		geom_text(color = "black", size = 1.7, lineheight = 0.8) +
		scale_fill_gradient2(low = "blue", mid = "white", high = "red",
		                     midpoint = 0, limits = c(-CLIP, CLIP),
		                     na.value = "grey85", name = legend_title) +
		BASE +
		theme(axis.title.x = element_blank())
}

####  Mean PDAL-Seq signal per state  ####

P.signal <- ggplot(df.summary, aes(y = State, x = empirical_mean)) +
	geom_bar(stat = "identity", fill = "grey40") +
	geom_errorbar(aes(xmin = empirical_mean - empirical_sd,
	                  xmax = empirical_mean + empirical_sd),
	              width = 0.3, linewidth = 0.2) +
	geom_vline(xintercept = 0, linewidth = 0.2) +
	BASE +
	# The only panel with state labels on it; the rest share this axis. No
	# axis title -- it collides with the labels at this width, and the figure
	# title already says these rows are HMM states.
	theme(axis.text.y = element_text(size = 6, color = "black")) +
	ggtitle("PDAL-Seq") +
	xlab("Mean signal\n(z-score)")

####  State size  ####

P.count <- ggplot(df.summary, aes(y = State, x = BP / 1e6,
                                  label = round(BP / 1e6, 1))) +
	geom_bar(stat = "identity", fill = "grey70") +
	geom_text(size = 1.8, hjust = -0.1) +
	scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
	BASE +
	ggtitle("State size") +
	xlab("Mbp")

####  non-B DNA motifs and functional annotations  ####

# No ggtitle on these two: the legend title names the data, and a panel with a
# legend puts its title higher than a panel without one, which left the row of
# titles stepped rather than aligned.
P.nonB <- heatmap_panel(df.nonB, "nonB", "non-B DNA motif\nlog2 enrichment")

P.ann <- heatmap_panel(df.ann, "ann", "Functional annotation\nlog2 enrichment")

####  RNA-Seq  ####

# Two libraries from the same cell line (P24 and P50), so the bar is their mean
# and the error bar is the spread between them, not a within-library estimate.
P.RNA <- ggplot(df.RNA %>% dplyr::filter(is.finite(Enrichment)),
                aes(y = State, x = Enrichment)) +
	geom_bar(stat = "summary", fun = mean, fill = "grey40") +
	stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, linewidth = 0.2) +
	geom_vline(xintercept = 1, linewidth = 0.2) +
	BASE +
	ggtitle("Transcription") +
	xlab("RNA-Seq read\nenrichment")

####  Hypomethylated CpGs  ####

P.hypo <- df.5mC %>%
	dplyr::filter(Range == "0_20") %>%
	ggplot(aes(y = State, x = Percent)) +
	geom_bar(stat = "identity", fill = "grey70") +
	BASE +
	ggtitle("Hypomethylated\nCpG sites") +
	xlab("Percent (%)")

####  Compile  ####

P <- cowplot::plot_grid(
	P.signal, P.count, P.nonB, P.ann, P.RNA, P.hypo,
	nrow = 1, align = "h", axis = "tb",
	rel_widths = c(3.2, 2.0, 5.0, 5.0, 2.4, 2.4)
)

P <- cowplot::plot_grid(
	cowplot::ggdraw() + cowplot::draw_label(
		paste0("Zebra finch CFS414 PDAL-Seq HMM, k = ", n_states, " states"),
		size = 9, hjust = 0.5
	),
	P,
	ncol = 1, rel_heights = c(0.06, 1)
)

dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)

# Height tracks the state count so the tile labels stay legible from k = 8 to
# k = 14 without hand tuning each panel.
ggsave(out_svg, P, width = 13, height = 2.2 + 0.30 * n_states,
       units = "in", limitsize = FALSE)

message("wrote ", out_svg)
