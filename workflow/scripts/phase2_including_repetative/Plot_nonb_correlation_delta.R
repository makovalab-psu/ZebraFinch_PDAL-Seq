# What does including the repetitive regions actually change?
#
#   delta rho = rho(All_windows) - rho(Blacklist_light)
#
# per motif class x window size x motif metric, on all contigs, CPM
# signal. A cell at zero means the blacklist was irrelevant to that
# correlation. A cell away from zero means the repeats were carrying
# signal that Phase 2's filter removed, and the sign says which way.
#
# The comparison uses Blacklist_light (blacklisted fraction <= the
# configured threshold) rather than Blacklist_free (fraction == 0).
# merge_blacklist runs `bedtools merge -d 1`, so blacklist intervals are
# long -- on the manuscript's human blacklist, 14% of autosomal sequence
# in intervals of mean 65 kbp touches 93% of 1 Mbp windows. Blacklist_free
# is therefore a tiny and strongly biased slice at coarse window sizes,
# which would make a difference against it meaningless. Blacklist_free is
# still in the CSV, where it is exactly Phase 2's filtered branch.
#
# Cells where either side is missing, or where the filtered side has too
# few windows to trust, are drawn grey and labelled n/a. They are not
# dropped -- a combination that cannot be compared is itself a result.
#
# Usage: Rscript Plot_nonb_correlation_delta.R <nonB_correlation_coefficient.csv> <light> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Plot_nonb_correlation_delta.R <correlation.csv> <light> <out.svg>")
}

in_file  <- args[1]
light    <- as.numeric(args[2])
out_file <- args[3]

if (!is.finite(light)) stop("light must be a numeric fraction, got '", args[2], "'")

# Below this many windows the filtered correlation is not worth
# differencing against.
MIN_WINDOWS <- 30

df <- read.csv(in_file, stringsAsFactors = FALSE)

df <- df %>% filter(
	Signal == "CPM",
	Chromosome_set == "All_contigs",
	Region_set %in% c("All_windows", "Blacklist_light")
)

if (nrow(df) == 0) stop("no rows left to plot after filtering")

key <- c("Experiment", "Label", "nonB", "Window", "Metric")

df.all <- df %>% filter(Region_set == "All_windows") %>%
	select(all_of(key), rho_all = Spearman, n_all = N_windows)

df.flt <- df %>% filter(Region_set == "Blacklist_light") %>%
	select(all_of(key), rho_flt = Spearman, n_flt = N_windows)

df.d <- inner_join(df.all, df.flt, by = key)

if (nrow(df.d) == 0) stop("All_windows and Blacklist_light do not share a single combination")
if (nrow(df.d) != nrow(df.all)) {
	stop("join is not one to one: ", nrow(df.all), " All_windows rows became ", nrow(df.d))
}

df.d$Delta <- df.d$rho_all - df.d$rho_flt
df.d$Delta[df.d$n_flt < MIN_WINDOWS] <- NA_real_

n_na <- sum(is.na(df.d$Delta))
if (n_na > 0) {
	message("note: ", n_na, " of ", nrow(df.d),
	        " cells cannot be compared (missing correlation, or fewer than ",
	        MIN_WINDOWS, " windows outside the blacklist)")
}
if (all(is.na(df.d$Delta))) {
	stop("no cell can be compared -- Blacklist_light is empty or degenerate at every window size")
}

nonb_levels <- c(setdiff(sort(unique(df.d$nonB)), "all"), "all")
df.d$nonB <- factor(df.d$nonB, levels = nonb_levels)

label_order <- c("WGS", "0mM", "20mM", "40mM", "combined")
df.d$Label <- factor(df.d$Label, levels = rev(intersect(label_order, unique(df.d$Label))))

window_label <- function(x) {
	n <- as.numeric(x)
	ifelse(n >= 1e6, paste0(n / 1e6, " Mbp windows"), paste0(n / 1e3, " kbp windows"))
}
df.d$Window_label <- factor(window_label(df.d$Window),
                            levels = window_label(sort(unique(as.numeric(df.d$Window)))))

df.d$Metric_label <- factor(
	ifelse(df.d$Metric == "count", "Motif count", "Motif base density"),
	levels = c("Motif count", "Motif base density")
)

lim <- max(abs(df.d$Delta), na.rm = TRUE)
lim <- ceiling(lim * 100) / 100
if (!is.finite(lim) || lim == 0) lim <- 0.01

# Same WCAG luminance rule as the other heatmaps: a diverging ramp is
# lightest in the middle, so text colour must come from the fill, not from
# a fraction of the ramp.
fill_rgb <- colorRamp(c("#2166AC", "#F2F2F0", "#B2182B"), space = "Lab")(
	pmin(pmax((df.d$Delta + lim) / (2 * lim), 0), 1)
)

relative_luminance <- function(m) {
	cs <- m / 255
	cs <- ifelse(cs <= 0.03928, cs / 12.92, ((cs + 0.055) / 1.055) ^ 2.4)
	0.2126 * cs[, 1] + 0.7152 * cs[, 2] + 0.0722 * cs[, 3]
}

L <- relative_luminance(fill_rgb)
df.d$Label_light <- (1.05 / (L + 0.05)) > ((L + 0.05) / 0.05)
df.d$Label_light[is.na(df.d$Delta)] <- FALSE
df.d$Cell_label <- ifelse(is.na(df.d$Delta), "n/a", sprintf("%+.2f", df.d$Delta))

p <- ggplot(df.d, aes(x = nonB, y = Label, fill = Delta)) +
	geom_tile(colour = "white", linewidth = 0.6) +
	geom_text(aes(label = Cell_label, colour = Label_light), size = 2.5) +
	facet_grid(Metric_label ~ Window_label) +
	scale_fill_gradient2(
		name = expression(Delta~rho),
		low = "#2166AC", mid = "#F2F2F0", high = "#B2182B",
		midpoint = 0, limits = c(-lim, lim), na.value = "#D9D9D9"
	) +
	scale_colour_manual(values = c("TRUE" = "grey97", "FALSE" = "grey10"), guide = "none") +
	labs(
		title = "What including the repetitive regions does to the correlation",
		subtitle = paste0(
			"Spearman rho on all windows minus rho on windows at most ",
			sprintf("%g%%", 100 * light), " blacklisted; ",
			"positive means the repetitive regions raise the correlation"
		),
		x = "Non-B DNA motif class", y = NULL
	) +
	theme_minimal(base_size = 9) +
	theme(
		panel.grid = element_blank(),
		strip.text = element_text(face = "bold"),
		plot.title = element_text(face = "bold"),
		plot.subtitle = element_text(size = 7.5),
		plot.background  = element_rect(fill = "white", colour = NA),
		panel.background = element_rect(fill = "white", colour = NA)
	)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, p, width = 11, height = 6, units = "in", bg = "white")

message("wrote ", out_file)
