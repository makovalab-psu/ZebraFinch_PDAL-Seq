# PDAL-Seq read enrichment in each repeat class, divided by the PCR-free
# control's enrichment in the same class. This is the figure the phase exists
# for.
#
# Form: a dot plot, not a heatmap. There are five PDAL-Seq libraries and the
# question is whether they AGREE -- two 20 mM replicates and two 40 mM
# replicates landing on the same value is the evidence that a class is really
# enriched, and a heatmap cell cannot show agreement. One row per class, one
# point per library, a vertical line at zero.
#
# Zero means the PDAL-Seq library is enriched in this class exactly as much as
# the PCR-free control is -- i.e. everything about the class that attracts
# reads is mapping or sequencing, not permanganate. Above zero is enrichment
# over and above that; below zero is depletion relative to it.
#
# The 0 mM library is drawn too and is the internal negative control: it went
# through the PDAL-Seq protocol without permanganate, so it should sit near
# zero wherever the treated libraries do not.
#
# Usage:
#   Rscript Plot_versus_control.R <PDALSeq_versus_control.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_versus_control.R <PDALSeq_versus_control.csv> <out.svg>")
}

in_csv  <- args[1]
out_svg <- args[2]

CLIP <- log2(8)

df <- readr::read_csv(in_csv, col_types = readr::cols())

if (nrow(df) == 0) stop("no rows in ", in_csv)

source_levels <- c("All TEs", "All satellites", "Centromere",
                   "Transposable element", "Satellite", "Other")
df$Source <- factor(df$Source, levels = intersect(source_levels, unique(df$Source)))

label_order <- c("0mM", "20mM-1", "20mM-2", "40mM-1", "40mM-2")
df$Label <- factor(df$Label, levels = intersect(label_order, unique(df$Label)))

# Same ordering as Plot_repeat_enrichment.R, computed the same way from the
# same file, so the two figures line up row for row within a panel.
repeat_order <- df %>%
	dplyr::mutate(v = dplyr::if_else(is.finite(log2Ratio), log2Ratio, NA_real_)) %>%
	dplyr::group_by(Repeat) %>%
	dplyr::summarise(rank_value = mean(v, na.rm = TRUE), .groups = "drop") %>%
	dplyr::arrange(rank_value)

df$Repeat <- factor(df$Repeat, levels = repeat_order$Repeat)

# A control enrichment of zero makes the ratio infinite, and a PDAL-Seq
# enrichment of zero makes it zero. Both are real results -- reads in only one
# of the two libraries -- but neither can be placed on a log axis, so they are
# pinned at the clip and drawn with an open symbol. Only 0/0, which says
# nothing, is dropped.
df <- df %>%
	dplyr::mutate(
		off_scale = !is.finite(log2Ratio) | abs(log2Ratio) > CLIP,
		plotted   = dplyr::case_when(
			is.na(Ratio) | is.nan(Ratio) ~ NA_real_,
			log2Ratio ==  Inf            ~  CLIP,
			log2Ratio == -Inf            ~ -CLIP,
			TRUE ~ pmax(pmin(log2Ratio, CLIP), -CLIP)
		)
	)

n_dropped <- sum(is.na(df$plotted))
if (n_dropped > 0) {
	message(n_dropped, " point(s) had an undefined ratio and are not drawn")
}

repeat_labels <- df %>%
	dplyr::distinct(Repeat, N_annotation) %>%
	dplyr::group_by(Repeat) %>%
	dplyr::summarise(N_annotation = max(N_annotation), .groups = "drop") %>%
	dplyr::arrange(Repeat)
label_lookup <- setNames(
	paste0(as.character(repeat_labels$Repeat),
	       " (", format(repeat_labels$N_annotation, big.mark = ",", trim = TRUE), ")"),
	as.character(repeat_labels$Repeat)
)

n_rows <- length(levels(df$Repeat))

p <- ggplot(df, aes(x = plotted, y = Repeat, colour = Label, shape = off_scale)) +
	geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
	geom_point(size = 1.5, alpha = 0.85, na.rm = TRUE) +
	facet_grid(rows = vars(Source), scales = "free_y", space = "free_y",
	           switch = "y") +
	scale_x_continuous(limits = c(-CLIP, CLIP),
	                   breaks = seq(-3, 3, 1)) +
	scale_y_discrete(labels = label_lookup) +
	scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
	                   labels = c(`FALSE` = "on scale", `TRUE` = "clipped"),
	                   name = NULL) +
	scale_colour_viridis_d(option = "C", end = 0.85, name = "library") +
	labs(
		x = expression(log[2] * " (PDAL-Seq enrichment / PCR-free enrichment)"),
		y = NULL,
		title = "PDAL-Seq read density in repeats, above the PCR-free control",
		subtitle = paste0("zero = the class attracts PDAL-Seq reads only as ",
		                  "much as it attracts PCR-free reads;\nopen symbols ",
		                  "are clipped at +/- 3; row label carries the number ",
		                  "of fragments")
	) +
	theme_bw(base_size = 8) +
	theme(
		axis.text.y = element_text(size = 5.5),
		strip.placement = "outside",
		strip.text.y.left = element_text(angle = 0, size = 7),
		panel.grid.minor = element_blank(),
		plot.title = element_text(size = 9, face = "bold")
	)

dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
ggsave(out_svg, p,
       width = 7,
       height = min(30, 2.2 + 0.115 * n_rows),
       limitsize = FALSE)

message("wrote ", out_svg, " (", n_rows, " repeat classes)")
