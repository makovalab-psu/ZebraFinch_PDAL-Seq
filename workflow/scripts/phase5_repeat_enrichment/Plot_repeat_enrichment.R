# Heatmap of read enrichment in every repeat class, one column per library.
#
# Form: a matrix of one statistic that carries POLARITY -- a class can be
# depleted as easily as enriched -- so the ramp is diverging with symmetric
# limits and neutral grey at zero, as in Phase 2's non-B heatmap. A sequential
# ramp would hide the sign, which is the result.
#
# The PCR-free column is the point of the figure and is drawn first, on the
# left. Any class where the WGS column is as red as the PDAL-Seq columns is a
# mapping or sequencing artefact, not permanganate reactivity; that is the
# comparison Plot_versus_control.R makes explicit.
#
# Rows are ordered by mean log2(PDAL-Seq / control) within each panel -- the
# same order Plot_versus_control.R uses -- so the two figures can be read side
# by side.
#
# Fill is clipped at +/- log2(8). Satellite arrays reach enrichments far
# outside that and would otherwise flatten every other row to grey; the
# unclipped numbers are in results/phase5/Repeat_enrichment.csv.
#
# Usage:
#   Rscript Plot_repeat_enrichment.R <Repeat_enrichment.csv> \
#           <PDALSeq_versus_control.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Plot_repeat_enrichment.R <Repeat_enrichment.csv> <PDALSeq_versus_control.csv> <out.svg>")
}

enrichment_csv <- args[1]
versus_csv     <- args[2]
out_svg        <- args[3]

CLIP <- log2(8)

df      <- readr::read_csv(enrichment_csv, col_types = readr::cols())
df.vs   <- readr::read_csv(versus_csv, col_types = readr::cols())

if (nrow(df) == 0) stop("no rows in ", enrichment_csv)

# Panels top to bottom: the two pooled classes and the centromere first, then
# the two long family lists.
source_levels <- c("All TEs", "All satellites", "Centromere",
                   "Transposable element", "Satellite", "Other")
df$Source <- factor(df$Source, levels = intersect(source_levels, unique(df$Source)))

# Libraries left to right: control first, then increasing permanganate.
label_order <- c("WGS", "0mM", "20mM-1", "20mM-2", "40mM-1", "40mM-2")
df$Label <- factor(df$Label, levels = intersect(label_order, unique(df$Label)))
if (any(is.na(df$Label))) {
	stop("a sample label is not in the expected set: ",
	     paste(unique(as.character(df$Label[is.na(df$Label)])), collapse = ", "))
}

# One order for both Phase 5 figures. Classes whose ratio is undefined (a
# control enrichment of zero) sort last rather than being dropped.
repeat_order <- df.vs %>%
	dplyr::mutate(log2Ratio = dplyr::if_else(is.finite(log2Ratio), log2Ratio, NA_real_)) %>%
	dplyr::group_by(Repeat) %>%
	dplyr::summarise(rank_value = mean(log2Ratio, na.rm = TRUE), .groups = "drop") %>%
	dplyr::arrange(rank_value)

df$Repeat <- factor(df$Repeat,
                    levels = c(repeat_order$Repeat,
                               setdiff(unique(df$Repeat), repeat_order$Repeat)))

df$fill <- pmax(pmin(df$log2Enrichment, CLIP), -CLIP)

stars <- function(p) {
	dplyr::case_when(
		is.na(p)   ~ "",
		p < 0.001  ~ "***",
		p < 0.01   ~ "**",
		p < 0.05   ~ "*",
		TRUE       ~ ""
	)
}
df$mark <- stars(df$P_value)

# A class with fewer annotations than the subsample size is the SAME sample
# every iteration, so its enrichment is one number rather than a mean over 100
# draws. The count is on the axis so a spectacular enrichment resting on three
# annotations cannot be read as if it rested on three thousand.
repeat_labels <- df %>%
	dplyr::group_by(Repeat) %>%
	dplyr::summarise(N_annotation = max(N_annotation), .groups = "drop") %>%
	dplyr::arrange(Repeat)
label_lookup <- setNames(
	paste0(as.character(repeat_labels$Repeat),
	       " (", format(repeat_labels$N_annotation, big.mark = ",", trim = TRUE), ")"),
	as.character(repeat_labels$Repeat)
)

n_rows <- length(unique(df$Repeat))

p <- ggplot(df, aes(x = Label, y = Repeat, fill = fill)) +
	geom_tile(colour = "white", linewidth = 0.2) +
	geom_text(aes(label = mark), size = 2.4, vjust = 0.75, colour = "black") +
	facet_grid(rows = vars(Source), scales = "free_y", space = "free_y",
	           switch = "y") +
	scale_fill_gradient2(
		low = "#2166AC", mid = "grey92", high = "#B2182B",
		midpoint = 0, limits = c(-CLIP, CLIP),
		breaks = seq(-3, 3, 1),
		labels = c("<= -3", "-2", "-1", "0", "1", "2", ">= 3"),
		name = expression(log[2] * " enrichment")
	) +
	scale_y_discrete(labels = label_lookup) +
	labs(
		x = NULL, y = NULL,
		title = "Read enrichment in repeat classes, whole genome, all reads",
		subtitle = paste0("summed per-base coverage over annotation fragments ",
		                  "vs the same fragments shuffled within chromosome;\n",
		                  "* p < 0.05, ** p < 0.01, *** p < 0.001; ",
		                  "row label carries the number of fragments")
	) +
	theme_bw(base_size = 8) +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1),
		axis.text.y = element_text(size = 5.5),
		strip.placement = "outside",
		strip.text.y.left = element_text(angle = 0, size = 7),
		panel.grid = element_blank(),
		legend.position = "right",
		plot.title = element_text(size = 9, face = "bold")
	)

dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
ggsave(out_svg, p,
       width = 6.5,
       height = min(30, 2.2 + 0.115 * n_rows),
       limitsize = FALSE)

message("wrote ", out_svg, " (", n_rows, " repeat classes x ",
        length(unique(df$Label)), " libraries)")
