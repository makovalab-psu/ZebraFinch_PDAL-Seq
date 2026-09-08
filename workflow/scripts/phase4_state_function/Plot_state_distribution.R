# Where each HMM state sits in the genome, one panel per chromosome.
#
# This is panel B of Figure_3_Ape_PDAL-Seq_MVGHMM.R (its P.percent): the
# percent of each chromosome held by each state, stacked to 100%. The
# manuscript's "nonsyntenic" remainder is "unsegmented" here -- blacklisted
# windows plus the contigs Phase 3 excluded outright (chrW, chrZ, chrMT and the
# rDNA morphs), which is why those contigs appear as a solid grey bar.
#
# Usage:
#   Rscript Plot_state_distribution.R <State_by_chromosome.csv> <n_states> <out.svg>

suppressPackageStartupMessages({
	library(tidyverse)
	library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	stop("Usage: Plot_state_distribution.R <State_by_chromosome.csv> <n_states> <out.svg>")
}

df <- readr::read_csv(args[1], col_types = readr::cols(State = readr::col_character()))
n_states <- as.integer(args[2])
out_svg  <- args[3]

# The file is written in .fai order; keep it rather than letting factor() sort
# chr10 before chr2.
CHR <- unique(df$chr)
df$chr <- factor(df$chr, levels = CHR)

State.levels <- c("unsegmented", as.character(seq_len(n_states) - 1L))
State.labels <- c("unsegmented", paste0("s", seq_len(n_states) - 1L))
df$State <- factor(df$State, levels = State.levels, labels = State.labels)

State_colors <- c("grey80", viridis::viridis(n_states, option = "turbo"))

P <- ggplot(df, aes(x = chr, y = Percent, fill = State)) +
	geom_bar(stat = "identity") +
	facet_wrap(~chr, nrow = 1, scales = "free_x") +
	scale_fill_manual(values = State_colors) +
	scale_y_continuous(limits = c(0, 100.01)) +
	theme_classic() +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.title.x = element_blank(),
	      legend.position = "top",
	      legend.title = element_text(size = 8, color = "black"),
	      legend.text = element_text(size = 6, color = "black"),
	      legend.key.size = unit(0.1, "in"),
	      axis.title.y = element_text(size = 8, color = "black"),
	      axis.text.y = element_text(size = 6, color = "black"),
	      plot.title = element_text(size = 8, color = "black"),
	      strip.text = element_text(size = 5, color = "black", angle = 90),
	      strip.background = element_blank(),
	      panel.spacing = unit(0, "lines")) +
	ggtitle(paste0("Distribution of Zebra finch PDAL-Seq HMM states (k = ",
	               n_states, ") across chromosomes")) +
	ylab("Percent of chromosome (%)") +
	guides(fill = guide_legend(nrow = 1))

dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
ggsave(out_svg, P, width = 13, height = 3.5, units = "in", limitsize = FALSE)

message("wrote ", out_svg)
