# k-means scree plot, the model-free companion to the BIC curve.
#
# The plotting half of workflow/scripts/09_Human_HMM/Human_Kmeans_scree.R in
# the PDAL-Seq manuscript pipeline, split into its own rule so that changing
# the figure does not re-run twelve k-means fits. The variance-explained
# panel is added: total within sum of squares falls without bound as k grows,
# which makes an elbow easy to talk yourself into, and the fraction of
# variance explained puts a ceiling on the same curve.
#
# Usage:
#   Rscript Plot_Kmeans_scree.R <Kmeans_scree.csv> <out.svg>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Usage: Plot_Kmeans_scree.R <Kmeans_scree.csv> <out.svg>")
}

df <- readr::read_csv(args[1], show_col_types = FALSE)

if (nrow(df) == 0) {
	stop("no models in ", args[1])
}

if (any(!df$converged)) {
	message("NOTE: k-means did not converge at k = ",
	        paste(df$n_clusters[!df$converged], collapse = ", "),
	        "; those wss values are upper bounds")
}

df.long <- dplyr::bind_rows(
	data.frame(n_clusters = df$n_clusters, value = df$wss,
	           panel = "Total within-cluster sum of squares"),
	data.frame(n_clusters = df$n_clusters, value = df$variance_explained,
	           panel = "Fraction of variance explained")
)
df.long$panel <- factor(
	df.long$panel,
	levels = c("Total within-cluster sum of squares",
	           "Fraction of variance explained")
)

P <- ggplot(df.long, aes(x = n_clusters, y = value)) +
	geom_line(colour = "#D95F02") +
	geom_point(colour = "#D95F02", size = 1.8) +
	facet_wrap(~ panel, ncol = 1, scales = "free_y") +
	scale_x_continuous(breaks = df$n_clusters) +
	theme_classic() +
	theme(strip.background = element_blank(),
	      strip.text = element_text(hjust = 0, face = "bold"),
	      panel.grid.major.y = element_line(colour = "grey93")) +
	xlab("number of clusters (k)") +
	ylab(NULL)

dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[2], P, width = 6, height = 6, units = "in")

message("wrote ", args[2])
