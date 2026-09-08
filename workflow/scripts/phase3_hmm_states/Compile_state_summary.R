# Collect the per-state parameters of every fitted model into one table.
#
# The twelve Parameters.csv files are written by extract_GHMM_parameters.py,
# which already computed both the Gaussian parameters and the empirical mean
# and sd of the windows each state holds -- the job the manuscript splits
# between 07_PhyloHGMP_model/extract_GHMM_parameters.py and
# 07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R. This rule only stacks
# them, so it stays cheap no matter how many models there are.
#
# Usage:
#   Rscript Compile_state_summary.R <out.csv> <Parameters.csv> [<Parameters.csv> ...]

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	stop("Usage: Compile_state_summary.R <out.csv> <Parameters.csv> ...")
}

out_csv <- args[1]
inputs  <- args[-1]

df <- dplyr::bind_rows(lapply(inputs, function(path) {
	readr::read_csv(path, show_col_types = FALSE)
}))

df <- df %>%
	dplyr::mutate(gaussian_sd = sqrt(gaussian_var)) %>%
	dplyr::arrange(n_states, feature, gaussian_mean)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(df, out_csv)

empty <- df %>% dplyr::filter(n_windows == 0)
if (nrow(empty) > 0) {
	message("states left unused by the variational fit:")
	print(as.data.frame(empty %>% dplyr::count(n_states, name = "unused_states")))
}

message("wrote ", nrow(df), " rows across ",
        length(unique(df$n_states)), " models to ", out_csv)
