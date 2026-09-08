library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

####Global variables####

theme_set(theme_classic(base_size = 6) +
            theme(text = element_text(color = "black", size = 6),
                  axis.text = element_text(color = "black", size = 6),
                  axis.title = element_text(color = "black", size = 6),
                  strip.text = element_text(color = "black", size = 6),
                  legend.text = element_text(color = "black", size = 6),
                  legend.title = element_text(color = "black", size = 6),
                  plot.title = element_text(color = "black", size = 6)))

nonB_motif_levels = c("all", "APR", "DR", "STR", "IR", "TRI", "G4", "Z") 

####Read data####

list.files("results/phase2_including_repetative")

df = readr::read_delim("results/phase2_including_repetative/nonB_coverage_versus_read.csv")
df
unique(df$Window)
unique(df$Experiment)
unique(df$Label)
unique(df$nonB)

####Winsorize data####


#split by experiment ("WGS", "0mM", "20mM", "40mM", or "combined")
l_df = df %>%
  group_by(Label) %>%
  group_split()

#split by nonB and Q95 winsorize
length(l_df)
l_df_out = {}
for (i in 1:length(l_df)){
  print(i)
  l_df_i = l_df[[i]] %>% group_by(nonB) %>% group_split()
  l_df_i_out = {}
  for (j in 1:length(l_df_i)){
    df_j = data.frame(l_df_i[[j]])
    Q95_motif = quantile(df_j$motif_base_density, 0.95)
    Q95_reads = quantile(df_j$CPM, 0.95)
    Norm_motif_base_density = df_j$motif_base_density/Q95_motif
    Norm_motif_base_density[which(Norm_motif_base_density >= 1)] = 1
    Norm_CPM = df_j$CPM/Q95_reads
    Norm_CPM[which(Norm_CPM >= 1)] = 1
    df_j$Norm_motif_base_density = Norm_motif_base_density
    df_j$Norm_CPM = Norm_CPM
    l_df_i_out[[j]] = df_j
  }
  l_df_out[[i]] = bind_rows(l_df_i_out)
}
df_out = bind_rows(l_df_out)

head(df_out)

df_out$`non-B motif` = factor(df_out$nonB, levels = nonB_motif_levels)

####Plot data####

P = df_out %>%
  ggplot(aes(x = Norm_motif_base_density, y = Norm_CPM, color = blacklist_fraction)) +
  geom_point(alpha = 0.1, size = 0.1) +
  scale_color_viridis() +
  facet_grid(cols = vars(`non-B motif`), rows = vars(`Label`)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/Jacobs_plots/Figure_js4025.4.svg", P, width = 7, height = 6)
