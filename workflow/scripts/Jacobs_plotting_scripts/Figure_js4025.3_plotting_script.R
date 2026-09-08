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

####Read data###

df = readr::read_delim("results/phase2/nonB_correlation_coefficient.csv")
head(df)

unique(df$Experiment)
unique(df$nonB)
df$`non-B motif` = factor(df$nonB, levels = nonB_motif_levels)
unique(df$Label)
unique(df$Library)
unique(df$Metric)
unique(df$Signal)
unique(df$Blacklist_filtered)

P = df %>%
  filter(Blacklist_filtered == FALSE) %>%
  filter(Metric == "count") %>%
  filter(Signal == "CPM") %>%
  ggplot(aes(x = `non-B motif`, y = Spearman)) +
    geom_bar(stat = "summary", fun = mean, fill = "dimgrey") +
    geom_point(mapping = aes(color = Label, shape = Label), size = 1) +
    facet_wrap(~Window + Library, nrow = 1) +
    scale_color_manual(values = viridis(5, option = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/Jacobs_plots/Figure_js4025.3.svg", P, width = 6, height = 2.5, units = "in")


