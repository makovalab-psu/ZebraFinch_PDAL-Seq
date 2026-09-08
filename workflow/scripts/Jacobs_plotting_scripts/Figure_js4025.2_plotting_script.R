library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

####Global variables####

theme_set(theme_classic(base_size = 2) +
            theme(text = element_text(color = "black", size = 2),
                  axis.text = element_text(color = "black", size = 2),
                  axis.title = element_blank(),
                  strip.text = element_text(color = "black", size = 2),
                  legend.text = element_text(color = "black", size = 2),
                  legend.title = element_text(color = "black", size = 2),
                  plot.title = element_text(color = "black", size = 2)))

####Read data####

#Chromosomes
list.files("")
Exclude = c("chrMT", "rDNA_morph_1", "rDNA_morph_2", "rDNA_morph_3")
df_chr = readr::read_delim("resources/genomes/bTaeGut7v0.4_mat_Z_MT_rDNA.chrom.sizes", col_names = c("chr", "length")) %>% filter(!chr %in% Exclude)
chr_levels = unique(df_chr$chr)
df_chr$chr = factor(df_chr$chr, levels = chr_levels)

#Centromeres
list.files("resources/centromere")
df_cen = readr::read_delim("resources/centromere/bTaeGut7v0.4_MT_rDNA.matZ.CEN.bed", col_names = c("chr", "start", "end"))  %>% filter(!chr %in% Exclude)
df_cen$chr = factor(df_cen$chr, levels = chr_levels)

#PDAL-Seq data
list.files()
df_PDALSeq = readr::read_delim("data/phase2/experiment_window_coverage/All_reads/10000_nucleotides/Tguttata_CFS414_PDALSeq.bed", col_names = c("chr", "start", "end", "reads")) %>% filter(!chr %in% Exclude)
#Clip PDAL-Seq data to Q99
range(df_PDALSeq$reads)
quantile(df_PDALSeq$reads, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
Q99 = quantile(df_PDALSeq$reads, 0.99)
df_PDALSeq$Norm_PDALSeq = df_PDALSeq$reads/Q99
df_PDALSeq$Norm_PDALSeq[which(df_PDALSeq$Norm_PDALSeq >= 1)] = 1
df_PDALSeq$chr = factor(df_PDALSeq$chr, levels = chr_levels)

#Control data
list.files("data/phase2/experiment_window_coverage/All_reads/10000_nucleotides/")
df_WGS = readr::read_delim("data/phase2/experiment_window_coverage/All_reads/10000_nucleotides/Tguttata_CFS414_WGS.bed", col_names = c("chr", "start", "end", "reads")) %>% filter(!chr %in% Exclude)
#Clip PDAL-Seq data to Q99
range(df_WGS$reads)
quantile(df_WGS$reads, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
Q99 = quantile(df_WGS$reads, 0.99)
df_WGS$Norm_WGS = df_WGS$reads/Q99
df_WGS$Norm_WGS[which(df_WGS$Norm_WGS >= 1)] = 1
df_WGS$chr = factor(df_WGS$chr, levels = chr_levels)

#Blacklist
list.files("data/phase2/blacklist/merged.bed")
df_bl = readr::read_delim("data/phase2/blacklist/merged.bed", col_names = c("chr", "start", "end"))  %>% filter(!chr %in% Exclude)
df_bl$chr = factor(df_bl$chr, levels = chr_levels)

#non-B DNA motifs

#APR
list.files("data/phase2/nonb_density/10000_nucleotides")
df_APR = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/APR.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_APR$chr = factor(df_APR$chr, levels = chr_levels)

#DR
list.files("data/phase2/nonb_density/10000_nucleotides")
df_DR = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/DR.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_DR$chr = factor(df_DR$chr, levels = chr_levels)

#STR
list.files("data/phase2/nonb_density/10000_nucleotides")
df_STR = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/STR.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_STR$chr = factor(df_STR$chr, levels = chr_levels)

#IR
list.files("data/phase2/nonb_density/10000_nucleotides")
df_IR = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/IR.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_IR$chr = factor(df_IR$chr, levels = chr_levels)

#TRI
list.files("data/phase2/nonb_density/10000_nucleotides")
df_TRI = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/TRI.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_TRI$chr = factor(df_TRI$chr, levels = chr_levels)

#G4
list.files("data/phase2/nonb_density/10000_nucleotides")
df_G4 = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/G4.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_G4$chr = factor(df_G4$chr, levels = chr_levels)

#Z
list.files("data/phase2/nonb_density/10000_nucleotides")
df_Z = readr::read_delim("data/phase2/nonb_density/10000_nucleotides/Z.bg", col_names = c("chr", "start", "end", "annotations", "length", "total", "fraction"))  %>% filter(!chr %in% Exclude)
df_Z$chr = factor(df_Z$chr, levels = chr_levels)

####Plot data####

labels_df <- data.frame(
  chr = factor("chr1_mat", levels = chr_levels),  # or whichever facet you want them on
  x = c(-0.75, -0.3, 1, 3, 4.2, 4.6, 5.0, 5.4, 5.8, 6.2, 6.6),
  y = 100000,
  label = c("cen","blacklist","PDAL-Seq","WGS","APR","DR","STR","IR","TRI","G4","Z")
)

P = ggplot() +
      geom_rect(data = df_chr, mapping = aes(xmin = -0.9, xmax = -0.6, ymax = 0, ymin = -length), color = "grey", fill = "white", linewidth = 0.25) +
      geom_rect(data = df_cen, mapping = aes(xmin = -0.95, xmax = -0.55, ymax = -start, ymin = -end), color = "dimgrey", fill = "dimgrey", linewidth = 0.25) +
      geom_rect(data = df_bl, mapping = aes(xmin = -0.5, xmax = -0.1, ymax = -start, ymin = -end), fill = "black") +
      geom_rect(data = df_PDALSeq, mapping = aes(xmin = 0, xmax =2*Norm_PDALSeq, ymax = -start, ymin = -end), fill = "darkred") +
      geom_rect(data = df_WGS, mapping = aes(xmin = 2.1, xmax = 2.1 + 2*Norm_WGS, ymax = -start, ymin = -end), fill = "dimgrey") +
      geom_rect(data = df_APR, mapping = aes(xmin = 4.0, xmax = 4.4, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_DR, mapping = aes(xmin = 4.4, xmax = 4.8, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_STR, mapping = aes(xmin = 4.8, xmax = 5.2, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_IR, mapping = aes(xmin = 5.2, xmax = 5.6, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_TRI, mapping = aes(xmin = 5.6, xmax = 6.0, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_G4, mapping = aes(xmin = 6.0, xmax = 6.4, ymax = -start, ymin = -end, fill = fraction)) +
      geom_rect(data = df_Z, mapping = aes(xmin = 6.4, xmax = 6.8, ymax = -start, ymin = -end, fill = fraction)) +
      geom_text(data = labels_df, aes(x, y, label = label), angle = 45, inherit.aes = FALSE, size = 1) +
      scale_fill_gradient(low = "white", high = "red", name = "non-B motif coverage") +
      facet_wrap(~chr, scales = "free_y") +
      theme_void() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom") 

ggsave("plots/Jacobs_plots/Figure_js4025.2.png", P, width = 7, height = 10)
