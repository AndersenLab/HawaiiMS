#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)
library(DT)
library(FSA)
library(scales)
library(ggrepel)
library(magick)
library(ggplotify)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

#load admixture proportions for K=11 LD8
admix <- data.table::fread("data/ADMIXTURE/full/BEST_K/K11_Processed_Ancestry.tsv",header = T) %>%
  dplyr::rename(isotype = samples) %>%
  tidyr::gather(pop, frac_pop, - isotype) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(max_pop_frac = max(frac_pop)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(pop, max_pop_frac) %>%
  dplyr::mutate(isotype = factor(isotype)) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(pop_assignment = ifelse(max_pop_frac == frac_pop, pop, NA)) %>%
  dplyr::arrange(isotype, pop_assignment) %>%
  tidyr::fill(pop_assignment) %>%
  tidyr::spread(pop, frac_pop) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Hawaiian = ifelse(isotype %in% hi_only_samples$isotype, "TRUE", "FALSE"))

# Load processed haps
load("data/HAPLOTYPE/haplotype_plot_df.Rda")

# join admixture info
hap_admix_df <- dplyr::left_join(plot_df, admix)

# Calculate % sharing of max haplotype among all populations in K=11
admix_sharing_all <- hap_admix_df  %>% 
  distinct(isotype, chromosome, .keep_all= TRUE) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(genome_swept_hap_length = sum(isotype_swept_haplotype_length),
                genome_max_length = sum(max_swept_haplotype_length),
                genome_frac_swept = genome_swept_hap_length/genome_max_length) %>%
  dplyr::ungroup()

# Calculate % sharing of max haplotype between global C population and Hawaiian poulations A, E, F, H. No filtering of sweep.
admix_sharing <- hap_admix_df  %>%
  dplyr::filter(pop_assignment %in% c("C","E", "A", "H", "F")) %>%
  distinct(isotype, chromosome, .keep_all= TRUE) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(genome_swept_hap_length = sum(isotype_swept_haplotype_length),
                genome_max_length = sum(max_swept_haplotype_length),
                genome_frac_swept = genome_swept_hap_length/genome_max_length) %>%
  dplyr::ungroup()

# Plot figure 7B
Figure7B <- ggplot(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T)) +
  aes(x = factor(pop_assignment, levels = c("A", "C", "E", "F", "H")), y = genome_frac_swept, fill = factor(pop_assignment, levels = c("A", "C", "E", "F", "H"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, shape = ifelse(admix_sharing %>% dplyr::distinct(isotype, .keep_all=T) %>% .$Hawaiian == T, 21, 25)) +
  labs(y = "Fraction of genome swept haplotype", x = "Ancestral population") +
  theme_bw() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=ifelse(genome_frac_swept > 0.2 & pop_assignment != "C" & Hawaiian == T,as.character(isotype),'')),hjust=-1,vjust=.2, size = 3)
Figure7B

# Analyze differences among populations 
options(scipen=999)
admix_sharing_anlaysis <- admix_sharing %>%
  distinct(isotype, .keep_all= TRUE) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")))

D_test <- dunnTest(data =admix_sharing_anlaysis, genome_frac_swept~pop_assignment, method = "bonferroni") 
D_test[[2]]

# Setup to plot Hawaii isotype haplotypes and global pop C haplotypes by admix population using filtered sweep ratio
mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

hap_admix_df_ordered <- hap_admix_df %>%
  dplyr::filter(Hawaiian == "TRUE" | pop_assignment == "C") %>%
  dplyr::arrange(desc(pop_assignment), plotpoint) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("C", "E",  "A", "H", "F")))

plotpoints <- hap_admix_df_ordered %>%
  dplyr::distinct(pop_assignment, isotype) %>%
  dplyr::mutate(pop_assignment = factor(pop_assignment, levels = c("C", "E",  "A", "H", "F"))) %>%
  dplyr::arrange(pop_assignment) %>%
  dplyr::mutate(plotpoint_hi = row_number())

hap_admix_df_ordered <- dplyr::left_join(hap_admix_df_ordered, plotpoints) 

strain_labels_hi <- hap_admix_df_ordered %>% 
  dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "C") %>%
  dplyr::select(isotype, plotpoint_hi)

# Plot figure 7A
Figure7A <- ggplot(hap_admix_df_ordered %>% dplyr::filter(Hawaiian == "TRUE"| pop_assignment == "C"),
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint_hi - 0.5, ymax = plotpoint_hi + 0.5,
           fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = unique(strain_labels_hi$plotpoint_hi),
                     labels = unique(strain_labels_hi$isotype),
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(factor(pop_assignment, levels=c("F", "H", "A", "E", "C"))~chromosome, scales="free", space="free") +
  theme(legend.position="none")
Figure7A

# Load Treemix migration and residual plot generated from Supplemental_figure_8.R script
Figure7C_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_5_migrations.pdf", density = 300))
Figure7C_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_5_migrations_RESIDUALS.pdf", density = 300))

Figure7C <- cowplot::plot_grid(Figure7C_1, Figure7C_2, ncol=2, rel_widths = c(1,.66))

# Make full Figure 7 plot
Figure7B_C <- cowplot::plot_grid(Figure7B, Figure7C, labels = c("B","C"), rel_widths = c(.5, 1), align = "hv")
Figure7 <- cowplot::plot_grid(Figure7A, Figure7B_C, labels = c("A", ""), nrow = 2, rel_heights = c(1, .4))
Figure7

ggsave('plots/Figure7.pdf', height = 15, width = 10, useDingbats = F)

# Write data for Figure 7
# #haplotype data for 7A
# hap_admix_df %>%
# readr::write_csv(., "data/elife_files/fig7-data1.csv")
# #Sweep data for 7B
# admix_sharing %>%
#   dplyr::distinct(isotype, .keep_all=T) %>%
#   dplyr::select(isotype, is_swept, pop_assignment, A, B, C, D, E, F, G, H, I, J, K,
#                 Hawaiian, genome_swept_hap_length, genome_max_length, genome_frac_swept) %>%
#   readr::write_csv(., "data/elife_files/fig7-data2.csv")
# #data for 7C is contained in .P and .Q files in data/ADMIXTURE/full/BEST_K directory

