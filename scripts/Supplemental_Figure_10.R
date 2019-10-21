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

pop_labels <- c("A" = "HI Low",
                "C" = "C",
                "E" = "HI Inv.",
                "F" = "HI Div.",
                "H" = "Vol.")
# % sharing broken up by chromosome for global C population and Hawaiian poulations A, E, F, H.
Supp_Figure10 <- ggplot(admix_sharing) +
  aes(x = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")), y = max_haplotype_shared, fill = factor(pop_assignment, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))) +
  scale_fill_manual(values=c(ancestry.colours)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = ifelse(admix_sharing$Hawaiian == T, 21, 25)) +
  scale_x_discrete(labels = pop_labels) +
  facet_wrap(~chromosome, scales = "free") +
  labs(y = "Fraction most common global haplotype", x = "") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_text_repel(aes(label=ifelse(max_haplotype_shared > 0.3 & pop_assignment != "C" & Hawaiian == T,
                                   as.character(isotype),'')),hjust=-2,vjust=.2, size = 2)
Supp_Figure10
ggsave(paste("plots/Supplemental Figure 10.pdf"), width = 7.5, height = 7.5, useDingbats = F)

# # Write data for supplemental figure 10
# admix_sharing %>%
#   dplyr::rename(fraction_swept_haplotype_on_chromosome = max_haplotype_shared) %>%
#   dplyr::select(-start, -stop, -haplotype, -isotype_has_swept_haplotype, -sat, -plotpoint, -segment, -cvalue, -hap_length, -color, -chrom_haplotype_sum,
#                 -swept_haplotype, -swept_haplotype_name, -isotypes_w_haplotype, -isotype_swept_haplotype_length, -max_swept_haplotype_length,
#                 -filtered_swept_haplotype_len, -filtered_sweep_len, -filtered_sweep_ratio, -is_swept, -max_pop_frac) %>%
#   readr::write_csv(., "data/elife_files/supp-fig10-data1.csv")
