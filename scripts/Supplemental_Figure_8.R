#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(ggtern)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

# load admix info for full 276_set
admix <- data.table::fread("data/ADMIXTURE/full/BEST_K/K11_Processed_Ancestry.tsv") %>%
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

#generate ternery plots
plotEH <- admix %>%
   dplyr::filter(!pop_assignment %in% c("A", "B", "D", "F", "G", "I", "J", "K")) %>%
   dplyr::select(isotype, pop_assignment, E, C, H)

Supp_Figure_8 <-  ggtern(data=plotEH,aes(E, C, H,fill=pop_assignment)) +
    scale_fill_manual(values = ancestry.colours) +
    theme_classic() +
    #theme_rgbw() +
    geom_point(size = 3, shape = 21) + 
    tern_limit(T = 1.1, L = 1.1, R = 1.1) +
    labs(x="HI Div.",y="Global C",z="Vol.",title="")
Supp_Figure_8

ggsave(paste("plots/Supplemental Figure 8.pdf"), width = 3.75, height = 3.75, useDingbats=FALSE)

# # Write data for ternary plot including all isotypes in ancestral populations C, E, H
# plotEH %>%
# readr::write_csv(., "data/elife_files/supp-fig8-data1.csv")