#!/usr/bin/env Rscript
#Load necessary packages
library(pophelper)
library(tidyverse)
library(ggthemes)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# hawaii strains
strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# get list of isotype names from samples.txt
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))

# plot panel A, K by CV summary
k_summary <- data.table::fread("data/ADMIXTURE/full/CV_Summary/admix_replicates_CV.tsv",header = T) 

ksum_plot <- ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  theme_bw()+
  labs(x = "K")

# generate K summary plot - Supplemental figure XX
admix_plots <- list()
for(kpops in 1:length(grep(".Q", list.files("data/ADMIXTURE/full/BEST_K/"), value = T))){
  K <- as.numeric(strsplit(grep(".Q", list.files("data/ADMIXTURE/full/BEST_K/"), value = T)[kpops], split = "\\.")[[1]][4])
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE/full/BEST_K/"))
  qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE/full/BEST_K/",qfile_name))[[1]]
  # add pop names
  colnames(qfile) <- LETTERS[1:K]
  
  qfile <- qfile %>%
    dplyr::mutate(samples = sample_names)
  
  write.table(qfile, file = glue::glue("data/ADMIXTURE/full/BEST_K/K{K}_Processed_Ancestry.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # make long and determin order of plotting
  long_admix_pops <- qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  # establish plot order of strains based on anc pop and max fraction
  plot_order <- long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  admix_plots[[kpops]] <- long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_blank(),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(angle = 90, vjust = .5),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))
  
  ggsave(admix_plots[[kpops]]  +
           theme(axis.text.x=element_text(angle = 90)), filename = glue::glue("plots/K{K}_names.pdf"), height = 4, width = 28)
  
  # extract representative strains from each ancesteral population for generating neighbor-net
  if(!exists("representative_K_strains")){
    representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999, !samples %in% names(strain_islands)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K)
  } else {
    representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999, !samples %in% names(strain_islands)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K) %>%
      dplyr::bind_rows(representative_K_strains, .)
  }
  
}

# Write representative strains file for use in Figure 5
write.table(representative_K_strains, 
            file = "data/ADMIXTURE/full/BEST_K/Hawaii_plus_RepStrains_by_K.tsv",
            quote = F,
            col.names = T, 
            row.names = F, 
            sep = "\t")

# make panel B
admixture_plots <- cowplot::plot_grid(admix_plots[[1]],
                   admix_plots[[2]],
                   admix_plots[[3]],
                   admix_plots[[4]],
                   admix_plots[[5]],
                   admix_plots[[6]],
                   ncol = 1)

# make final figure
ksummary_plot <- cowplot::plot_grid(ksum_plot,
                   admixture_plots,
                   ncol = 2,
                   labels = c("A", "B"),
                   rel_widths = c(0.5, 1))

# save
ggsave(ksummary_plot, filename = "plots/Supplemental Figure 5.pdf", height = 8, width = 12, useDingbats=FALSE)

# save supplemental data for manuscript
# # Write admixture files
# qfile7 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.7.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 7)
# qfile8 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.8.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 8)
# qfile9 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.9.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 9)
# qfile10 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.10.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 10)
# qfile11 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.11.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 11)
# qfile12 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.12.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 12)
# qfile13 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.13.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 13)
# qfile14 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.14.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 14)
# qfile15 <- pophelper::readQ(files = "data/ADMIXTURE/full/BEST_K/LD_0.8_MAF_0.004.15.Q")[[1]] %>%
#   dplyr::mutate(samples = sample_names) %>%
#   dplyr::mutate(K_value = 15)
# full_admixture <- full_join(qfile15, qfile14) %>%
#   full_join(., qfile13) %>%
#   full_join(., qfile12) %>%
#   full_join(., qfile11) %>%
#   full_join(., qfile10) %>%
#   full_join(., qfile9) %>%
#   full_join(., qfile8) %>%
#   full_join(., qfile7) %>%
#   dplyr::rename(A = Cluster1,
#                 B = Cluster2,
#                 C = Cluster3,
#                 D = Cluster4,
#                 E = Cluster5,
#                 F = Cluster6,
#                 G = Cluster7,
#                 H = Cluster8,
#                 I = Cluster9,
#                 J = Cluster10,
#                 K = Cluster11,
#                 L = Cluster12,
#                 M = Cluster13,
#                 N = Cluster14,
#                 O = Cluster15) %>%
#   readr::write_csv(., "data/elife_files/supp-fig5-data2.csv")

#write data for cross validation error plot
#k_summary %>% readr::write_csv(., "data/elife_files/supp-fig5-data1.csv")