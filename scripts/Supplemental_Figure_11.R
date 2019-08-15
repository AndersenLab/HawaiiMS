#!/usr/bin/env Rscript
#Load necessary packages
library(pophelper)
library(tidyverse)
library(ggthemes)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#set VCF path
VCF_PATH="data/ANNOTATE_VCF/LD_0.8.vcf.gz"

# hawaii strains
strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# pie_chart define a color pallette with maximal contrast
ancestry.colours.pie <- c("perc_a"="gold2", "perc_b"="plum4","perc_c"= "darkorange1", 
                          "perc_d"="lightskyblue2", "perc_e"="firebrick","perc_f"= "burlywood3",
                          "perc_g"="gray51",
                          "perc_h" = "springgreen4",
                          "perc_i" = "lightpink2",
                          "perc_j" = "deepskyblue4",
                          "perc_k" = "mediumpurple4",
                          "perc_l" = "orange")

# strain names
# get sample information
sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples_HI_ONLY.txt", header = F) %>% dplyr::pull(V1)) 

# plot panel A, K by CV summary
k_summary <- data.table::fread("data/ADMIXTURE/Hawaii/CV_Summary/admix_replicates_CV.tsv",header = T) 

ksum_plot <- ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  theme_bw()+
  labs(x = "K")

# generate K summary plot - Supplemental figure XX
admix_plots <- list()
for(kpops in 1:length(grep(".Q", list.files("data/ADMIXTURE/Hawaii/BEST_K/"), value = T))){
  K <- as.numeric(strsplit(grep(".Q", list.files("data/ADMIXTURE/Hawaii/BEST_K/"), value = T)[kpops], split = "\\.")[[1]][4])
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE/Hawaii/BEST_K/"))
  qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE/Hawaii/BEST_K/",qfile_name))[[1]]
  # add pop names
  colnames(qfile) <- LETTERS[1:K]
  
  qfile <- qfile %>%
    dplyr::mutate(samples = sample_names)
  
  write.table(qfile, file = glue::glue("data/ADMIXTURE/Hawaii/BEST_K/K{K}_Processed_Ancestry.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
  
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
           theme(axis.text.x=element_text(angle = 90)), filename = glue::glue("plots/{K}_Hawaii_names.pdf"), height = 4, width = 28)
  ggsave(admix_plots[[kpops]]  +
           theme(axis.text.x=element_text(angle = 90)), filename = glue::glue("plots/{K}_Hawaii_names.pdf"), height = 4, width = 28,dpi=300)
}

# make panel B
admixture_plots <- cowplot::plot_grid(admix_plots[[1]],
                                      admix_plots[[2]],
                                      admix_plots[[3]],
                                      ncol = 1)

# make final figure
ksummary_plot <- cowplot::plot_grid(ksum_plot,
                                    admixture_plots,
                                    ncol = 2,
                                    labels = c("A", "B"),
                                    rel_widths = c(0.5, 1))

# get big legend
admix_legend <- cowplot::plot_grid(cowplot::get_legend(admix_plots[[3]]))

# save
ggsave(ksummary_plot, filename = "plots/Supplmental Figure 11.pdf", height = 2, width = 7.5, useDingbats = F)

# generate input for running SplitsTree
strain_vector <- paste(sort(names(strain_islands)), sep = ",", collapse = ",")
system(glue::glue("bcftools view -s {strain_vector} {VCF_PATH} -Oz -o data/ANNOTATE_VCF/Hawaii.vcf.gz"))
# Generate nexus file, .nexus file is used for SplitsTree
system(glue::glue("python scripts/vcf2phylip.py -i data/ANNOTATE_VCF/Hawaii.vcf.gz -m 43 --fasta --nexus --nexus-binary"))

# Write data for submission
# k3 <- admix_plots[[1]]
# k3df <- as.data.frame(k3[1]) %>%
#   dplyr::rename(isotypes = data.samples, cluster = data.cluster, frac_cluster = data.frac_cluster, max_frac = data.max_frac) %>%
#   dplyr::select(-data.ordered_samples) %>%
#   dplyr::mutate(K = 3)
# k4 <- admix_plots[[1]]
# k4df <- as.data.frame(k4[1]) %>%
#   dplyr::rename(isotypes = data.samples, cluster = data.cluster, frac_cluster = data.frac_cluster, max_frac = data.max_frac) %>%
#   dplyr::select(-data.ordered_samples) %>%
#   dplyr::mutate(K = 4)
# k5 <- admix_plots[[1]]
# k5df <- as.data.frame(k5[1]) %>%
#   dplyr::rename(isotypes = data.samples, cluster = data.cluster, frac_cluster = data.frac_cluster, max_frac = data.max_frac) %>%
#   dplyr::select(-data.ordered_samples) %>%
#   dplyr::mutate(K = 5)
# #bind all dataframes
# rbind(k3df, k4df, k5df) %>%
#   readr::write_csv('data/elife_files/supp-fig11-data1.csv')
# #ADMIXTURE subsampling data for ploting cross validation (CV) error with Hawaiian isolates only
# k_summary %>%
#   readr::write_csv('data/elife_files/supp-fig11-data2.csv')