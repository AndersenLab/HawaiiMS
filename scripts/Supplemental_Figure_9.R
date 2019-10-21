#!/usr/bin/env Rscript
#Load necessary packages
#treeMix v. 1.13 is required to run this script
library(tidyverse)
library(ggplotify)
library(magick)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#define outgroup as XZ1516
outgroup_strain <- "ECA191"

sample_names <- sort(data.table::fread("data/ANNOTATE_VCF/samples.txt", header = F) %>% dplyr::pull(V1))

# define K to use for analysis
K = 11

# load P files from ADMIXTURE analysis
pfile_name <- grep(pattern = glue::glue("{K}\\.P$"), value = T, x = list.files("data/ADMIXTURE/full/BEST_K/"))
pfile <- pophelper::readQ(files = paste0("data/ADMIXTURE/full/BEST_K/",pfile_name))[[1]]

# label P file rownames and colnames
colnames(pfile) <- LETTERS[1:K]

# Make treemix input
treemix_input <- apply(pfile, MARGIN = c(1,2), function(x){
  f <- round(x*length(sample_names),digits = 0)
  paste(length(sample_names)-f, f, sep = ",")
})

# find outgroup population
# load Q files
qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files("data/ADMIXTURE/full/BEST_K/"))
qfile <- pophelper::readQ(files = paste0("data/ADMIXTURE/full/BEST_K/",qfile_name))[[1]]
# add pop names
colnames(qfile) <- LETTERS[1:K]
qfile$strain <- sample_names

outgroup_population <- qfile%>%
  dplyr::filter(strain == outgroup_strain)%>%
  tidyr::gather(Population, Frequency, -strain)%>%
  dplyr::filter(Frequency == max(Frequency))

write.table(x = treemix_input,
            file = glue::glue("data/ADMIXTURE/full/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
            quote = F,
            col.names = T,
            row.names = F)

treemix_input_name <- strsplit(glue::glue("data/ADMIXTURE/full/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                               split = "/")[[1]][4]

system(glue::glue("gzip data/ADMIXTURE/full/{treemix_input_name}"))

outgroup_pop <- outgroup_population$Population[1]

output_name <- strsplit(strsplit(glue::glue("data/ADMIXTURE/full/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                                 split = "/")[[1]][4], split = "\\.txt")[[1]][1]

#Running treemix with migration up to five events
for(m in 1:5){
  system(glue::glue("treemix -i data/ADMIXTURE/full/{treemix_input_name}.gz -o data/ADMIXTURE/full/{output_name}_{m} -root {outgroup_pop} -k 500 -m {m} -se"))
}

# load functions for ploting TreeMix output
source("scripts/PLOT_TREEMIX.R")

#Plot all treemix plots and residual plots
for(migration in gsub("\\.llik","",grep(glue::glue("K-{K}"), grep("llik",list.files("data/ADMIXTURE/full"),value = T), value = T))){
  
  system(paste0(paste0("printf \"", paste(LETTERS[1:K], collapse = '\n')), '\n\" > data/ADMIXTURE/full/poporder'))
  
  pdf(glue::glue("plots/TREEMIX_phylogeny_{migration}_migrations.pdf"))
  plot_tree(stem =  glue::glue("data/ADMIXTURE/full/{migration}"))
  dev.off()
  
  pdf(glue::glue("plots/TREEMIX_phylogeny_{migration}_migrations_RESIDUALS.pdf"))
  plot_resid(stem = glue::glue("data/ADMIXTURE/full/{migration}"), "data/ADMIXTURE/full/poporder")
  dev.off()
  # 
}

# Plot Supplemental Figure 9
Figure9A_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_migrations.pdf", density = 300))
Figure9A_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_migrations_RESIDUALS.pdf", density = 300))

Figure9B_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_1_migrations.pdf", density = 300))
Figure9B_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_1_migrations_RESIDUALS.pdf", density = 300))

Figure9C_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_2_migrations.pdf", density = 300))
Figure9C_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_2_migrations_RESIDUALS.pdf", density = 300))

Figure9D_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_3_migrations.pdf", density = 300))
Figure9D_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_3_migrations_RESIDUALS.pdf", density = 300))

Figure9E_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_4_migrations.pdf", density = 300))
Figure9E_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_4_migrations_RESIDUALS.pdf", density = 300))

Figure9F_1 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_5_migrations.pdf", density = 300))
Figure9F_2 <- as.ggplot(image_read("plots/TREEMIX_phylogeny_K-11_Outgroup=H=TREEMIX_input_5_migrations_RESIDUALS.pdf", density = 300))

Supp_Figure_9 <- cowplot::plot_grid(Figure9A_1, Figure9A_2, Figure9B_1, Figure9B_2, 
                   Figure9C_1, Figure9C_2, Figure9D_1, Figure9D_2,
                   Figure9E_1, Figure9E_2, Figure9F_1, Figure9F_2,
                   ncol = 4, nrow = 3, labels = c("A", "", "B", "",
                                                  "C", "", "D", "",
                                                  "E", "", "F", ""))
Supp_Figure_9

ggsave('plots/Supplemental Figure 9.pdf', height = 5, width = 7.5, useDingbats = F)