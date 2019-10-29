#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(ggtree)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# load in wi strain info for lat long coordinates for map
df <- data.table::fread("data/WI_strain_list.csv")
strain_colors <- dplyr::filter(df, state== "Hawaii")

tree_275 <- ape::read.tree(glue::glue("data/tree/249-hawaii_genome.raxml.bestTree"))

tree_pt<- ggtree(tree_275,
                 branch.length="rate",layout="equal_angle")+xlim(NA,0.3)

hawaii_strains <- df%>%
  dplyr::mutate(strain_loc = ifelse(strain %in% strain_colors$isotype, "Hawaii", "Not")) %>%
  dplyr::filter(release %in% c("20160408","20170531") | strain %in% strain_colors$isotype) %>%
  dplyr::filter(reference_strain == 1) %>%
  dplyr::select(isotype, lat = latitude, long = longitude, strain_loc) %>%
  dplyr::mutate(color = isotype ) %>%
  as.data.frame()


colored_tree <- tree_pt %<+% hawaii_strains + 
  geom_tiplab(aes(color = strain_loc)) +
  scale_color_manual(values = c("#D7263D", "gray40")) + theme_tree2()

branch_strains <- list(CONNECT = dplyr::filter(hawaii_strains, strain_loc == "Hawaii") %>% dplyr::pull(isotype))

tree_pt_h <- ggtree::groupOTU(tree_275, branch_strains)

# define colors
highlight_color <- "#D7263D"
background_color <- "#000F08"

#plot tree with Hawaiian isotypes colored red
hi_tree <- ggtree(tree_pt_h,
       branch.length="rate", 
       layout="equal_angle",
       aes(color=group)) + 
  geom_tiplab() +
  scale_color_manual(values=c(background_color, highlight_color), 
                     name = "Presence of TALT", 
                     breaks=c("0", "TALT"),
                     labels=c("FALSE", "TRUE")) + 
  theme(legend.position="right")+
  # geom_hilight(1, "steelblue") +
  theme_tree2() 
hi_tree

ggsave(hi_tree, filename = "plots/Supplemental Figure 4.pdf", height = 5, width = 7.5)