#!/usr/bin/env Rscript
# load necessary packages
library(tidyverse)
library(cowplot)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

#====================#
# Gridsect Functions #
#====================#
isotype_palette <- c("ECA760" = "#654522" ,
                     "ECA778" =  "#8DB600",
                     "ECA768" = "#882D17" , 
                     "ECA812" = "#DCD300" ,
                     "ECA730" = "#B3446C" , 
                     "ECA712" =  "#F6A600", 
                     "ECA777" = "#604E97" , 
                     "ECA807" = "#F99379")

species_palette <- c("C. elegans" = "#BE0032", #7
                     "C. sp. 53" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     "Panagrolaimus sp." = "#C2B280", #8
                     "Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1", #6 
                     "Other PCR +" = "#008856", #10
                     "PCR -" = "#848482", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "No Worm" = "#222222") #2

substrate_shapes <- c("Leaf litter" = 21,
                      "Fruit"= 24,
                      "Flower"= 22,
                      "Fungus"= 23,
                      "Invertebrate"= 25)

adjust_x <- function(x, n, rn) {
  
  if (n == 1) {
    return(x)
  } else if ((n == 2 && rn == 1) || (n == 3 && rn == 2) || (n == 4 && rn == 1)  || (n == 4 && rn == 3) || (n %in% c(6, 7) && rn %in% c(3,5)) ) {
    return(x - 0.09)
  } else if (n == 2 && rn == 2  || (n == 3 && rn == 3) || (n == 4 && rn == 2) || (n == 4 && rn == 4) || (n %in% c(6, 7) && rn %in% c(2,4)) ) {
    return(x + 0.09)
  } else if (n == 3 && rn >= 2) {
    return(x + 0.09)
  } else if (n == 5 && rn == 2) {
    return(x - 0.09)
  } else if (n == 5 && rn == 3) {
    return(x + 0.09)
  } else if (n == 5 && rn == 4) {
    return (x - 0.09)
  } else if (n == 5 && rn == 5) {
    return (x + 0.09)
  } else if (n == 7 && rn == 7) {
    return (x + 0.25)
  }
  return(x)
}

adjust_y <- function(y, n, rn) {
  if (n == 1) {
    return(y)
  } else if ( (n == 3 && rn == 1 ) || (n == 4 && rn <= 2) || (n == 5 && rn %in% c(4, 5)) || (n %in% c(6, 7) && rn %in% c(4,5))) {
    return(y - 0.09)
  } else if ( (n == 3 && rn >= 2 ) || (n == 4 && rn >= 3) || (n == 5 && rn %in% c(2, 3))  || (n %in% c(6, 7) && rn %in% c(2,3)) ) {
    return (y + 0.09)
  } else if ((n == 5 && rn == 1) || (n == 6 && rn == 1)) {
    return (y + 0.25)
  } else if ( (n %in% c(6, 7) && rn == 6)) {
    return (y - 0.25)
  }
  return(y)
}

set_size <- function(n) {
  if (n == 1) {
    return(5)
  } else {
    return(2.4)
  }
}


plot_gridsect <- function(gn, data1 = data1) {
  angles = list(A = 1,
                B = 2,
                C = 3,
                D = 4,
                E = 5,
                F = 6)
  
  mm <- data1 %>%
    dplyr::filter(gridsect_number == gn, !is.na(gridsect_number)) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(x = gridsect_radius * (sin( (angles[[gridsect_direction]] - 1) * (pi/3) )),
                  y = gridsect_radius * (cos( (angles[[gridsect_direction]] - 1) * (pi/3) )),
                  label = paste0(gridsect_direction, gridsect_radius)) %>%
    dplyr::group_by(gridsect_radius, gridsect_direction) %>%
    dplyr::mutate(n = n(),
                  rn = row_number()) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(x = adjust_x(x, n, rn),
                  y = adjust_y(y, n, rn),
                  size_c = set_size(n)) %>%
    dplyr::ungroup()
  
  ggplot(mm, aes(x = x, y = y, label = label)) +
    annotate("path",
             x=0+3.5*cos(seq(0,2*pi,length.out=100)),
             y=0+3.5*sin(seq(0,2*pi,length.out=100))) +
    theme_void() +
    theme(legend.position = "none") +
    scale_size_area(guide = "none")
}

#=======================#
# Supplemental Figure 1 # 
#=======================#
sp = data1 %>%
  dplyr::filter(!is.na(gridsect_number)) %>%
  #dplyr::mutate(spp_id = ifelse(is.na(spp_id), "Not genotyped", spp_id)) %>%
  dplyr::mutate(plot_type = ifelse(worms_on_sample %in% c("No","?"), "No Worm",
                                   ifelse(is.na(species_id), "Not genotyped",
                                          ifelse(pcr_positive == 0, "PCR -",
                                                 ifelse(species_id %in% c("Chabertia ovina",
                                                                      "Choriorhabditis cristata",
                                                                      "Choriorhabditis sp.",
                                                                      "Heterhabditis zealandica",
                                                                      "Mesorhabditis sp.",
                                                                      "no match",
                                                                      "C. kamaaina",
                                                                      "Rhabditis terricola",
                                                                      "Rhanditis tericola",
                                                                      "Teratorhabditis sp.",
                                                                      "Unknown",
                                                                      "unknown",
                                                                      NA), "Other PCR +", species_id)))))

gridsect_list_full <- lapply(c(1:15, 17:21), function(x){ 
  if (x == 21) {
    plot_gridsect(5, sp) + 
      geom_text(size= 8)
  } else {
    plot_gridsect(x, sp) +
      geom_point(aes(fill = plot_type,  size = size_c, shape = fixed_substrate), stroke = 0.3) +  # scale_shape_manual(values=c(3, 16, 17))+
      scale_fill_manual("Species", values = c(species_palette), drop = FALSE) +
      scale_shape_manual("Substrate", values = c(substrate_shapes), drop = FALSE)
  }
})

# Make legends for Supplemental Figure 1
legend_species <- get_legend(ggplot(sp %>% dplyr::filter(plot_type %in% names(species_palette)) %>%
                                      dplyr::mutate(plot_type = factor(plot_type, levels = names(species_palette)))) +
                               geom_bar(aes(x = strain, fill = plot_type)) +
                               scale_fill_manual(values=c(species_palette)) +
                               labs(fill = "species") +
                               theme(legend.position = "bottom"))
  

legend_substrates <- get_legend(ggplot(sp %>% dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = names(substrate_shapes)))) +
                                  geom_point(aes(x = substrate_temperature, y = substrate_temperature, shape = fixed_substrate)) +
                                  scale_shape_manual("substrates", values = c(substrate_shapes), drop = FALSE) +
                                  labs(shape = "substrates") +
                                  theme(legend.position = "bottom"))

Supp_Fig_1 <- cowplot::plot_grid(plot_grid(plotlist=gridsect_list_full, hjust = -1, label_size = 20,
                    nrow = 4, ncol=5, rel_heights = c(5,5,5,5), labels=c(1:15, 17:20)), legend_species, legend_substrates,
          rel_heights = c(20, 1, 1), nrow=3)

cowplot::ggsave("plots/Supplemental Figure 1.pdf", width = 28.56, height = 25.2)

# # Make data files for submission
# # Collection categories and substrate types for gridsect samples
# sp %>% dplyr::filter(plot_type %in% names(species_palette)) %>%
#   dplyr::mutate(plot_type = factor(plot_type, levels = names(species_palette))) %>%
#   dplyr::select(collection_category = plot_type, fixed_substrate, everything()) %>%
#   dplyr::arrange(collection_category) %>%
#   dplyr::select(-substrate) %>%
#   readr::write_csv('data/elife_files/supp-fig1-data1.csv')