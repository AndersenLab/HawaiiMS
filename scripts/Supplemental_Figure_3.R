#!/usr/bin/env Rscript
# load necessary packages
library(tidyverse)
library(ggnetwork)
library(intergraph)
library(igraph)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

###############################
# Gridsect Functions          #
###############################
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
                      "Fruit/nut/vegetable"= 24,
                      "Rotting flower"= 22,
                      "Fungus"= 23,
                      "Isopod"= 25)

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

# Supplemental Figure 3A: Gridsects 1 and 3 by isotypes
isotype = data1 %>%
  dplyr::mutate(isotype = ifelse(is.na(isotype), "No isotype", isotype))

gridsect_list <- lapply(c(1:15, 17:21), function(x){ 
  if (x == 21) {
    plot_gridsect(5, isotype) + 
      geom_text(size= 8)
  } else {
    plot_gridsect(x, isotype) +
      geom_point(aes(fill = isotype,  size = size_c,  shape = substrate), stroke = 0.3) +
      scale_fill_manual("Isotypes", values = c(isotype_palette), drop = FALSE) +
      scale_shape_manual("Substrate", values = c(substrate_shapes), drop = FALSE)
  }
})

#Plot Supplemental Figure 3A
legend_isotypes <- get_legend(ggplot(isotype %>% dplyr::filter(isotype %in% names(isotype_palette))) +
                                geom_bar(aes(x = strain, fill = isotype)) +
                                scale_fill_manual(values=c(isotype_palette)) +
                                labs(fill = "isotypes"))

legend_substrates <- get_legend(ggplot(isotype %>% dplyr::filter(!is.na(substrate)) %>% dplyr::mutate(substrate = factor(substrate, levels = names(substrate_shapes)))) +
                                  geom_point(aes(x = substrate_temperature, y = substrate_temperature, shape = substrate)) +
                                  scale_shape_manual("substrates", values = c(substrate_shapes), drop = FALSE) +
                                  labs(shape = "substrates"))

Supp_fig_3A <- cowplot::plot_grid(gridsect_list[[1]],
                                  gridsect_list[[3]],
                                  legend_isotypes,
                                  legend_substrates,
                                  hjust = -0.09,
                                  labels = c("A", "B"),
                                  ncol = 4)

###############################
# Supplemental Figure 3B      #
###############################

# function to build input for network from relational observations
# df contains a grouping column, relationship column, observation number column
# ~ ~ #  grouping column - is used to group by and look for relationships
# ~ ~ #  relationship column - end up being the nodes of the network (can have multiple relationship columns)
# repeat_observations - count observations of the same individual

generate_relational_network_df <- function(df,  
                                           group_col, 
                                           relationship_col, 
                                           repeat_observations) { 
  
  # process input variables 
  define_range <- df %>%
    dplyr::group_by(!!group_col) %>% 
    dplyr::mutate(observation = 1:n())%>%
    dplyr::ungroup() %>%
    dplyr::mutate(dimensions = length(unique((!!relationship_col)))) %>% 
    dplyr::select( .,
                   groupings = !!group_col, #group_col
                   rel = !!relationship_col, #relationship_col
                   observation, 
                   dimensions)
  
  if(!repeat_observations){
    define_range <- define_range %>%
      dplyr::distinct(groupings, rel, .keep_all =T)
  }
  
  # initialize matrix to fill with relationships
  empty_net <- matrix(nrow=length(unique(define_range$rel)),
                      ncol=length(unique(define_range$rel)),
                      dimnames = list(unique(define_range$rel),
                                      unique(define_range$rel)))
  
  # iterate through unique observations in relationship variable 
  for(i in unique(define_range$rel)){
    
    subet_groups <- dplyr::filter(define_range, rel == i)%>%
      dplyr::distinct(groupings, rel,.keep_all=T)%>%
      dplyr::mutate(id = paste(groupings, rel, observation, sep = "_"))
    
    found_together <- dplyr::filter(define_range, groupings %in% subet_groups$groupings)%>%
      dplyr::mutate(id = paste(groupings, rel, observation, sep = "_"))%>%
      dplyr::filter(!(id %in% subet_groups$id))
    
    for(j in found_together$rel){
      if(is.na(empty_net[j,i])){
        empty_net[j,i] <- 1
      } 
      else {
        empty_net[j,i] <- empty_net[j,i] + 1
      }
    }
  }
  
  empty_net[is.na(empty_net)] <-0
  
  full_net <- empty_net
  
  return(full_net)
}

data1_isotypes <- data1 %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::select(isotype, c_label, s_label)%>%
  dplyr::group_by(c_label)%>%
  dplyr::summarise(isotype_ct = length(unique(isotype)))%>%
  dplyr::filter(isotype_ct > 1)

isotypes_found_together <- data1%>%
  dplyr::filter(c_label %in% data1_isotypes$c_label)%>%
  dplyr::select(isotype, c_label, s_label)%>%
  dplyr::filter(!is.na(isotype))%>%
  dplyr::group_by(c_label)

function_net <- generate_relational_network_df(isotypes_found_together, 
                                               quo(c_label), 
                                               quo(isotype), 
                                               repeat_observations = F)

net=graph.adjacency(function_net,
                    mode="undirected",
                    weighted=TRUE,
                    diag=FALSE)

# plot supplemental figure 3B
Supp_figure_3b <- ggplot(ggnetwork::ggnetwork(net), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(color = "black", size = 8) +
  geom_nodelabel_repel(aes( label = vertex.names ),
                       fontface = "bold.italic", box.padding = unit(1, "lines"))+
  theme_blank()+
  geom_edgetext(aes(label = weight), color = "grey25")

###################################
# Plot Full Supplemental Figure 3 #
###################################
# plot Supplemental Figure 3
Supp_figure_3 <- cowplot::plot_grid(Supp_fig_3A,
                                    Supp_figure_3b,
                                    labels = c("", "C"),
                                    ncol = 2)

cowplot::ggsave(filename = "plots/Supplemental Figure 3.pdf", plot = Supp_figure_3, height = 5, width = 37.5) 

###################################
# Describe gridsects and isotypes #
###################################
# number of isotypes on a c_plate (2 collections are excluded because no isotypes were recovered from them due to extinction before sequencing)
num_isotypes_per_sample <- data1 %>%
  dplyr::filter(!is.na(s_label)) %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(c_label, isotype, .keep_all = TRUE) %>%
  dplyr::select(isotype, species_id, c_label, s_label, substrate, latitude, longitude) %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(distinct_isotypes = n()) %>%
  dplyr::distinct(c_label, .keep_all=TRUE)

# Write data files for submission
# Collection data for gridsects with C. elegans isotypes 
# isotype %>%
#   dplyr::filter(gridsect_number %in% c(1,3)) %>%
#   readr::write_csv('data/elife_files/supp-fig3-data1.csv')
# Collection data for Instances where two distinct isotypes were isolated from the same sample
# data1 %>%
#   dplyr::filter(s_label %in% isotypes_found_together$s_label) %>%
#   readr::write_csv('data/elife_files/supp-fig3-data2.csv')