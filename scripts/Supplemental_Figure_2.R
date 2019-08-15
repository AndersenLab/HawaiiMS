#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(cowplot)
library(ggnetwork)
library(intergraph)
library(igraph)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

data1_species <- data1 %>%
  dplyr::filter(!is.na(species_id), species_id != "no match", species_id != "Unknown", species_id != "unknown")%>%
  dplyr::select(species_id, c_label)%>%
  dplyr::group_by(c_label)%>%
  dplyr::summarise(species_ct = length(unique(species_id)))%>%
  dplyr::filter(species_ct > 1)

species_found_together <- data1 %>%
  dplyr::mutate(species_id = ifelse(species_id == "Heterhabditis zealandica", "Heterohabditis zealandica", species_id)) %>%
  dplyr::filter(c_label %in% data1_species$c_label)%>%
  dplyr::select(c_label, species_id)%>%
  dplyr::filter(!is.na(species_id), species_id != "Unknown")%>%
  dplyr::group_by(c_label)

species_found_together$species_id <- gsub(" ", 
                                      "\n",
                                      species_found_together$species_id)

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

function_net <- generate_relational_network_df(species_found_together, 
                                               quo(c_label), 
                                               quo(species_id), 
                                               repeat_observations = F)

net=graph.adjacency(function_net,
                    mode="undirected",
                    weighted=TRUE,
                    diag=FALSE)


ggplot(ggnetwork::ggnetwork(net), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(color = "black", size = 8) +
  geom_nodelabel_repel(aes( label = vertex.names ),
                       fontface = "bold.italic", box.padding = unit(1, "lines"))+
  theme_blank()+
  geom_edgetext(aes(label = weight), color = "grey25") 

ggsave('plots/Supplemental Figure 2.pdf', width = 7.5, height = 5)

# Write data for submission
# Instances where two distinct species were isolated from the same sample
# data1 %>%
#   dplyr::filter(c_label %in% species_found_together$c_label) %>%
#   dplyr::distinct(c_label, .keep_all = T) %>%
#   dplyr::select(-species_id, -s_label, -isotype, -strain, -pcr_positive) %>%
#   left_join(species_found_together, .) %>%
#   readr::write_csv('data/elife_files/supp-fig2-data1.csv')