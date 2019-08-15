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
library(rcompanion)

# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "K"="mediumpurple4","L"= "orange","M"= "maroon","N"= "yellow3","O"= "brown4", 
                      "P"="yellow4", "Q"="sienna4", "R"="chocolate", "S"="gray19")

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

# assign Hawaii isotypes
hi_only_samples <- read.csv(file = "data/fulcrum/hawaii_isotypes.csv") 

#load admixture proportions for LD8
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

# read in climate data (Katie 2018) (these data are for isotype ref strains only)
clim <- data.table::fread("data/NOAA_12mo_climate_variables.csv") %>%
  dplyr::select(-V1) %>%
  dplyr::rename(strain = isotype) %>%
  dplyr::filter(strain %in% admix$isotype)

# get wild isolate info for all 276 isotypes.
isotype_data <- data.table::fread("data/WI_strain_list.csv") %>%
  dplyr::filter(isotype %in% c(as.character(admix$isotype))) %>%
  dplyr::select(isotype, strain, latitude, longitude, isolation_date)
  
# combine admix data with climate data (filter out strains with station distance greter than 100km)
pop_clim <- full_join(admix, clim  %>% dplyr::rename(isotype = strain)) %>% # %>% dplyr::filter(station_distance <= 100)
  tidyr::spread(trait, value)

# plot sd for populations 
Supp_Figure10 <- ggplot(pop_clim) +
  aes(x = factor(pop_assignment, levels = c("B", "C", "D", "G", "I", "J", "K", "A", "E", "F", "H")), y = as.double(sd.temp), fill = pop_assignment) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .25, shape = ifelse(pop_clim$Hawaiian == TRUE, 21, 24)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(y = "12-month temperature SD (°C)", x = "Ancestral population", fill = "") +
  theme_bw()
Supp_Figure10

# save climate data plot
ggsave('plots/Supplemental Figure 10.pdf', width = 7.5, height = 5, useDingbats=FALSE)

# Write data file
# pop_clim %>%
#   dplyr::select(isotype, max_pop_frac, pop_assignment, A, B, C, D, E, F, G, H, I, K, K, Hawaiian,
#                 latitude, longitude, isolation_date, start_date, end_date, nearest_station, station_distance, station_lat, station_lon,
#                 mean.temp, sd.temp) %>%
# readr::write_csv(., "data/elife_files/supp-fig10-data1.csv")
