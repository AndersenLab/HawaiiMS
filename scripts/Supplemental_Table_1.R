#!/usr/bin/env Rscript
#Load necessary packages
library(janitor)
library(tidyverse)
library(rcompanion)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

#======================#
# Supplemental Table 1 #
#======================#

# Species by island
Supp_table1 <- data1 %>%
  dplyr::mutate(plot_type = ifelse(worms_on_sample %in% c("No","?"), "No Worm",
                                   ifelse(is.na(species_id), "Not genotyped",
                                          ifelse(pcr_positive == 0, "PCR -",
                                                 ifelse(species_id == "Unknown", "Not genotyped",
                                                        species_id))))) %>%
  dplyr::distinct(c_label, plot_type, .keep_all = T) %>% 
  dplyr::mutate(plot_type = forcats::as_factor(plot_type),
                plot_type = forcats::fct_relevel(plot_type,
                                                 "C. elegans",
                                                 "C. oiwi",
                                                 "C. tropicalis",
                                                 "C. kamaaina",
                                                 "C. briggsae",
                                                 "Panagrolaimus sp.",
                                                 "Oscheius sp.",
                                                 "Teratorhabditis sp.",
                                                 "Rhabditis terricola",
                                                 "Choriorhabditis sp.",
                                                 "Mesorhabditis sp.",
                                                 "Chabertia ovina",
                                                 "Heterhabditis zealandica",
                                                 "PCR -",
                                                 "Not genotyped",
                                                 "No Worm")) %>%
  dplyr::arrange(desc(plot_type)) %>%
  dplyr::rename(collection_category = plot_type) %>%
  dplyr::group_by(collection_category, island) %>%
  dplyr::summarize(worm_isolates=as.integer(n())) %>%
  tidyr::spread(island, worm_isolates, fill=0) %>%
  janitor::adorn_totals('row') %>%
  janitor::adorn_totals('col') %>%
  readr::write_tsv("data/elife_files/supp-table1.tsv")

fractions <- Supp_table1 %>%
  dplyr::mutate(all_samples = 2263,
                perc_total = (Total/all_samples)*100)

# Island enrichment analysis
# shape data for analysis for each Caenorhabditis species
FE_isl_ce <- Supp_table1 %>%
  dplyr::rename(species = collection_category) %>%
  dplyr::filter(species == "C. elegans") %>%
  tidyr::gather(island, count, -Total, -species) %>%
  dplyr::select(island, ce = count) %>%
  dplyr::mutate(isl_total = case_when(
    island == "Big Island" ~ 701,
    island == "Maui" ~ 466,
    island == "Molokai" ~ 86,
    island == "Oahu" ~ 101,
    island == "Kauai" ~ 909),
    no_ce = isl_total-ce) %>%
  dplyr::select(island, no_ce, ce) %>%
  column_to_rownames(var = "island") %>%
  as.matrix(.)

FE_isl_cb <- Supp_table1 %>%
  dplyr::rename(species = collection_category) %>%
  dplyr::filter(species == "C. briggsae") %>%
  tidyr::gather(island, count, -Total, -species) %>%
  dplyr::select(island, cb = count) %>%
  dplyr::mutate(isl_total = case_when(
    island == "Big Island" ~ 701,
    island == "Maui" ~ 466,
    island == "Molokai" ~ 86,
    island == "Oahu" ~ 101,
    island == "Kauai" ~ 909),
    no_cb = isl_total-cb) %>%
  dplyr::select(island, no_cb, cb) %>%
  column_to_rownames(var = "island") %>%
  as.matrix(.)

FE_isl_ct <- Supp_table1 %>%
  dplyr::rename(species = collection_category) %>%
  dplyr::filter(species == "C. tropicalis") %>%
  tidyr::gather(island, count, -Total, -species) %>%
  dplyr::select(island, ct = count) %>%
  dplyr::mutate(isl_total = case_when(
    island == "Big Island" ~ 701,
    island == "Maui" ~ 466,
    island == "Molokai" ~ 86,
    island == "Oahu" ~ 101,
    island == "Kauai" ~ 909),
    no_ct = isl_total-ct) %>%
  dplyr::select(island, no_ct, ct) %>%
  column_to_rownames(var = "island") %>%
  as.matrix(.)

FE_isl_co <- Supp_table1 %>%
  dplyr::rename(species = collection_category) %>%
  dplyr::filter(species == "C. oiwi") %>%
  tidyr::gather(island, count, -Total, -species) %>%
  dplyr::select(island, co = count) %>%
  dplyr::mutate(isl_total = case_when(
    island == "Big Island" ~ 701,
    island == "Maui" ~ 466,
    island == "Molokai" ~ 86,
    island == "Oahu" ~ 101,
    island == "Kauai" ~ 909),
    no_co = isl_total-co) %>%
  dplyr::select(island, no_co, co) %>%
  column_to_rownames(var = "island") %>%
  as.matrix(.)

# Post-hoc pairwise Fisher tests with pairwise.table
options(scipen = 999)
pairwiseNominalIndependence(FE_isl_ce,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)
pairwiseNominalIndependence(FE_isl_cb,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)
pairwiseNominalIndependence(FE_isl_ct,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)
pairwiseNominalIndependence(FE_isl_co,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

#==================================#
# Number of worm isolates / sample #
#==================================#

data1 %>%
  dplyr::group_by(c_label) %>%
  dplyr::filter(!is.na(s_label)) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::summarize(m = mean(n))

# Species by sample

data1 %>% 
  dplyr::group_by(c_label) %>%
  dplyr::filter(!is.na(spp_id)) %>%
  dplyr::distinct(c_label, spp_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n >= 2) %>%
  tidyr::nest(spp_id, .key = "species") %>%
  dplyr::mutate(has_unknown = purrr::map_lgl(species,  ~ "Unknown" %in% .x$spp_id)) %>%
  dplyr::mutate(n_species = purrr::map_int(species, ~ length(.x$spp_id))) %>% 
  dplyr::select(has_unknown) %>% table()


#=======================#
# Isotypes by substrate #
#=======================#

data1 %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(isotype, substrate) %>%
  dplyr::group_by(substrate) %>%
  dplyr::summarize(n=n())

#=======================#
# Isotypes by landscape #
#=======================#

data1 %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::distinct(isotype, landscape) %>%
  dplyr::group_by(landscape) %>%
  dplyr::summarize(n=n())
