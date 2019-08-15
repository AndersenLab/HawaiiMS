#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(grid)
library(gridExtra)
library(imager)
library(jpeg)
library(pals)
library(rcompanion)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

# Color palettes
species_palette <- c("C. elegans" = "#BE0032", #7
                     "C. oiwi" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     "Panagrolaimus sp." = "#C2B280", #8
                     "Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1", #6 
                     "Other PCR +" = "#008856", #10
                     "Caenorhabditis" = "#C2B280", #8
                     "PCR -" = "#848482", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "No Worm" = "#222222", #2
                     "multiple" = "#FFFFFF") 

substrate_palette <- c("Leaf litter" = "#E68FAC",
                       "Fruit/nut/veg" = "#0067A5",
                       "Flower" = "#DCD300",
                       "Fungus" = "#604E97",
                       "Compost" = "#F6A600",
                       "Other" = "#B3446C")

allPalette <- c("Caenorhabditis" = "#C2B280",
                "Other PCR +" = "#008856",
                "PCR -" = "#848482",
                "Not genotyped" = "#F2F3F4",
                "Tracks only" = "#b3b3b3",
                "No nematode" ="#222222")

#######################################
# Figure 2A                           #
#######################################
#Figure 2A define rh postitive c-labels
rh_positive_c_labels <- data1 %>%
  dplyr::filter(pcr_positive == 1) %>%
  dplyr::distinct(c_label) 

#Figure 2A df for all collections broken into braod collection categories (i.e., No worm, tracks, not genptyped, pcr-, pcr+)
worms <- data1 %>%
  dplyr::filter(worms_on_sample != "?") %>%
  dplyr::mutate(plot_type = ifelse(worms_on_sample == "No", "No nematode",
                                         ifelse(worms_on_sample == "Tracks", "Tracks only",
                                                ifelse(worms_on_sample == "Yes" & is.na(pcr_positive), "Not genotyped",
                                                  ifelse(worms_on_sample == "Yes" & pcr_positive == 0, "PCR -",
                                                        ifelse(worms_on_sample == "Yes" & pcr_positive == 1 & !(species_id %in% c("C. elegans",
                                                                                                                                "C. briggsae",
                                                                                                                                "C. oiwi",
                                                                                                                                "C. tropicalis",
                                                                                                                                "C. kamaaina")), "Other PCR +", "Caenorhabditis")))))) %>%
  dplyr::mutate(plot_type = factor(plot_type, levels = c("Caenorhabditis",
                                                         "Other PCR +",
                                                         "PCR -",
                                                         "Not genotyped",
                                                         "Tracks only",
                                                         "No nematode"))) %>%
  dplyr::arrange(plot_type) %>%
  dplyr::distinct(c_label, .keep_all = TRUE) %>%
  dplyr::group_by(plot_type, fixed_substrate) %>%
  dplyr::mutate(worm_per_substrate = n()) %>% 
  dplyr::distinct(plot_type, fixed_substrate, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_substrates = sum(worm_per_substrate), perc_worm_sub = worm_per_substrate / total_substrates * 100) %>%
  dplyr::arrange(total_substrates) %>%
  dplyr::select(fixed_substrate, worm_per_substrate, total_substrates, perc_worm_sub, plot_type) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(plot_type = factor(plot_type, levels = rev(names(allPalette)))) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(c("Leaf litter","Fruit/nut/veg", "Flower", "Fungus",
                                                                         "Invertebrate", "Other"))))

# plot for all collections
Fig2A <- ggplot(worms) +
  geom_bar(stat = "identity", aes(x = factor(fixed_substrate), y = perc_worm_sub, fill = plot_type), colour = "black") + 
  scale_fill_manual(values=c(allPalette))+
  coord_flip() + 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of all collections") +
  geom_text(aes(x=fixed_substrate, y=113, label=paste0("n=", total_substrates)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 120))

#######################################
# Figure 2B                           #
#######################################
newdf <- data1 %>%
  dplyr::filter(species_id %in% c("C. elegans",
                              "C. briggsae",
                              "C. oiwi",
                              "C. tropicalis")) %>%
  dplyr::select(c_label, species_id, fixed_substrate) %>%
  dplyr::mutate(species_family = factor(species_id, levels = c("C. elegans", "C. oiwi", "C. tropicalis", "C. briggsae"))) %>%
  dplyr::arrange(species_family) %>%
  dplyr::distinct(c_label, .keep_all=T) %>%
  dplyr::group_by(species_family, fixed_substrate) %>%
  dplyr::mutate(worm_per_substrate = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(fixed_substrate) %>%
  dplyr::mutate(total_substrates = n(),
                perc_worm_sub = worm_per_substrate / total_substrates * 100) %>%
  dplyr::arrange(total_substrates) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(species_family, fixed_substrate, .keep_all = T) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(names(substrate_palette)))) %>%
  dplyr::mutate(species_family = factor(species_family, levels = rev(c("C. elegans", "C. oiwi", "C. tropicalis", "C. briggsae"))))

# Fig2B plot for rhabditida positive collections
newdf_2 <- newdf %>%
  dplyr::mutate(tot_sub = case_when(
    fixed_substrate == "Leaf litter" ~ 1480,
    fixed_substrate == "Fruit/nut/veg" ~ 333,
    fixed_substrate == "Flower" ~ 202,
    fixed_substrate == "Fungus" ~ 121,
    fixed_substrate == "Other" ~ 73),
    perc_worm_sub2 = (worm_per_substrate/tot_sub)*100) %>%
  dplyr::select(species_family, fixed_substrate, perc_worm_sub2)

species_family <- c("C. tropicalis", "C. tropicalis", "C. oiwi", "C. oiwi", "C. oiwi", "C. elegans")
fixed_substrate <- c("Fungus", "Other", "Other", "Leaf litter", "Fungus", "Fungus")
temp <- data.frame(species_family, fixed_substrate)

# join to fill out dataframe
newdf_2 <- full_join(newdf_2, temp) %>%
  dplyr::mutate(perc_worm_sub2 = ifelse(is.na(perc_worm_sub2), 0, perc_worm_sub2)) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = rev(names(substrate_palette)))) %>%
  dplyr::mutate(species_family = factor(species_family, levels = rev(names(species_palette))))

Fig2B_v2 <- ggplot(data = newdf_2) +
  geom_bar(stat = "identity", position = "dodge", aes(x = factor(fixed_substrate), y = perc_worm_sub2, fill = species_family), colour = "black") +
  scale_fill_manual(values=c(species_palette)) +
  coord_flip() + 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of all collections") +
  guides(fill = guide_legend(reverse = T)) +
  scale_y_continuous(breaks = c(2.5, 5.0, 7.5, 10), limits = c(0, 10.5))

##################################
# Plot 2A and 2B together        #
##################################
# Put fig 2A and fig 2B_v2 together
figure2AB <- cowplot::plot_grid(Fig2A, Fig2B_v2, labels = c("A","B"), ncol=1, align = "v")
#ggsave('plots/Fig2AB_substrate_panel_v3.pdf', height = 6.59, width = 5)

##################################
# C Image gallery                #
##################################
# create function to load and rasterize images from a set of urls in the cso or df
rplots <- function(x) {
  ggplot() + annotation_raster(raster = imager::load.image(x), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
}

# Collection plots with intersting images. 
mdf <- data1 %>%
  dplyr::filter(pcr_positive == 1) %>%
  dplyr::distinct(c_label, .keep_all = T) %>%
  dplyr::arrange(species_id, substrate) %>%
  dplyr::mutate(., plot_labels = paste0(c_label,"\n",substrate,"\n",species_id)) %>%
  dplyr::select(c_label, substrate, species_id, photo_url_thumb, plot_labels )

spp <- rep(c("C. elegans", "C. briggsae", "C. tropicalis", "C. oiwi"), each = 3)
type <- rep(c("Leaf litter", "Flower", "Fruit/nut/veg"), 4)
c_label <- c("C-3133", "C-0736", "C-2830", NA,"C-0257", "C-0072",
             "C-1083", "C-2909", "C-2906", "C-2735", "C-0846", "C-2845")

# Make df with dummy grey box image for plotting blanks
sub_df <- as.data.frame(cbind(spp, type, c_label)) %>%
  dplyr::left_join(., mdf) %>%
  dplyr::mutate(photo_url_thumb = ifelse(is.na(photo_url_thumb),"https://drive.google.com/uc?export=download&id=10D7zXBStCPcdHHk_UE2r6QBZYqHiBNL6", photo_url_thumb))

pl <- lapply(sub_df$photo_url_thumb, rplots)
margin = theme(plot.margin = unit(c(.5,.5,.5,.5), "mm"))
gl <- gridExtra::grid.arrange(grobs = lapply(pl, "+", margin), ncol = 3)

z <- textGrob("")
a <- textGrob("C. elegans")
b <- textGrob("C. briggsae")
c <- textGrob("C. tropicalis")
d<- textGrob("C. oiwi")

sa <- textGrob("Leaf litter")
sb <- textGrob("Flower")
sc <- textGrob("Fruit/nut/veg")

yaxis <- grid.arrange(z,a,d,c,b, ncol = 1)
xaxis <- grid.arrange(sa,sb,sc, nrow = 1)
lay <- rbind(c(1,2,2,2),
             c(1,3,3,3),
             c(1,3,3,3),
             c(1,3,3,3),
             c(1,3,3,3))

figure2C <- grid.arrange(yaxis, xaxis, gl, layout_matrix = lay, ncol = 4, nrow =5)

##################################
# Plot and Save Figure 2         #
##################################
Figure_2 <- cowplot::plot_grid(figure2AB, figure2C, labels = c("","C"), rel_widths = c(1, .6067))
ggsave('plots/Figure2.pdf', plot = Figure_2, height = 5, width = 7.5)
# the height and width ratio accounts for image dimesions and column numbers to avoid transforming image dimensions. Height = 1.648 times Width
#ggsave('plots/20190219_Reduced_Substrate_images_test.pdf', plot = full, height = 8.24, width = 5)

##################################
# Occurence frequency statistics #
##################################
# shape data for analysis at Caenorhabditis level
chi_sub_test <- worms %>%
  dplyr::filter(plot_type == "Caenorhabditis") %>%
  dplyr::mutate(no_c = total_substrates-worm_per_substrate,
                c = worm_per_substrate) %>%
  dplyr::select(fixed_substrate, no_c, c) %>%
  column_to_rownames(var = "fixed_substrate") %>%
  as.matrix(.)

# show p-values without scientific notation
options(scipen=999)

# full chi squared test to see if there are differernces 
chisq.test(chi_sub_test, simulate.p.value = TRUE)

# Post-hoc pairwise Fisher tests with pairwise.table
pairwiseNominalIndependence(chi_sub_test,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

##########################################################
# Analysis within Caenorhabditis for substrate enrichment#
##########################################################
# shape data for analysis for each Caenorhabditis species
FE_test_ce <- newdf %>%
  dplyr::filter(species_family == "C. elegans") %>%
  dplyr::mutate(tot_sub = case_when(
                                    fixed_substrate == "Leaf litter" ~ 1480,
                                    fixed_substrate == "Fruit/nut/veg" ~ 333,
                                    fixed_substrate == "Flower" ~ 202,
                                    fixed_substrate == "Fungus" ~ 121,
                                    fixed_substrate == "Other" ~ 73),
                no_ce = tot_sub-worm_per_substrate,
                ce = worm_per_substrate) %>%
  dplyr::select(fixed_substrate, no_ce, ce) %>%
  column_to_rownames(var = "fixed_substrate") %>%
  as.matrix(.)

FE_test_co <- newdf %>%
  dplyr::filter(species_family == "C. oiwi") %>%
  dplyr::mutate(tot_sub = case_when(
    fixed_substrate == "Leaf litter" ~ 1480,
    fixed_substrate == "Fruit/nut/veg" ~ 333,
    fixed_substrate == "Flower" ~ 202,
    fixed_substrate == "Fungus" ~ 121,
    fixed_substrate == "Other" ~ 73),
    no_co = tot_sub-worm_per_substrate,
    co = worm_per_substrate) %>%
  dplyr::select(fixed_substrate, no_co, co) %>%
  column_to_rownames(var = "fixed_substrate") %>%
  as.matrix(.)

FE_test_ct <- newdf %>%
  dplyr::filter(species_family == "C. tropicalis") %>%
  dplyr::mutate(tot_sub = case_when(
    fixed_substrate == "Leaf litter" ~ 1480,
    fixed_substrate == "Fruit/nut/veg" ~ 333,
    fixed_substrate == "Flower" ~ 202,
    fixed_substrate == "Fungus" ~ 121,
    fixed_substrate == "Other" ~ 73),
    no_ct = tot_sub-worm_per_substrate,
    ct = worm_per_substrate) %>%
  dplyr::select(fixed_substrate, no_ct, ct) %>%
  column_to_rownames(var = "fixed_substrate") %>%
  as.matrix(.)

FE_test_cb <- newdf %>%
  dplyr::filter(species_family == "C. briggsae") %>%
  dplyr::mutate(tot_sub = case_when(
    fixed_substrate == "Leaf litter" ~ 1480,
    fixed_substrate == "Fruit/nut/veg" ~ 333,
    fixed_substrate == "Flower" ~ 202,
    fixed_substrate == "Fungus" ~ 121,
    fixed_substrate == "Other" ~ 73),
    no_cb = tot_sub-worm_per_substrate,
    cb = worm_per_substrate) %>%
  dplyr::select(fixed_substrate, no_cb, cb) %>%
  column_to_rownames(var = "fixed_substrate") %>%
  as.matrix(.)

# Run analysis
# Post-hoc pairwise Fisher tests with pairwise.table
pairwiseNominalIndependence(FE_test_ce,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

pairwiseNominalIndependence(FE_test_co,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

pairwiseNominalIndependence(FE_test_ct,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

pairwiseNominalIndependence(FE_test_cb,
                            fisher = TRUE,
                            gtest  = FALSE,
                            chisq  = FALSE,
                            method = "bonferroni", simulate.p.value = TRUE)

#### Spearman's rho statistic is used to estimate a rank-based measure of association between substrate and population size
cor_test_df <- data1 %>%
  dplyr::filter(species_id %in% c("C. elegans", "C. tropicalis", "C. briggsae",  "C. oiwi"))  %>%
  dplyr::distinct(c_label, species_id, .keep_all = T) %>%
  dplyr::select(species_id, fixed_substrate, approximate_number_of_worms) %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = c("Leaf litter", "Flower", "Fruit/nut/veg", "Fungus",  "Other")),
                fixed_substrate1 = as.numeric(fixed_substrate),
                approx_pop = factor(approximate_number_of_worms, levels = c("Very Few (1-3)", "Few (4-10)", "Some (11-25)", "Proliferating (25+)")),
                approx_pop1 = as.numeric(approx_pop),
                approx_pop2 = ifelse(approximate_number_of_worms == "Proliferating (25+)", 2,1))

# setup test levels for cb and ct and do test inside of dataframe
cor_test_df_flvsll_test <- cor_test_df %>%
  dplyr::filter(fixed_substrate %in% c("Flower", "Leaf litter")) %>%
  dplyr::group_by(species_id) %>%
  dplyr::do(broom::tidy(cor.test(.$fixed_substrate1, .$approx_pop2, alternative = "g")))
cor_test_df_flvsll_test

# setup test levels for co and do test inside of dataframe
cor_test_df_fnvvsfl_test <- cor_test_df %>%
  dplyr::mutate(fixed_substrate = factor(fixed_substrate, levels = c("Fruit/nut/veg", "Flower", "Leaf litter", "Fungus",  "Other")), # reorder for cor test between fruit and flower 
                fixed_substrate1 = as.numeric(fixed_substrate)) %>%
  dplyr::filter(fixed_substrate %in% c("Fruit/nut/veg", "Flower"))%>%
  dplyr::group_by(species_id) %>%
  dplyr::do(broom::tidy(cor.test(.$fixed_substrate1, .$approx_pop2, alternative = "g")))
cor_test_df_fnvvsfl_test

# Write data for submission
# #2A
# worms %>%
#   readr::write_csv('data/elife_files/fig2-data1.csv')
# #2B
# newdf_2 %>%
#   readr::write_csv('data/elife_files/fig2-data2.csv')