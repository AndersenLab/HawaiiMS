#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(DT)
library(cowplot)
library(FSA)
library(GGally)
library(ggplotify)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

# set palettes
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
####################################
# tukey box plots and stats        #
####################################
# Kruskal-Wallis non-parametric tests on distinct collections, i.e. if a c-label had three c. elegans and 2 C. briggsae isolated from it we are counting it once for C. elegaans and once for C. briggsae.
stat_df <- data1 %>%
  dplyr::filter(pcr_positive == 1) %>%
  dplyr::mutate(species_family = ifelse(species_id %in% c("Panagrolaimus sp.",
                                                      "Oscheius sp.",
                                                      "Chabertia ovina", 
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
                                                      NA), 
                                        "PCR +", species_id)) %>%
  dplyr::filter(species_family != "PCR +") %>%
  dplyr::distinct(c_label, species_id, .keep_all = T) %>%
  dplyr::select(c_label, s_label, species_family, substrate_temperature, substrate_moisture, ambient_humidity, 
                ambient_temperature, altitude) %>%
  dplyr::mutate(species_family = forcats::as_factor(species_family),
                species_family = forcats::fct_relevel(species_family,
                                                      "C. elegans",
                                                      "C. oiwi",
                                                      "C. tropicalis",
                                                      "C. briggsae")) %>%
  tidyr::gather(env_par, value, substrate_temperature, substrate_moisture, ambient_humidity, ambient_temperature, altitude)%>%
  dplyr::group_by(env_par) %>%
  dplyr::mutate(KM_pvalue = kruskal.test(value ~ species_family)[[3]])

# Kruskal–Wallis test then multiple comparisions test using Dunn's test with pvalues adjusted with Bonferroni method
Dunn_list <- list()

for (e in 1:length(unique(stat_df$env_par))){
  KM_df <- stat_df %>%
    dplyr::filter(env_par == (unique(stat_df$env_par)[e]))
  
  D_test <- dunnTest(KM_df$value ~ KM_df$species_family, method = "bonferroni") 
  Dunn_list[[unique(stat_df$env_par)[e]]] <- D_test
}

# Plot data as box plots for each environmental parameter
plot_atemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "ambient_temperature")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(species_palette)) +
  geom_jitter(aes(x = species_family, y = value, fill = species_family), size = 1.5, width = .25, shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  labs(fill = "Species", x = "", y = "Ambient temperature (°C)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_stemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "substrate_temperature")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(species_palette)) +
  geom_jitter(aes(x = species_family, y = value, fill = species_family), size = 1.5, width = .25, shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  labs(fill = "Species", x = "", y = "Substrate temperature (°C)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_ahum <- ggplot(data = stat_df %>% dplyr::filter(env_par == "ambient_humidity")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(species_palette)) +
  geom_jitter(aes(x = species_family, y = value, fill = species_family), size = 1.5, width = .25, shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  labs(fill = "Species", x = "", y = "Ambient humidity (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_smoist <- ggplot(data = stat_df %>% dplyr::filter(env_par == "substrate_moisture")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(species_palette)) +
  geom_jitter(aes(x = species_family, y = value, fill = species_family), size = 1.5, width = .25, shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  labs(fill = "Species", x = "", y = "Substrate moisture (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_elev <- ggplot(data = stat_df %>% dplyr::filter(env_par == "altitude")) +
  geom_boxplot(aes(x = species_family, y = value, fill = species_family), outlier.color = NA) +
  scale_fill_manual(values=c(species_palette)) +
  geom_jitter(aes(x = species_family, y = value, fill = species_family), size = 1.5, width = .25, shape = 21) +  
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  labs(fill = "Species", x = "", y = "Elevation (m)") +
  theme(axis.title.x= element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

# create nice legend
legend <- get_legend(ggplot(stat_df %>% dplyr::filter(env_par == "altitude")) +
  geom_bar(aes(x = species_family, fill = species_family)) +
  scale_fill_manual(values=c(species_palette)) +
  labs(fill = "") +
  theme(legend.position = "bottom"))

####################################
# Correlation matrix and gg pairs  #
####################################

# Create correlation matrix dataframe
env_corr <- data1 %>%
  dplyr::filter(species_id %in% c("C. elegans", "C. oiwi", "C. tropicalis", "C. briggsae")) %>%
  dplyr::distinct(c_label, species_id, .keep_all = T) %>%
  dplyr::select(species_id, s_label, substrate_temperature, substrate_moisture, ambient_humidity, 
                ambient_temperature, altitude) %>%
  tibble::column_to_rownames ('s_label')

# define labels for correlation matrix plot
env_pars_labels <- c(NA, expression("substrate temperature"~degree~"C"), "substrate moisture (%)", "ambient humidity (%)", expression("ambient temperature"~degree~"C"), "altitude (m)")

# generate person correlation matrix plot
figure3F <- env_corr_matrix <- ggcorr(env_corr, low = "#440154", mid = "#218F8B", high = "#FDE725",
                          label_size=3, label_color='white', label = T, label_round = 2,
                          size = 3, color = "black",
                          layout.exp = 1) +
  theme(axis.line=element_blank())

# put all of figure 3 together
Figure3 <- cowplot::plot_grid(plot_atemp, plot_stemp, plot_ahum, plot_smoist, plot_elev, figure3F, legend, labels = c("A", "B", "C", "D", "E", "F", "", ""), ncol = 2, nrow = 4, align = "hv", axis = "l", rel_heights = c(1,1,1,.05))
ggsave('plots/Figure3.pdf', width = 7.5, height = 7.5, useDingbats = F)

# mean temps for species
par_means <- stat_df %>%
  dplyr::group_by(env_par, species_family) %>%
  dplyr::mutate(means = mean(value))

#save data for submission
#Figure 3A-E Environmental parameters for all Caenorhabditis collections
# stat_df %>%
#   readr::write_csv('data/elife_files/fig3-data1.csv')
# #Figure 3F pearson correlation matrix for environmental parameters
# as.data.frame(env_corr_matrix[1]) %>%
#   dplyr::select(-data.label) %>%
#   dplyr::rename(x = data.x, y = data.y, pearson_correlation = data.coefficient) %>%
#   dplyr::arrange(x) %>%
#   readr::write_csv('data/elife_files/fig3-data2.csv')
