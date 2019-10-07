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

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
data1 <- data.table::fread("data/Supplemental Data 1.tsv")

species_palette <- c("C. elegans" = "#BE0032", #7
                     "C. oiwi" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     "Panagrolaimus sp." = "#C2B280", #8
                     "Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1", #6 
                     "Other PCR +" = "#008856", #10
                     "PCR -" = "#848482", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "Tracks only" = "#b3b3b3", #manual
                     "No Nematode" = "#222222")  #2

species_palette <- c("C. elegans" = "#222222", #7
                     "C. oiwi" =  "#222222", #4 
                     "C. tropicalis" = "#222222", #5 
                     "Panagrolaimus sp." = "#222222", #8
                     "Oscheius sp." = "#222222", #3
                     "C. briggsae" = "#222222", #6 
                     "Other PCR +" = "#222222", #10
                     "PCR -" = "#222222", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "Tracks only" = "#BE0032", #manual
                     "No Nematode" = "#222222")  #2

island_palette <- c("Kauai" = "#E69F00",
                    "Oahu" = "#56B4E9",
                    "Molokai" = "#009E73",
                    "Maui" = "#F0E442",
                    "Big Island" = "#D55E00")

####################################################
#  A: Define Functions                             #
####################################################
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

islands = list(  "Kauai" = c(-159.830818,21.750571,-159.230003,22.350076),
  "Oahu" = c(-158.323116,21.112767,-157.623081,21.814254),
  "Molokai" = c(-157.3515,20.793,-156.6515,21.4956),
  "Maui" = c(-156.745977,20.405495,-155.942774,21.207099),
  "Big Island" = c(-156.3651,18.8049,-154.765,20.4064)
)

gtmap <- function(loc) {
  get_map(location = loc,
          maptype = "terrain-background",
          source = "stamen",
          scale = "auto")
}

mget_map <- memoise(gtmap)

####################################################
#           Make overview plot  dataframe          # 
####################################################
# setup overview plot groups
plot_spp_ids <- data1 %>%
  dplyr::mutate(collection_type = ifelse(worms_on_sample %in% c("No", "?"), "No Nematode",
                                   ifelse(worms_on_sample == "Tracks", "Tracks only",
                                          ifelse(worms_on_sample == "Yes" & is.na(pcr_positive), "Not genotyped",
                                                 ifelse(worms_on_sample == "Yes" & pcr_positive == 0, "PCR -",
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
                                                                                "Oscheius sp.",
                                                                                "Panagrolaimus sp.",
                                                                                NA),
                                                              "Other PCR +", species_id)))))) %>%
  dplyr::select(c_label, collection_type, island, species_id, pcr_positive, worms_on_sample, longitude, latitude) %>%
  dplyr::distinct(c_label, collection_type, .keep_all=T) %>% 
  dplyr::group_by(c_label) %>%
  dplyr::mutate(multiple_type = ifelse(n() > 1, "yes", "no")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(collection_type = forcats::as_factor(collection_type),
                collection_type = forcats::fct_relevel(collection_type,
                                                       "C. elegans",
                                                       "C. oiwi",
                                                       "C. tropicalis",
                                                       "C. briggsae",
                                                       "Other PCR +",
                                                       "PCR -",
                                                       "Not genotyped",
                                                       "Tracks only",
                                                       "No Nematode")) %>%
  dplyr::arrange(collection_type) %>% # arrange sets order for c-labels with multiples so highest priority collection type is on top
  dplyr::distinct(c_label, .keep_all = T) %>% # selects highest priority collection type from a c-label with multiple collection types on it
  dplyr::arrange(desc(collection_type)) # reorders c-labels so highest priority collections are plotted on top

####################################################
#  Bar chart inset                   # 
####################################################
# Plot stacked bar chart
bar_chart <- plot_spp_ids %>%
  dplyr::group_by(collection_type, island) %>%
  dplyr::mutate(collections_per_island = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(island) %>%
  dplyr::mutate(total_collections = n(), perc_class_island = collections_per_island / total_collections * 100) %>%
  dplyr::arrange(total_collections) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(collection_type, island, .keep_all = T) %>%
  dplyr::mutate(island = factor(island, levels = names(island_palette))) %>%
  dplyr::mutate(collection_type = factor(collection_type, levels = c("No Nematode", "Tracks only", "Not genotyped", "PCR -", "Other PCR +", "C. briggsae", "C. tropicalis", "C. oiwi", "C. elegans")))

# Fig2B plot for rhabditida positive collections
plot_bar_chart <- ggplot(data = bar_chart) +
  geom_bar(stat = "identity", aes(x = factor(island), y = perc_class_island, fill = collection_type), colour = "black") +
  scale_fill_manual(values=c(species_palette)) +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), units = "cm")) + 
  #legend.text = element_text(size = 8, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of all collections") +
  geom_text(aes(x=island, y=102, label=paste0("n=",total_collections)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 102))

plot_bar_chart_no_legend <- plot_bar_chart +
  theme(legend.position="none",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

plot_bar_chart_legend <- cowplot::plot_grid(cowplot::get_legend(plot_bar_chart))

# Figure 1 plotting function
map_overview <- function(F, label, cso, geoms, face = "plain") {
  island_set = lapply(names(islands), function(i) {
    
    l_position = "none"
    island_size = 2
    imap <- cso %>% dplyr::filter(island == i)
    rects = element_blank()
    map = mget_map(islands[[i]])
    
    # Calculate scalebar
    bb <- attr(map,"bb")
    sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                       lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                       lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                       lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))
    
    sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                      c(sbar$lon.end,sbar$lat.end))
    
    scalebar.length <- 20
    sbar$lon.end <- sbar$lon.start +
      ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000
    ptspermm <- 2.83464567
    
    base_map <- ggplot(imap) +
      ggmap::inset_ggmap(map) +
      #rects +
      geoms +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(c(0,0,0,0), "lines"),
            axis.line = element_blank(),
            plot.title = element_text(lineheight=.8, face="bold", vjust=1),
            plot.margin = unit(c(0,0,0,0), "lines"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = l_position,
            legend.background = element_rect(fill="white"),
            legend.text=element_text(size=12, color = "black", face = face)) +
      coord_equal(ratio=1) +
      scale_x_continuous(limits = islands[[i]][c(1,3)], expand = c(0, 0)) +
      scale_y_continuous(limits = islands[[i]][c(2,4)], expand = c(0, 0)) +
      geom_segment(data = sbar,
                   aes(x = lon.start,
                       xend = lon.end,
                       y = lat.start,
                       yend = lat.end),
                   arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                               ends = "both", type = "open")) +
      geom_text(data = sbar,
                aes(x = (lon.start + lon.end)/2,
                    y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
                    label = paste(format(scalebar.length),
                                  'km')),
                hjust = 0.5,
                vjust = 0,
                size = 8/ptspermm)  +
      coord_map(projection = "mercator",
                xlim=c(bb$ll.lon, bb$ur.lon),
                ylim=c(bb$ll.lat, bb$ur.lat)) +
      scale_radius(range = c(island_size, island_size), guide = "none") #+
    #scale_shape_manual(values = shape)
    
    base_map
    
  })
  
  island_set[[6]] <- F
  
  without_label <- plot_grid(plotlist = island_set,
            labels = c("A - Kauai",
                       "B - O'ahu",
                       "C - Moloka'i",
                       "D - Maui",
                       "E - Island of Hawai'i",
                       "F", ""
            ),
            label_y = 0.98,
            hjust = 0,
            label_x = 0.06,
            align = "vh")
  
  cowplot::plot_grid(without_label, label, nrow = 2, rel_heights = c(1, .05))
}



# Make map overview plot
Figure1 <- map_overview(plot_bar_chart_no_legend, plot_bar_chart_legend, plot_spp_ids,
                               c(geom_point(aes(x=longitude,
                                                y=latitude,
                                                fill=collection_type,
                                                size = 1),
                                            color="black",
                                            shape=21,
                                            stroke = 0.5
                               ),
                               scale_fill_manual("species", values = species_palette)
                               ),
                               face="italic"
)
Figure1

ggsave('plots/Figure1_REVISED.pdf', width = 12, height = 9)

# # Write data used for this plot
# plot_spp_ids %>%
# readr::write_csv("data/elife_files/fig1-data1_Revised.csv")