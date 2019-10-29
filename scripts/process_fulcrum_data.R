#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(lubridate)
library(geosphere)
library(googlesheets)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

LL <- c("longitude", "latitude")
LLA <- c("longitude", "latitude", "altitude")

# Team
RAPTORS <- c("dec@u.northwestern.edu",
             "daehan.lee@northwestern.edu",
             "erik.andersen@northwestern.edu",
             "stefanzdraljevic2018@u.northwestern.edu")


filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

FtoC <- function(F) {
  (F - 32)*(5/9)
}

# Read in C-labels
sc <- readr::read_csv("data/fulcrum/sample_collection.csv") %>%
  dplyr::mutate(c_label = stringr::str_to_upper(c_label)) %>%
  dplyr::filter(project != "Worm Meeting Collection") %>%
  dplyr::rename(sampled_by = created_by) %>%
  dplyr::select(-updated_at,
                -system_created_at,
                -system_updated_at,
                -date) %>%
  dplyr::mutate(datetime = lubridate::ymd_hms(created_at, tz = "HST")) %>%
  dplyr::mutate(date = lubridate::date(created_at)) %>%
  dplyr::select(-created_at) %>%
  # Label substrate moisture issue (moisture meters were on potentially wrong settings at times before 2017-08-09)
  dplyr::mutate(substrate_moisture = ifelse(substrate_moisture == -1, NA, substrate_moisture)) %>%
  dplyr::mutate(substrate_moisture_issue = !(lubridate::ymd(date) %within% lubridate::interval("2017-08-09", "2017-08-31"))) %>%
  # Two observations had a C > 50; Clearly wrong.
  dplyr::mutate(substrate_temperature = ifelse(substrate_temperature == 100, NA, substrate_temperature)) %>%
  # Fix Fahrenheit observations
  dplyr::mutate(substrate_temperature = ifelse(substrate_temperature > 35,
                                               FtoC(substrate_temperature),
                                               substrate_temperature)) %>%
  # Fix substrate_temp_c when C < 9
  dplyr::mutate(substrate_temperature = ifelse(substrate_temperature < 9,
                                               NA,
                                               substrate_temperature)) %>%
  # Fix ambient temp F to C
  dplyr::mutate(ambient_temperature = ifelse(ambient_temperature > 50,
                                             FtoC(ambient_temperature),
                                             ambient_temperature)) %>%
  # Fix mispelling
  dplyr::mutate(substrate = ifelse(substrate == "Millipeed",
                                   "Millipede",
                                   substrate))



# Read in S-labels
po <- readr::read_csv("data/fulcrum/plating_out.csv") %>%
  dplyr::select(c_label_id = c_label,
                po_id = fulcrum_id,
                po_created_at = system_created_at,
                po_created_by = created_by,
                worms_on_sample,
                approximate_number_of_worms,
                males_observed,
                dauers_on_sample,
                approximate_number_of_worms,
                po_date = date,
                po_time = time,
                po_latitude = latitude,
                po_longitude = longitude)


# Add data from collection photos lat, long, elevation. Only need to run once to extract data.
# Read in data from photos. Need to install using ‘brew install exiftool’ in terminal.
# comm <- paste0("exiftool -coordFormat '%+.6f' -csv -ext jpg ",
#                getwd(),
#                "/data/fulcrum/photos/id/*")

# Exif Data
# exif <- readr::read_csv(pipe(comm)) %>%
#   dplyr::mutate(SourceFile = stringr::str_replace(basename(SourceFile), ".jpg", "")) %>%
#   dplyr::select(sample_photo = SourceFile,
#                 altitude = GPSAltitude,
#                 latitude = GPSLatitude,
#                 longitude = GPSLongitude,
#                 ExposureTime,
#                 Artist,
#                 Aperture,
#                 BrightnessValue,
#                 PhotoDate = DateCreated,
#                 FOV) %>%
#   dplyr::mutate(altitude =  as.numeric(stringr::str_replace(altitude, " m", ""))) %>%
#   dplyr::mutate(FOV =  as.numeric(stringr::str_replace(FOV, " deg", ""))) %>%
#   dplyr::group_by(sample_photo) %>%
#   # Only retain data from one sample photo.
#   dplyr::distinct(.keep_all=T)
# save(file = "data/fulcrum/exif.Rda", exif)
# load data from images already processed by Exif
load("data/fulcrum/exif.Rda")

# Join Data
df <- dplyr::full_join(po, sc, by = c("c_label_id" = "fulcrum_id")) %>%
  dplyr::rename(record_latitude = latitude, record_longitude = longitude) %>%
  dplyr::select(c_label,
                everything(),
                -c_label_id,
                -sample_photo_url) %>%
  dplyr::left_join(exif) %>%
  # In rare cases, lat/lon not with photo; fix.
  dplyr::mutate(latitude = ifelse(is.na(latitude), record_latitude, latitude)) %>%
  dplyr::mutate(longitude = ifelse(is.na(longitude), record_longitude, longitude)) %>%
  dplyr::mutate(ambient_temperature = as.numeric(ambient_temperature)) %>%
  dplyr::mutate(ambient_temperature = ifelse(ambient_temperature > 70,
                                             ((5/9)*(ambient_temperature-32)),
                                             ambient_temperature)) %>%
  dplyr::mutate_at(.vars = vars(dplyr::starts_with("gps")),
                   .funs = funs(as.numeric)) %>%
  dplyr::mutate(team = ifelse(sampled_by %in% RAPTORS, "RAPTORS", "MOANA")) %>%
  dplyr::mutate(worms_on_sample = ifelse(is.na(worms_on_sample), "?", worms_on_sample)) %>%
  dplyr::filter(!is.na(c_label)) %>%
  dplyr::select(-assigned_to,
                -status,
                -Artist) %>%
  # Calculate the Haversine distance
  dplyr::rowwise() %>%
  dplyr::mutate(gps_err = geosphere::distHaversine(c(longitude, latitude),
                                                   c(record_longitude, record_latitude))) %>%
  dplyr::ungroup()

# Generate dataset mapping C-labels to S-labels
po_slabels <- readr::read_csv("data/fulcrum/plating_out_s_labeled_plates.csv") %>%
  dplyr::select(fulcrum_parent_id, s_label) %>%
  dplyr::left_join(df, by = c("fulcrum_parent_id" = "po_id")) %>%
  dplyr::select(c_label,
                s_label,
                worms_on_sample,
                males_observed,
                dauers_on_sample,
                approximate_number_of_worms,
                po_date,
                po_time,
                longitude,
                latitude,
                substrate_temperature,
                substrate_moisture,
                ambient_humidity,
                ambient_temperature,
                sampled_by)

po_slabels <- po_slabels %>% 
  dplyr::filter(!is.na(c_label), !is.na(s_label))

# Add Nested S-labels
df <- dplyr::left_join(df, 
                       po_slabels %>%
                         dplyr::select(c_label, s_label) %>%
                         dplyr::group_by(c_label) %>%
                         dplyr::summarize(s_label_cnt = length(s_label), 
                                          s_label = paste0(s_label, collapse = ","))) %>%
  dplyr::select(c_label, s_label, s_label_cnt, everything(), -po_id)


# Samples collected that were never processed.
c_labels_never_processed = df %>% 
  dplyr::filter(worms_on_sample == "?") %>%
  dplyr::select(c_label)


# Keep track of points that are corrected
df$GPS_corrected = F

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Part 1 Average GPS location of off or missing points    #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

df <-dplyr::arrange(df, team, date, time)
wrong_pos <- c("C-0804",
               "C-0240",
               "C-0742",
               "C-0549",
               "C-1066",
               "C-0642",
               "C-1203",
               "C-1132",
               "C-1308",
               "C-1319",
               "C-2598",
               "C-2599")

sapply(wrong_pos, function(c_label) {
  row_num <- which(df$c_label == c_label)
  p1 <- unlist(df[row_num-1, LL])
  p2 <- unlist(df[row_num+1, LL])
  # The <<- goes to the global environment.
  df[row_num, LL] <<- gcIntermediate(p1, p2, n = 1, addStartEnd=FALSE)
  df[row_num, "altitude"] <<- mean(unlist(df[c(row_num -1, row_num + 1), "altitude"]))
  df[row_num, c("GPS_corrected")] <<- T
})

# * * * * * * * * * * * * * * * * * * * * #
# Part 2: Use Record GPS to correct point #
# * * * * * * * * * * * * * * * * * * * * #
wrong_pos <- c("C-1079",
               "C-1072",
               "C-0561",
               "C-0588",
               "C-0772",
               "C-0642",
               "C-1330",
               "C-1231",
               "C-0869",
               "C-0339",
               "C-1365",
               "C-2527",
               "C-2543",
               "C-2547",
               "C-1462",
               "C-1460",
               "C-1466",
               "C-2591")

sapply(wrong_pos, function(c_label) {
  row_num <- which(df$c_label == c_label)
  # The <<- goes to the global environment.
  df[row_num, LL] <<- df[row_num, c("record_longitude", "record_latitude")]
  df[row_num, c("GPS_corrected")] <<- T
})

#=============================#
# Part 3: Manual Corrections #
#===========================#

# Part 3: Manual Corrections

positions <- df %>% dplyr::filter(c_label %in% c("C-0447", "C-0446")) %>% dplyr::select(longitude, latitude)
p1 <- positions[1,]
p2 <- positions[2,]
corrected_lat_lon <- gcIntermediate(p1, p2, n=1, addStartEnd=F)
df[df$c_label == "C-1077", LL] <- corrected_lat_lon
df[df$c_label == "C-1078", LL] <- corrected_lat_lon
df[df$c_label %in% c("C-1078","C-1077"), "altitude"] <- df[df$c_label %in% c("C-0447") , "altitude"]
df[df$c_label %in% c("C-1078","C-1077"), "GPS_corrected"] <- T

#====================================#
# Part 4: Fix Day 3 sample from Maui #
#====================================#

c_label_update <- (df[df$latitude == unlist(df[df$c_label == "C-1222","latitude"]), "c_label"])$c_label
df[df$c_label %in% c_label_update, "longitude"] <- -156.150537
df[df$c_label %in% c_label_update, "latitude"] <- 20.854776
df[df$latitude == unlist(df[df$c_label == "C-1222","latitude"]), "GPS_corrected"] = T

#====================================#
# Part 5: Fix Day 3 sample from Maui #
#====================================#

#Change C-1226, 1230, and 1234 to average between samples C-1238 and C-1222
positions <- df[df$c_label %in% c("C-1238", "C-1222"), LL]
p1 <- positions[1,]
p2 <- positions[2,]
corrected_lat_lon <- gcIntermediate(p1, p2, n=1, addStartEnd=F)
df[df$c_label == "C-1226", LL] <- corrected_lat_lon
df[df$c_label == "C-1230", LL] <- corrected_lat_lon
df[df$c_label == "C-1234", LL] <- corrected_lat_lon
df[df$c_label %in% c("C-1226", "C-1230", "C-1234"), "GPS_corrected"] = T

#====================================#
# Part 6: Fix Day 3 sample from Maui #
#====================================#

#Check maui day 4
#13 samples on Big Island, should be at Waihou Spring Forrest Reserve. Change to estimated gps
errant_latitude <- (df[which(df$c_label == "C-1487"), "latitude"])$latitude
errant_c_labels <- (df[df$latitude == errant_latitude,"c_label"])$c_label
corrected_lat_lon <- df[df$c_label == "C-2627", LL]
sapply(errant_c_labels, function(c_label) {
  row_num <- which(df$c_label == c_label)
  # The <<- goes to the global environment.
  df[row_num, LL] <<- corrected_lat_lon
  df[row_num, c("GPS_corrected")] <<- T
})

#====================================#
# Part 7: Fix Day 5 sample from Maui #
#====================================#

#Check Maui day 5
df[df$c_label == "C-2521", LL] <- df[df$c_label == "C-2554", LL]

#=============================#
# Part 8: Gridsect Positions #
#===========================#

# Correct a group by another position from the same gridsect
wrong_pos <- c("C-0708", "C-1064", "C-0617", "C-1063", "C-0616")
df[df$c_label %in% wrong_pos, LLA] <- df[df$c_label == "C-1061", LLA]
df[df$c_label %in% wrong_pos, "GPS_corrected"] <- T

sapply(wrong_pos, function(c_label) {
  row_num <- which(df$c_label == c_label)
  p1 <- unlist(df[row_num-1, LL])
  p2 <- unlist(df[row_num+1, LL])
  # The <<- goes to the global environment.
  df[row_num, LL] <<- gcIntermediate(p1, p2, n = 1, addStartEnd=FALSE)
  df[row_num, c("GPS_corrected")] <<- T
})


#=============================#
# Part 9: Fix gridsect values #
#=============================#

gridsect_modified = c("C-0567",
                      "C-0266",
                      "C-0282",
                      "C-0271",
                      "C-0777",
                      "C-1133",
                      "C-2605",
                      "C-2782")

# 20170807 - C-0567 : C2 -> D2
df[df$c_label == "C-0567", "gridsect_direction"] = "D"
# 20170809 - C-0266,0282,0271 : B123 -> D123
df[df$c_label %in% c("C-0266", "C-0282", "C-0271"), "gridsect_direction"] = "D"
# 20170812 - C-0777 : E2 -> F2
df[df$c_label == "C-0777", "gridsect_direction"] = "F"
# 20170813 - C-1133 : B1 -> C1
df[df$c_label == "C-1133", "gridsect_direction"] = "C"
# 20170817 - C-2605 : E3 -> F3
df[df$c_label == "C-2605", "gridsect_direction"] = "F"
#20170817, multiple B3s
df[df$c_label == "C-2782", "gridsect_radius"] = "1"

df$gridsect_corrected = F
df[df$c_label %in% gridsect_modified, "gridsect_corrected"] <- T

# Convert radius to numeric
df <- df %>% dplyr::rowwise() %>%
  dplyr::mutate(gridsect_radius = as.integer(
    substr(gridsect_radius,1,1)
  ))


#=======================#
# Reclassify substrates #
#=======================#

# * Merge fruit/nut/vegetable
# * Reclassify rotting in substrate other as rotting wood
# * Create vegetation category.

substrate_merge <- c("Fruit",
                     "Rotting fruit",
                     "Nut",
                     "Rotting nut")

df <- df %>% dplyr::ungroup() %>%
  dplyr::mutate(substrate = ifelse(substrate %in% substrate_merge, "Fruit/nut/vegetable", substrate)) %>%
  dplyr::mutate(substrate = ifelse(substrate %in% c("Rotting fungus"), "Fungus", substrate)) %>%
  dplyr::mutate(substrate = ifelse(grepl("rotting|Rotting", substrate_other), "Rotting wood", substrate)) %>%
  dplyr::mutate(substrate = ifelse(is.na(substrate), "Vegetation", substrate)) %>%
  dplyr::mutate(substrate = ifelse(substrate == "Soil", "Leaf litter", substrate)) %>%
  dplyr::mutate(substrate = ifelse(substrate == "Rotting vegetable", "Vegetation", substrate))

#===============================#
#  Add flags for runs of values #
#===============================#

df <- dplyr::arrange(df, team, date, time) %>%
  dplyr::group_by(team) %>%
  dplyr::mutate(ambient_run_flag = (ambient_humidity == dplyr::lag(ambient_humidity)) &
                  (ambient_temperature == dplyr::lag(ambient_temperature))
                & (gridsect == "no")) 

#==================#
# Update altitudes #
#==================#
# only need to run once to get altitudes.
# library(geonames)
# options(geonamesUsername="katiesevans")
# altitudes <- df %>% dplyr::ungroup() %>%
#        dplyr::select(c_label, latitude, longitude, altitude) %>%
#        dplyr::rowwise() %>%
#        dplyr::mutate(altitude = ifelse(is.na(altitude),
#                                       geonames::GNsrtm3(latitude, longitude)$srtm3,
#                                       altitude)
#                      ) %>%
#         dplyr::ungroup()
# save(altitudes, file = "data/fulcrum/altitude.Rda")

load("data/fulcrum/altitude.Rda")

df <- df %>% dplyr::ungroup() %>%
  dplyr::select(-altitude) %>%
  dplyr::left_join(altitudes, by = c("c_label", "longitude", "latitude"))

#=========================================#
# Remove dauer estimates in Evanston (NA) #
#=========================================#

evanston_collectors <- c("joost.vanderzwaag@wur.nl",
                         "steffen.hahnel@northwestern.edu",
                         "robyn.tanny@northwestern.edu",
                         "tcrombie@northwestern.edu",
                          NA)

df <- df %>% 
  dplyr::mutate(po_in_evanston = ifelse(po_created_by %in% evanston_collectors, T, F)) %>%
  dplyr::mutate(dauers_on_sample = ifelse(po_in_evanston, NA, dauers_on_sample))

#=============================#
# Part X: Set the islands!    #
#=============================#

# Create Island Column
df$island <- "?"
df[filter_box(df$longitude, df$latitude, c(-158.3617,21.1968,-157.5117,21.7931)), "island"] <- "Oahu"
df[filter_box(df$longitude, df$latitude, c(-159.9362, 21.6523, -159.1782, 22.472)), "island"] <- "Kauai"
df[filter_box(df$longitude, df$latitude, c(-157.327, 21.0328, -156.685, 21.2574)), "island"] <- "Molokai"
df[filter_box(df$longitude, df$latitude, c(-156.7061, 20.4712, -155.9289, 21.0743)), "island"] <- "Maui"
df[filter_box(df$longitude, df$latitude, c(-156.1346, 18.6619, -154.6985, 20.4492)), "island"] <- "Big Island"

# Fix errant GPS locations from team Moana and Erik
df[df$island == "BIG_ISLAND" & df$team == "MOANA" , c("latitude", "longitude")] <- NA
df[df$island == "BIG_ISLAND" & df$team == "MOANA" & !is.na(df$island), "island"] <- "MAUI"

# Create Trail Column
df$location <- NA

df[filter_box(df$longitude, df$latitude, c(-157.72537,21.303309,-157.71919,21.32122)), "location"] <- "Kuliouou Ridge Trail"
df[filter_box(df$longitude, df$latitude, c(-158.0192352613,21.5014265529,-158.0145925283,21.5041245046)), "location"] <- "Wahiawa Botanical Garden"
df[filter_box(df$longitude, df$latitude, c(-157.8598800302,21.3149311581,-157.855797708,21.3182194587)), "location"] <- "Foster Community Garden"
df[filter_box(df$longitude, df$latitude, c(-157.7829487403,21.3569863645,-157.7752268314,21.3655295525)), "location"] <- "Maunawili Demonstration Trail"
df[filter_box(df$longitude, df$latitude, c(-157.8014534712,21.3322593,-157.798127532,21.3427719396)), "location"] <- "Manoa Falls Trail"
df[filter_box(df$longitude, df$latitude, c(-157.8135502338,21.3779082884,-157.7915561199,21.3970691079)), "location"] <- "Ho'omaluhia Botanical Garden"
df[filter_box(df$longitude, df$latitude, c(-159.613624,22.167098,-159.575601,22.226422)), "location"] <- "Na Pali Coast State Wilderness Park"

#==========================#
# Setup gridsect variables #
#==========================#

grids <- df %>% dplyr::filter(gridsect == 'yes') %>%
  select(c_label, gridsect_direction, gridsect_radius, island, datetime, team, date) %>%
  dplyr::group_by(c_label) %>%
  dplyr::distinct()

missing_gridsect <- tibble::tibble(c_label = NA,
                                   gridsect_direction = "E",
                                   gridsect_radius = 3,
                                   island = "Maui",
                                   datetime = as.POSIXct("2017-08-16 15:56:00 HST"),
                                   team = "MOANA",
                                   date = as.Date("2017-08-16"))

grids <- dplyr::bind_rows(grids, missing_gridsect) %>%
  dplyr::arrange(island, date, team, datetime)

grids$grid_num <- sort(rep(1:20, 19))

df <- df %>% dplyr::left_join(grids %>% dplyr::select(c_label, grid_num), 
                              by = c("c_label"))


#===============#
# Add photo URL #
#===============#
df <-df %>% dplyr::rowwise() %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(photo = paste0(c_label,
                               ".",
                               stringr::str_to_lower(stringr::str_replace_all(substrate, "[^[:alnum:]]", "_")),
                               ".1.jpg"),
                photo_url = paste0("https://storage.googleapis.com/elegansvariation.org/photos/hawaii2017/",
                                   c_label,
                                   ".jpg"),
                photo_url_thumb =  paste0("https://storage.googleapis.com/elegansvariation.org/photos/hawaii2017/",
                                          c_label,
                                          ".thumb.jpg")) %>%
  dplyr::ungroup()
# 
# photo_comms <- df %>% dplyr::mutate(sample_photo = str_split(sample_photo, ",")) %>%
#   dplyr::select(-s_label) %>%
#   dplyr::select(c_label, sample_photo, substrate) %>%
#   tidyr::unnest() %>%
#   dplyr::group_by(c_label) %>%
#   dplyr::mutate(comm = paste0("cp ../data/photos/id/",
#                               sample_photo,
#                               ".jpg",
#                               " ",
#                               "../data/photos/c/",
#                               c_label,
#                               ".",
#                               stringr::str_to_lower(str_replace_all(substrate, "[^[:alnum:]]", "_")),
#                               ".",
#                               dplyr::row_number(c_label),
#                               ".jpg")) %>%
# 
#   dplyr::select(-c_label, comm)
# 
# writeLines(photo_comms$comm, con = file("scripts/rename_photos.sh"))

# redefine
cso <- po_slabels

# Fold in variables from df to the cso data frame
cso <- cso %>% dplyr::left_join(
  df %>% dplyr::select(c_label,
                       substrate,
                       landscape,
                       sky_view,
                       photo_url,
                       photo_url_thumb,
                       altitude,
                       team,
                       island,
                       location,
                       date,
                       time,
                       FOV,
                       po_created_by,
                       dplyr::starts_with("grid")),
  by = "c_label"
)

# Merge in blast data; Take top hit
blast_results <- readr::read_tsv("data/sanger/blast_results.tsv") %>%
  dplyr::group_by(s_plate) %>%
  dplyr::filter(row_number() == 1)

#==============================#
# Load manual curation results #
#==============================#
cso <- cso %>% dplyr::left_join(blast_results, by = c("s_label" = "s_plate")) %>%
  dplyr::left_join(
    googlesheets::gs_key("1bavR10CEyvWt2zBSNBz-ADXx06b1mDFmuvaaM8Uobi4") %>%
      googlesheets::gs_read("Full", na = c("#N/A", "NA", ""),
    by = c("c_label", "s_label"))
  ) %>%
  dplyr::mutate_at(vars("pcr_rhpositive"), funs(as.numeric))

# Fix C-3295
df[df$c_label == "C-3295", c("longitude", "latitude", "record_longitude", "record_latitude", "date", "datetime")] <- NA

# Merge in missing c-labels that don't exist in cso but do exist in df
cso <- cso %>% dplyr::bind_rows(df[colnames(df) %in% colnames(cso)] %>%
                                           dplyr::filter(!(c_label %in% cso$c_label)))  %>%
  dplyr::arrange(c_label, s_label) %>%
  dplyr::mutate(spp_id = ifelse(
    (pcr_rhpositive == 0) | 
      ((pcr_rhpositive == 1) & is.na(spp_id)),
    "Unknown",
    spp_id)
  )

#=================#
# Data Correction #
#=================#

# remove C-3295 because we don't know what happened with that collection
remove_cso <- cso$c_label == "C-3295"
cso <- cso[!remove_cso,]
remove_df <- df$c_label == "C-3295"
df <- df[!remove_df,]

# correct substrate moisture values less than or equal to zero
cso <- cso%>%
  dplyr::mutate(substrate_moisture = ifelse(substrate_moisture <= 0, NA, substrate_moisture))
df <- df%>%
  dplyr::mutate(substrate_moisture = ifelse(substrate_moisture <= 0, NA, substrate_moisture))

#====================#
# Integrate isotypes #
#====================#

wi_info_sheet <- gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
                 gs_read() %>%
                 dplyr::select(strain, isotype, c_label, s_label)

cso <- cso %>%
          dplyr::left_join(wi_info_sheet,
                           by = c("c_label", "s_label")) %>%
          dplyr::select(isotype, strain, dplyr::everything())

df <- df %>% dplyr::left_join(wi_info_sheet %>%
                        dplyr::select(-s_label),
                        by = c('c_label')) %>%
       dplyr::select(isotype, strain, dplyr::everything())

save(file = "data/fulcrum/df.Rda", df, cso)

#==============================#
# Generate Supplemental Data 1 #
#==============================#
cso %>%
  dplyr::select(isotype,
                strain,
                c_label,
                s_label,
                worms_on_sample,
                approximate_number_of_worms,
                latitude,
                longitude,
                island,
                substrate,
                landscape,
                sky_view,
                gridsect,
                gridsect_number = grid_num,
                gridsect_direction,
                gridsect_radius,
                substrate_temperature,
                substrate_moisture,
                ambient_humidity, 
                ambient_temperature,
                altitude,
                date,
                time,
                pcr_positive = pcr_rhpositive,
                species_id = spp_id,
                photo_url_thumb) %>%
  dplyr::mutate(species_id = ifelse(species_id == "C. sp. 53", "C. oiwi", species_id)) %>%
  dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit",
                                         ifelse(substrate == "Rotting flower", "Flower",
                                                ifelse(substrate == "Rotting fungus", "Fungus",
                                                       ifelse(substrate %in% c("Rotting wood",
                                                                               "Compost",
                                                                               "Soil",
                                                                               "Grass"), 
                                                              "Vegetation",
                                                              ifelse(substrate %in% c("Isopod", "Millipede", "Slug"), "Invertebrate", substrate)))))) %>%
  readr::write_tsv("data/Supplemental Data 1.tsv")

# # Write elife supplemental data 1
# cso %>%
#   dplyr::select(isotype,
#                 strain,
#                 c_label,
#                 s_label,
#                 worms_on_sample,
#                 approximate_number_of_worms,
#                 latitude,
#                 longitude,
#                 island,
#                 substrate,
#                 landscape,
#                 sky_view,
#                 gridsect,
#                 gridsect_number = grid_num,
#                 gridsect_direction,
#                 gridsect_radius,
#                 substrate_temperature,
#                 substrate_moisture,
#                 ambient_humidity, 
#                 ambient_temperature,
#                 altitude,
#                 date,
#                 time,
#                 pcr_positive = pcr_rhpositive,
#                 species_id = spp_id,
#                 photo_url_thumb) %>%
#   dplyr::mutate(species_id = ifelse(species_id == "C. sp. 53", "C. oiwi", species_id)) %>%
#   dplyr::mutate(fixed_substrate = ifelse(substrate == "Fruit/nut/vegetable", "Fruit",
#                                          ifelse(substrate == "Rotting flower", "Flower",
#                                                 ifelse(substrate == "Rotting fungus", "Fungus",
#                                                        ifelse(substrate %in% c("Rotting wood",
#                                                                                "Compost",
#                                                                                "Soil",
#                                                                                "Grass"), 
#                                                               "Vegetation",
#                                                               ifelse(substrate %in% c("Isopod", "Millipede", "Slug"), "Invertebrate", substrate)))))) %>%
#   dplyr::select(-substrate) %>%
#   readr::write_csv("data/elife_files/Supplemental_Data_1.csv")
