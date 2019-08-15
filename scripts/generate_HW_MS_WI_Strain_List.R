#Load necessary packages
library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#####################################
# generate wild isolate strain info #
#####################################

# filter full wild isolate strain list
df <- data.table::fread("data/original_WI_strain_list_remove_when_cleaning.csv")

# List of HI isotypes in manuscript
hi_manuascript <- c("CB4856", "DL238", "ECA189", "ECA191", "ECA347", "ECA363", "ECA369", "ECA372", "ECA396", "ECA701",
                    "ECA703", "ECA705", "ECA706", "ECA710", "ECA712", "ECA722", "ECA723", "ECA724", "ECA730", "ECA732",
                    "ECA733", "ECA738", "ECA740", "ECA741", "ECA742", "ECA743", "ECA744", "ECA745", "ECA746", "ECA760",
                    "ECA768", "ECA777", "ECA778", "ECA807", "ECA812", "QX1791", "QX1792", "QX1793", "QX1794", "XZ1513",
                    "XZ1514", "XZ1515", "XZ1516")

# list of isotypes in manuscript
hi_manuscript_276 <- c("AB1", "BRC20067", "BRC20263", "CB4851", "CB4852", "CB4853", "CB4854", "CB4855", "CB4856", "CB4857",
                       "CB4858", "CB4932", "CX11254", "CX11262", "CX11264", "CX11271", "CX11276", "CX11285", "CX11292",
                       "CX11307", "CX11314", "CX11315", "DL200", "DL226","DL238", "ECA189", "ECA191", "ECA252", "ECA347",
                       "ECA348", "ECA349", "ECA36", "ECA363", "ECA369", "ECA372", "ECA396", "ECA701", "ECA703", "ECA705",
                       "ECA706", "ECA710", "ECA712", "ECA722", "ECA723", "ECA724", "ECA730", "ECA732", "ECA733", "ECA738",
                       "ECA740", "ECA741", "ECA742", "ECA743", "ECA744", "ECA745", "ECA746", "ECA760", "ECA768", "ECA777",
                       "ECA778", "ECA807", "ECA812", "ED3005", "ED3011", "ED3012", "ED3017", "ED3040", "ED3046", "ED3048",
                       "ED3049", "ED3052", "ED3073", "ED3077", "EG4347","EG4349", "EG4724", "EG4725", "EG4946", "GXW1",
                       "JT11398", "JU1088", "JU1172", "JU1200", "JU1212", "JU1213", "JU1242", "JU1246", "JU1249", "JU1395",
                       "JU1400", "JU1409", "JU1440", "JU1491", "JU1530", "JU1543", "JU1568", "JU1580", "JU1581", "JU1586",
                       "JU1652", "JU1666", "JU1792", "JU1793", "JU1808", "JU1896", "JU1934", "JU2001", "JU2007", "JU2016",
                       "JU2017", "JU2106", "JU2131", "JU2141", "JU2234", "JU2250", "JU2257", "JU2316", "JU2464", "JU2466",
                       "JU2478", "JU2513", "JU2519", "JU2522", "JU2526","JU2534", "JU2565", "JU2566", "JU2570", "JU2572",
                       "JU2575", "JU2576", "JU2578", "JU258", "JU2581", "JU2586", "JU2587", "JU2592", "JU2593", "JU2600", 
                       "JU2610", "JU2619", "JU2800", "JU2811", "JU2825", "JU2829", "JU2838", "JU2841", "JU2853", "JU2862",
                       "JU2866", "JU2878", "JU2879", "JU2906", "JU2907", "JU310", "JU311", "JU3125", "JU3127", "JU3128",
                       "JU3132", "JU3134", "JU3135", "JU3137", "JU3140", "JU3144", "JU323", "JU346", "JU360", "JU367",
                       "JU393", "JU394", "JU397", "JU406", "JU440", "JU561","JU642", "JU751", "JU774", "JU775", "JU778",
                       "JU782", "JU792", "JU830", "JU847", "KR314", "LKC34", "LSJ1", "MY1", "MY10", "MY16", "MY18", "MY2147",
                       "MY2212", "MY23", "MY2453", "MY2530", "MY2535", "MY2573", "MY2585", "MY2693", "MY2713", "MY2741",
                       "MY518","MY679", "MY772", "MY795", "MY920", "N2", "NIC1", "NIC1049", "NIC1107", "NIC166", "NIC195",
                       "NIC199", "NIC2", "NIC207", "NIC231", "NIC236", "NIC242", "NIC251", "NIC252", "NIC255", "NIC256",
                       "NIC258", "NIC259", "NIC260", "NIC261", "NIC262", "NIC265", "NIC266", "NIC267", "NIC268", "NIC269",
                       "NIC271", "NIC272", "NIC274", "NIC275", "NIC276", "NIC277", "NIC3", "NIC501", "NIC511", "NIC513",
                       "NIC514", "NIC515", "NIC522", "NIC523", "NIC526", "NIC527", "NIC528", "NIC529", "PB303", "PB306",
                       "PS2025", "PX179","QG2075", "QG536", "QG556", "QG557", "QW947", "QX1211", "QX1212", "QX1233",
                       "QX1791", "QX1792", "QX1793", "QX1794", "RC301", "WN2001", "WN2002", "WN2033", "WN2050", "XZ1513",
                       "XZ1514", "XZ1515", "XZ1516")

# filter WI sheet to 276 isotypes in manuscript
hi_276_isotypes <- df %>%
  dplyr::filter(isotype %in% hi_manuscript_276)

# remove strains from 20180413 release
cendr_249_isotypes <- hi_276_isotypes %>%
  dplyr::filter(release != "20180413" & sequenced == 1)

# select only hawaiian strains from 20180413 release
cendr_20180413_hi_isotypes <- hi_276_isotypes %>%
  dplyr::filter(release == "20180413" & state == "Hawaii")

# join 249 set with hawaiian strains from 20180413 release and write new wild isolate strain list
HW_MS_wild_isolate_list <- full_join(cendr_249_isotypes, cendr_20180413_hi_isotypes) %>%
  #dplyr::filter(sequenced == 1)%>%
  dplyr::select(-sequenced, -warning_message, -set_divergent, -set_1, -set_2, -set_3, -set_4, -set_5) %>%
  readr::write_csv("data/WI_strain_list.csv")

# check that isotypes are correct number
num_isotypes <- length(unique(HW_MS_wild_isolate_list$isotype))
