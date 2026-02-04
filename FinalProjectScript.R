require(librarian)
librarian::shelf(tidyverse,here, janitor, googlesheets4, lubridate, splitstackshape,
                 googledrive,httpuv,dplyr,ggplot2,pwr2,tidyr,broom,ggpubr, paletteer)

load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

readRDS("quad_build3.rds")
