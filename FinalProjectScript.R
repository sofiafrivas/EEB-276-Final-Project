require(librarian)
librarian::shelf(tidyverse, here, janitor, googlesheets4, lubridate, splitstackshape,
                 dplyr,ggplot2,pwr2,tidyr, broom, ggpubr, paletteer)

#don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

gonadindex_raw <- readRDS("quad_build3.rds")


gi_predictors <- gonadindex_raw %>%
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  select(-c(patch_id, latitude, longitude, survey_date, site, site_type,
            pred_patch,patch_cat, zone, year, sd_biomass, sd_gi, se_gonad_mass_g, geometry))


model <- glm(mean_gi ~ ., data = gi_predictors, family = gaussian())
summary(model)

glm(mean_gi ~ macr, data = gi_predictors, family = gaussian())
