# Notes -------------------------------------------------------------------
#stepAIC function to compare models 
#compare AIC scores and r squared (pseudo r squared) 
#subsample data and compare to whole data 
#total mean of all recruits (add columns for total mean) 
#look into different family for glm models

# Load First --------------------------------------------------------------

require(librarian)
librarian::shelf(tidyverse, here, janitor, googlesheets4, lubridate, splitstackshape,
                 dplyr,ggplot2,pwr2,tidyr, broom, ggpubr, paletteer, performance)

install.packages("MASS")   # run once if not installed
library(MASS)
library(ggeffects)


#don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

# Read in Dataframe(s) ----------------------------------------------------

gonadindex_raw <- readRDS("quad_build3.rds")
##write.csv(gonadindex_raw, "gonadindex_raw.csv", row.names = FALSE)

gi_predictors <- gonadindex_raw %>%
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  select(-c(patch_id, latitude, longitude, survey_date, site, site_type, patch_cat, 
            zone, year, sd_biomass, sd_gi, se_gonad_mass_g, geometry)) %>% 
  mutate(juveniles = macj + nerj + ptej + lsetj + eisj) %>% 
  mutate (n_macro_plants_m2 = n_macro_plants_20m2/20)
##write.csv(gonadindex_raw, "gi_predictors.csv", row.names = FALSE

# Modeling ----------------------------------------------------------------
m1 <- glm(mean_gi ~ mean_gonad_mass_g 
              + total_biomass_g +
                purple_urchin_densitym2 + 
                juveniles
              , data = gi_predictors, family = gaussian)
summary(m1) #AIC: 384.43
r2(m1) #0.599

m2 <- glm(mean_gi ~ mean_gonad_mass_g +
                        total_biomass_g +
                        purple_urchin_densitym2 + 
                        juveniles + 
                        macr + 
                        cov_crustose_coralline + 
                        cov_mac_holdfast_live +
                        lamr +
                        n_macro_plants_m2 
              , data = gi_predictors, family = gaussian)
summary(m2) #AIC: 375.79
r2(m2) #0.674
step_model <- stepAIC(m2,
                      direction = "both", 
                      trace = TRUE)

m3 <- glm(mean_gi ~ purple_urchin_densitym2 + #ecological drivers only 
                    juveniles +
                    lamr +
                    n_macro_plants_m2 +
                    macr +
                    cov_mac_holdfast_live +
                    cov_crustose_coralline,
                  family = gaussian,
                  data = gi_predictors)
summary(m3)#AIC: 421.07
r2(m3) #0.436
stepAIC(m3, direction="both")

m4 <- glm(mean_gi ~ purple_urchin_densitym2 +
            juveniles +
            lamr +
            n_macro_plants_m2 +
            macr +
            red_urchin_densitym2, 
          family = gaussian,
          data = gi_predictors)
summary(m4)#AIC: 421.07
r2(m4) #0.436
stepAIC(m4, direction="both")

# Practice Plots ----------------------------------------------------------

ggplot(gi_predictors, aes(x = pred_patch, y = mean_gi, fill = pred_patch)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = 21, fill = "white") +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c("BAR" = "#d73027", "INCIP" = "#fee090", "FOR" = "#1a9641")) +
  labs(title = "Gonad Index by Habitat Type",
       x = "Habitat Type", y = "Mean Gonad Index (GI)", fill = NULL) +
  theme_classic() +
  theme(legend.position = "none")


ggplot(gonadindex_raw, aes(x = zone, y = mean_gi)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8),
               outlier.size = 1.5) +
  labs(title = "GI by Depth Zone",
       x = "Zone", y = "Mean Gonad Index (%)") +
  theme_classic()  


