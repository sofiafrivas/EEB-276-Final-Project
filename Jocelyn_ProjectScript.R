# Jocelyn's script


# Loading Relevant R Packages ---------------------------------------------

require(librarian)
librarian::shelf(tidyverse, here, janitor, googlesheets4, lubridate, splitstackshape,
                 dplyr,ggplot2,pwr2,tidyr, broom, ggpubr, paletteer, performance)

# Read in Dataframe(s) ----------------------------------------------------

# Don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

gonadindex_raw <- readRDS("quad_build3.rds")
## write.csv(gonadindex_raw, "gonadindex_raw.csv", row.names = FALSE)

# Cleaning the Dataframe --------------------------------------------------

# Renaming raw dataframe to gi_predictors
gi_predictors <- gonadindex_raw %>%
  # Grouping site details into one variable called "site_id"
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  # Removing listed varibles from the tibble
  select(-c(patch_id, latitude, longitude, survey_date, site, site_type, patch_cat, 
            zone, year, sd_biomass, sd_gi, se_gonad_mass_g, geometry)) %>% 
  # Grouping all juvenile kelp data into one variable called "juveniles"
  mutate(juveniles = macj + nerj + ptej + lsetj + eisj) %>% 
  # Changing area data to 1 m2
  mutate (n_macro_plants_m2 = n_macro_plants_20m2/20)

# New Tibble (Sofia's Data Wrangling Section)
gi_predictors <- gonadindex_raw %>%
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  dplyr::select(-c(patch_id, latitude, longitude, survey_date, site, site_type, patch_cat, 
                   zone, year, sd_biomass, sd_gi, se_gonad_mass_g, geometry)) %>% 
  mutate(juveniles = macj + nerj + ptej + lsetj + eisj,
         n_macro_plants_m2 = n_macro_plants_20m2/20,
         macro_stipe_density_m2 = macro_stipe_density_20m2/20,
         adults = (density20m2_ptecal + density20m2_eisarb + density20m2_nerlue + density20m2_lamset)/20,
         marine_snails = tegula_densitym2 + pomaulax_densitym2,
         density_cancer = density20m2_cancer_spp/20,
         density_lamstump = density20m2_lamstump/20,
         density_macstump = density20m2_macstump/20)

predictors_simple <- gi_predictors %>% 
  dplyr::select(mean_gi, n_macro_plants_m2, macro_stipe_density_m2, adults, juveniles,
                density_cancer, density_lamstump, density_macstump, relief_cm, risk_index,
                purple_urchin_densitym2, purple_urchin_conceiledm2, red_urchin_densitym2,
                red_urchin_conceiledm2, marine_snails, lamr, macr, cov_crustose_coralline,
                cov_desmarestia_spp, cov_lam_holdfast_live, cov_articulated_coralline, 
                cov_fleshy_red, cov_bare_rock, cov_green_algae, cov_barnacle, cov_bare_sand,
                cov_colonial_tunicate, cov_dictyota_dictyopteris, cov_sponge, cov_dead_kelp_holdfast_any,
                cov_mac_holdfast_live, cov_red_turf_2_cm, cov_corynactis_californica,
                cov_dictyoneurum_spp)

# Distribution Testing ----------------------------------------------------

# Plotting the distribution of "mean_gi" data
hist(gi_predictors_2$mean_gi)

# Produces a Q-Q plot to assess whether the data follows a normal distribution
# If data (more or less) follows a straight line, it is normally distributed
qqnorm(gi_predictors_2$mean_gi)

# Adds a line to Q-Q plot
qqline(gi_predictors_2$mean_gi)

# Plotting Data Across Sites ----------------------------------------------

boxplot(mean_gi ~ site, data = gonadindex_raw)
boxplot(mean_gi ~ site, data = gonadindex_raw)


# LMER Model Exploration --------------------------------------------------

model_1 <- lmer(mean_gi ~ var1 +
                         var 2 +
                         var 3 +
                         zone +
                         year +
                         (1 | site),
                       data = gonadindex_raw)


# BRMS Exploration --------------------------------------------------------

library(brms)

# model with site_type
mean_gi_model <- brm(data = gonadindex_raw,
             family = gaussian(),
             mean_gi ~ 0 + lamr + macr + site_type + (1 | site),
             iter = 2000, warmup = 1000, chains = 4, cores = 4, 
             seed = 4, 
             file = "output/mean_gi_model")

summary(mean_gi_model)
plot(mean_gi_model)

# model with pred_patch  
mean_gi_model1 <- brm(data = gonadindex_raw,
                     family = gaussian(),
                     mean_gi ~ 0 + lamr + macr + pred_patch + (1 | site),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4, 
                     seed = 4, 
                     file = "output/mean_gi_model1")

summary(mean_gi_model1)
plot(mean_gi_model1)

# Modeling ----------------------------------------------------------------

model1 <- glm(mean_gi ~ purple_urchin_densitym2 + purple_urchin_conceiledm2 + 
                cov_encrusting_red + macro_stipe_density_20m2
              , data = gi_predictors, family = gaussian())

model1 <- glm(mean_gi ~ mean_gonad_mass_g 
              + total_biomass_g +
                purple_urchin_densitym2 + 
                juveniles
              #+ purple_urchin_conceiledm2  
              #+ cov_encrusting_red 
              #+ macro_stipe_density_20m2
              , data = gi_predictors, family = gaussian)

summary(model1)
r2(model1)


model2 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                purple_urchin_densitym2 + 
                juveniles + 
                macr + 
                cov_crustose_coralline + 
                cov_mac_holdfast_live +
                lamr +
                n_macro_plants_m2 +
              , data = gi_predictors, family = gaussian)

summary(model2)
r2(model2)

model2_summary <- summary(model2)
r_squared <- 1 - (model2_summary$deviance / model2_summary$null.deviance)


glm_gamma <- glm(mean_gi ~ total_biomass_g, data =gi_predictors, family = Gamma(link = "log"))
summary(glm_gamma)

#stepAIC function to compare models 
#compare AIC scores and r squared (pseudo r squared) 
#subsample data and compare to whole data 
#total mean of all recruits (add columns for total mean) 
#look into different family for glm models

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

# Modeling Exploration ----------------------------------------------------

#command, shift, r to make a section in R

#new tibble (JR)

gi_predictors_2 <- gonadindex_raw %>%
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  select(-c(patch_id, latitude, longitude, survey_date, site, site_type, patch_cat, 
            zone, year, sd_biomass, sd_gi, se_gonad_mass_g, geometry)) %>% 
  mutate(juveniles = macj + nerj + ptej + lsetj + eisj) %>% 
  mutate (n_macro_plants_m2 = n_macro_plants_20m2/20) %>%
  mutate (densitym2_purps_on_kelp = density20m2_purps_on_kelp/20)

#model 3
model3 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                purple_urchin_densitym2 + 
                juveniles + 
                densitym2_purps_on_kelp
              , data = gi_predictors_2, family = gaussian)

summary(model3)
r2(model3)

#model 4
model4 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                purple_urchin_densitym2 + 
                red_urchin_densitym2 +
                red_urchin_conceiledm2 +
                juveniles + 
                densitym2_purps_on_kelp +
                cov_bare_sand
              , data = gi_predictors_2, family = gaussian)

summary(model4)
r2(model4)

#model 5
model5 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                purple_urchin_densitym2 + 
                risk_index +
                juveniles + 
                relief_cm +
                cov_bare_sand
              , data = gi_predictors_2, family = gaussian)

summary(model5)
r2(model5)

#model 6
model6 <- glm(mean_gi ~ n_biomass +
                total_biomass_g +
                purple_urchin_densitym2 +
                purple_urchin_conceiledm2 +
                risk_index +
                relief_cm +
                cov_bare_sand
              , data = gi_predictors_2, family = gaussian)

summary(model6)
r2(model6)

#model 7
model7 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                lamr +
                macr +
                juveniles +
                cov_mac_holdfast_live +
                cov_lam_holdfast_live +
                n_macro_plants_m2 
              , data = gi_predictors_2, family = Gamma (link = "log"))

summary(model7)
r2(model7)

#model 8
model8 <- glm(mean_gi ~ mean_gonad_mass_g +
                total_biomass_g +
                purple_urchin_densitym2 +
                purple_urchin_conceiledm2 +
                risk_index +
                relief_cm +
                cov_bare_sand
              , data = gi_predictors_2, family = Gamma (link = "log"))

summary(model8)
r2(model8)

#model 9
model 9 <- glm(mean_gi ~ mean_gonad_mass_g +
                 total_biomass_g +
                 purple_urchin_densitym2 +
                 purple_urchin_conceiledm2 +
                 risk_index +
                 relief_cm +
                 cov_bare_sand
               , data = gi_predictors_2, family = Gamma (link = "log"))


