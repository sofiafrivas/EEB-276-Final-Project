# Notes -------------------------------------------------------------------
#stepAIC function to compare models 
#compare AIC scores and r squared (pseudo r squared) 
#subsample data and compare to whole data? 
#total mean of all recruits (add columns for total mean) 
#look into different family for glm models

# Load First --------------------------------------------------------------

require(librarian)
librarian::shelf(tidyverse, here, janitor, googlesheets4, lubridate, splitstackshape,
                 dplyr,ggplot2,pwr2,tidyr, broom, ggpubr, paletteer, performance)

library(MASS) #stepAIC
library(ggeffects)
library(brms) #brm modeling 
library(cmdstanr)

#don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

# Read in Dataframe(s) ----------------------------------------------------

gonadindex_raw <- readRDS("quad_build3.rds")

# Data Wrangling ----------------------------------------------------------

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

# Modeling ----------------------------------------------------------------

#PCA 
priors_horseshoe <- c(
  prior(horseshoe(df = 3), class = "b"),
  prior(normal(6.5, 2), class = "Intercept"),
  prior(student_t(3, 0, 2), class = "sigma"))

model_horseshoe <- brm(mean_gi ~ .,
                data = predictors_simple,
                family = gaussian(),
                prior = priors_horseshoe,
                chains = 4, iter = 6000, warmup = 3000,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                cores = 4)
summary(model_horseshoe)

loo_horseshoe  <- loo(model_horseshoe)

model_pca <- brm(mean_gi ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8,
                 data = pc_data,
                 family = gaussian(),
                 prior = priors,
                 chains = 4, iter = 4000, warmup = 2000,
                 cores = 4)

loo_pca <- loo(model_brms_final)

loo_compare(loo_hs, loo_pca)



#GLM 
model_glm <- glm(mean_gi ~ ., data = predictors_simple, family = gaussian())

summary(model_glm)
stepAIC(model_glm)

