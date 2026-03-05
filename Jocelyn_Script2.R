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

# Standardize Continuous Variables ----------------------------------------

predictors_simple2 <- predictors_simple %>%
  mutate(across(c(n_macro_plants_m2: cov_dictyoneurum_spp),
                ~ scale(.)[,1],
                .names = "{.col}_z"))

# Setting Priors ----------------------------------------------------------

priors <- c(
  prior(normal(0, 1), class = "b"),        # slopes
  prior(normal(6, 3), class = "Intercept"), # based on GI range
  prior(exponential(1), class = "sd")      # random effect SD
)

# BRMS Exploration --------------------------------------------------------

colnames(predictors_simple2)
library(brms)

# model with site_type
model_01 <- brm( mean_gi ~ 
    kelp_cover_z +
    red_urchin_density_z +
    purple_urchin_density_z +
    zone +
    year +
    (1 | site),
  data = predictors_simple2,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  seed = 123,
  file = "output/model_01"
)

summary(mean_gi_model)
plot(mean_gi_model)