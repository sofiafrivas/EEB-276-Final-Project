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
                cov_dictyoneurum_spp) %>% 
  na.omit()


# Modeling ----------------------------------------------------------------

#dont really understand but this is the first step 
pca <- prcomp(predictors_simple %>% dplyr::select(-mean_gi), scale. = TRUE)
pc_scores <- as.data.frame(pca$x[, 1:4])
pc_data <- pc_scores %>%
  mutate(mean_gi = predictors_simple$mean_gi)

<<<<<<< HEAD
install.packages("factoextra")
library(factoextra)

=======
#setting priors 
priors_pca <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(6.5, 2), class = "Intercept"),
  prior(student_t(3, 0, 2), class = "sigma"))

#scree plot
library(factoextra)
>>>>>>> 725e90424ac751b01807966226641092588c17d1
fviz_eig(pca)

model_pca <- brm(mean_gi ~ PC1 + PC2 + PC3 + PC4,
                 data = pc_data,
                 family = gaussian(),
                 prior = priors_pca,
                 chains = 4, iter = 4000, warmup = 2000,
                 cores = 4)
summary(model_pca) #PC1 and PC3 are significant

mcmc_areas(model_pca,
           pars = vars(starts_with("b_"), -b_Intercept),
           prob = 0.95,
           prob_outer = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_classic()

pp_check(model_pca, ndraws = 100)

for(pc in c("PC1","PC2","PC3")) {
  cat("\n---", pc, "---\n")
  cat("TOP POSITIVE:\n")
  print(head(sort(pca$rotation[,pc], decreasing = TRUE), 5))
  cat("TOP NEGATIVE:\n")
  print(head(sort(pca$rotation[,pc], decreasing = FALSE), 5))
}

pc_data$fitted <- fitted(model_pca)[, "Estimate"]
pc_data$lower  <- fitted(model_pca)[, "Q2.5"]
pc_data$upper  <- fitted(model_pca)[, "Q97.5"]

ggplot(pc_data, aes(x = fitted, y = mean_gi)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted mean_gi", y = "Observed mean_gi") +
  theme_classic()


## glm 
glm(mean_gi ~ cov_mac_holdfast_live + 
      macro_stipe_density_m2 + 
      n_macro_plants_m2 + 
      cov_bare_sand + 
      purple_urchin_conceiledm2 + 
    , data = predic, family = gaussian)
summary(m3) #AIC: 371.21
cov_mac_holdfast_live    macro_stipe_density_m2         n_macro_plants_m2 
0.2891270                 0.2880061                 0.2824447 
cov_bare_sand purple_urchin_conceiledm2 
0.2411657                 0.2153091 

