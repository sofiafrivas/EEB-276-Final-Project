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

library(MASS)
library(ggeffects)
library(brms)
library(cmdstanr)

#don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

# Read in Dataframe(s) ----------------------------------------------------

gonadindex_raw <- readRDS("quad_build3.rds")
##write.csv(gonadindex_raw, "gonadindex_raw.csv", row.names = FALSE)

gi_predictors <- gonadindex_raw %>%
  mutate(site_id = paste(site, pred_patch, zone, year, sep = " ")) %>%
  dplyr::select(-c(patch_id, latitude, longitude, survey_date, site, site_type, patch_cat, 
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

brm_formula <- bf(
  mean_gi ~ purple_urchin_densitym2 +
    juveniles +
    lamr +
    macr +
    cov_mac_holdfast_live +
    cov_crustose_coralline
)

priors <- c(
  prior(normal(0, 2), class = "b"),        # slopes
  prior(normal(8, 5), class = "Intercept"),
  prior(student_t(3, 0, 5), class = "sigma")
)

brm_model <- brm(
  formula = brm_formula,
  data = gi_predictors,
  family = gaussian(),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 42
)
summary(brm_model)
plot(brm_model)
# Practice Plots ----------------------------------------------------------

library(ggplot2)
install.packages("ggfortify")
library(ggfortify)

# Run PCA (on numeric predictors only, no NAs)
pca <- prcomp(na.omit(gi_predictors %>% dplyr::select(where(is.numeric))), 
              scale. = TRUE)

# Quick plot
autoplot(pca, data = na.omit(gi_predictors), 
         colour = "pred_patch",  # color by a grouping variable
         loadings = TRUE, 
         loadings.label = TRUE)

# PC1 vs PC3
autoplot(pca, x = 1, y = 3, data = na.omit(gi_predictors),
         colour = "site_id",
         loadings = TRUE,
         loadings.label = TRUE)

# PC2 vs PC3
autoplot(pca, x = 2, y = 3, data = na.omit(gi_predictors),
         colour = "site_id",
         loadings = TRUE,
         loadings.label = TRUE)

pca$rotation[, 1:3]
sort(abs(pca$rotation[,1]), decreasing = TRUE) |> head(10)

pc_scores <- as.data.frame(pca$x[, 1:8]) 

pc_data <- pc_scores %>%
  mutate(mean_gi = na.omit(gi_predictors)$mean_gi)

pc_data$site_id <- na.omit(gi_predictors)$site_id

model_brms <- brm(mean_gi ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                    (1 | site_id),
                  data = pc_data,
                  family = gaussian(),
                  chains = 4, iter = 2000)

#this one looks best rn 
model_brms <- brm(mean_gi ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8,
                  data = pc_data,
                  family = gaussian(),
                  chains = 4, iter = 4000)

summary(model_brms)
r2(model_brms)

sort(abs(pca$rotation[,2]), decreasing = TRUE) |> head(10)

library(bayesplot)
mcmc_areas(model_brms, pars = vars(starts_with("b_")))
plot(model_brms)

pp_check(model_brms)


loadings <- as.data.frame(pca$rotation[, c(1,3,5,6,8)])

for(pc in c("PC1","PC3","PC5","PC6","PC8")) {
  cat("\n---", pc, "---\n")
  cat("TOP POSITIVE:\n")
  print(head(sort(loadings[,pc], decreasing = TRUE), 5))
  cat("TOP NEGATIVE:\n")
  print(head(sort(loadings[,pc], decreasing = FALSE), 5))
}

library(corrplot)
corrplot(cor(gi_predictors %>% dplyr::select(where(is.numeric)), use = "complete.obs"))


install.packages("factoextra")
library(factoextra)

fviz_eig(pca)
