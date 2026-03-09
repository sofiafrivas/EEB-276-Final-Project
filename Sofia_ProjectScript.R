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

predictors <- gi_predictors %>% 
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


summary(pca)$importance[3,] 

pc_scores1 <- as.data.frame(pca$x[, 1:12])
pc_data1 <- pc_scores1 %>%
  mutate(mean_gi = predictors_simple$mean_gi)
glm(mean_gi ~ ., data = pc_data)

summary(pca)
pve <- pca$sdev^2 / sum(pca$sdev^2)
round(pve, 3)



pc_scores <- as.data.frame(pca$x[, 1:12])

pc_data <- pc_scores %>%
  mutate(mean_gi = predictors$mean_gi)
glm(mean_gi ~ ., data = pc_data)

biplot(pca, scale = 0)
library(factoextra)

fviz_pca_var(pca,
             col.var = "contrib",
             gradient.cols = c("lightblue","blue","darkblue"),
             repel = TRUE) ## look into this 

#combines pca and horseshoe
posterior <- as_draws_df(brm_hs)
coef_cols <- grep("^b_", colnames(posterior), value = TRUE)
coef_cols <- coef_cols[coef_cols != "b_Intercept"]
coef_summary <- posterior %>%
  summarise(across(all_of(coef_cols), ~ median(.))) %>%
  pivot_longer(everything(), names_to = "predictor", values_to = "median") %>%
  mutate(
    predictor_name = gsub("^b_", "", predictor),
    median = as.numeric(median)
  )
# Select only meaningful predictors
top_predictors <- coef_summary %>%
  filter(abs(median) > 0.05) %>%
  pull(predictor_name)
# 2. Create color vector for PCA variable plot
var_colors <- ifelse(rownames(pca$rotation) %in% top_predictors, "Top", "Other")
# 3. Plot PCA variables with top predictors highlighted
fviz_pca_var(
  pca,
  col.var = var_colors,      # discrete values: "Top" vs "Other"
  repel = TRUE
) +
  scale_color_manual(values = c("Top" = "red", "Other" = "grey70")) +
  labs(title = "PCA variable plot highlighting top horseshoe predictors")




install.packages("factoextra")
library(factoextra)


#setting priors 
priors_pca <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(6.5, 2), class = "Intercept"),
  prior(student_t(3, 0, 2), class = "sigma"))

#scree plot
library(factoextra)
fviz_eig(pca)
plot(pca, type = "l")



predictors_scaled <- predictors %>%
  mutate(across(-mean_gi, scale))

brm_hs <- brm(
  mean_gi ~ .,
  data = predictors_scaled,
  prior = prior(horseshoe(1), class = "b"),
  control = list(adapt_delta = 0.99)
)
summary(brm_hs)
plot(brm_hs)
fixef(brm_hs)
pp_check(brm_hs)
conditional_effects(brm_hs)
bayes_R2(brm_hs)
library(bayesplot)
posterior <- as_draws_df(brm_hs)
coef_cols <- grep("^b_", colnames(posterior), value = TRUE)
coef_cols <- coef_cols[coef_cols != "b_Intercept"]  # remove intercept
mcmc_intervals(posterior, pars = coef_cols)
coef_summary <- posterior %>%
  summarise(across(all_of(coef_cols), ~ median(.))) %>%
  pivot_longer(everything(), names_to = "predictor", values_to = "median")
top_predictors <- coef_summary %>%
  filter(abs(median) > 0.05) %>%  # threshold can be adjusted
  pull(predictor) %>%
  gsub("^b_", "", .)  # remove b_ prefix to match PCA loadings
loadings <- as.data.frame(pca$rotation)
top_loadings <- loadings[top_predictors, 1:4]  # first 4 PCs
top_loadings$predictor <- rownames(top_loadings)
top_loadings
long_loadings <- top_loadings %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading")

# Merge with effect size from horseshoe
long_loadings <- long_loadings %>%
  left_join(coef_summary %>% mutate(predictor = gsub("^b_", "", predictor)), by = "predictor")

# nope
ggplot(long_loadings, aes(x = PC, y = loading, fill = median)) +
  geom_bar(stat = "identity") +
  facet_wrap(~predictor, scales = "free_y") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme_classic() +
  labs(
    title = "Top predictor loadings on PCs with horseshoe effect size",
    y = "PCA loading",
    fill = "Median posterior"
  )


#finding r2 from only significant horseshoe variables
sig_vars <- posterior_summary(brm_hs) %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  filter(grepl("^b_", parameter),
         parameter != "b_Intercept") %>%
  filter(Q2.5 > 0 | Q97.5 < 0) %>%
  pull(parameter) %>%
  gsub("^b_", "", .)
# Build formula dynamically
sig_formula <- as.formula(paste("mean_gi ~", paste(sig_vars, collapse = " + ")))

model_hs_reduced <- brm(sig_formula,
                        data = predictors_scaled,
                        family = gaussian(),
                        prior = c(
                          prior(normal(0, 1), class = "b"),
                          prior(normal(6.5, 2), class = "Intercept"),
                          prior(student_t(3, 0, 2), class = "sigma")
                        ),
                        chains = 4, iter = 4000, warmup = 2000,
                        cores = 4)
bayes_R2(model_hs_reduced)
loo_hs_reduced <- loo(model_hs_reduced)
loo_compare(loo_skew, loo_hs_reduced)



fviz_pca_var(
  pca,
  col.var = var_colors,      # discrete values: "Top" vs "Other"
  repel = TRUE) +
  scale_color_manual(values = c("Top" = "red", "Other" = "grey70")) +
  labs(title = "PCA variable plot highlighting top horseshoe predictors")




## glm 
glm_pca <- glm(mean_gi ~ cov_mac_holdfast_live + 
      macro_stipe_density_m2 + 
      n_macro_plants_m2 + 
      cov_bare_sand + 
      purple_urchin_conceiledm2  
    , data = predictors_simple, family = gaussian)
summary(glm_pca)
r2(glm_pca)

glm(mean_gi ~ PC1 + PC2 + PC3 + PC4, data = pc_data)
brm(mean_gi ~ PC1 + PC2 + PC3 + PC4, data = pc_data)


# predictors 2.0  ---------------------------------------------------------

predictors2 <- predictors %>%
  dplyr::select(mean_gi, cov_mac_holdfast_live, macro_stipe_density_m2, n_macro_plants_m2, cov_bare_sand, purple_urchin_conceiledm2, cov_crustose_coralline, purple_urchin_densitym2, relief_cm, marine_snails, red_urchin_densitym2, cov_corynactis_californica, cov_barnacle, cov_articulated_coralline, cov_desmarestia_spp, cov_colonial_tunicate, density_lamstump, cov_dead_kelp_holdfast_any, cov_dictyoneurum_spp, red_urchin_conceiledm2)

pca2 <- prcomp(predictors2 %>% 
                dplyr::select(-mean_gi), 
              scale. = TRUE) 

#visualize 
fviz_eig(pca2)

pc_scores2 <- as.data.frame(pca2$x[, 1:10])

#add response variable (mean_gi) into new dataframe
pc_data2 <- pc_scores2 %>%
  mutate(mean_gi = predictors2$mean_gi)

#set priors for brm 
priors_pca2 <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(6.5, 2), class = "Intercept"), #mean mean_gi is ~6.5
  prior(student_t(3, 0, 2), class = "sigma"))

model_pca_skew2 <- brm(mean_gi ~ .,  
                      data = pc_data2,
                      family = skew_normal(), 
                      prior = priors_pca2,
                      chains = 4, iter = 4000, warmup = 2000,
                      cores = 4)
summary(model_pca_skew2)
fviz_pca_var(
  pca2,
  col.var = "contrib",  # color by contribution to the PC axes
  repel = TRUE,
  gradient.cols = c("grey80", "steelblue", "darkblue")
) +
  labs(title = "PCA variable plot")

fviz_pca_var(pca2, axes = c(1, 3),
             col.var = "contrib",
             repel = TRUE,
             gradient.cols = c("grey80", "steelblue", "darkblue")) +
  labs(title = "PCA variable plot - PC1 vs PC3")

fviz_pca_var(pca2, axes = c(1, 2),
             select.var = list(contrib = 10),  # only top 10 contributors
             col.var = "contrib",
             repel = TRUE,
             gradient.cols = c("grey80", "steelblue", "darkblue")) +
  labs(title = "PCA variable plot")+
  theme_classic()

fviz_pca_var(pca, axes = c(1, 3),
             select.var = list(contrib = 10),  # only top 10 contributors
             col.var = "contrib",
             repel = TRUE,
             gradient.cols = c("grey80", "steelblue", "darkblue")) +
  labs(title = "PCA variable plot - PC1 vs PC3")+
  theme_classic()


for(pc in c("PC1","PC2","PC3")) { #starts a for loop 
  cat("\n---", pc, "---\n") #syntax for output: headers for each PC
  cat("TOP POSITIVE:\n") #syntax for output: label for the variables with  largest positive loadings for the current PC
  print(head(sort(pca2$rotation[,pc], decreasing = TRUE), 5)) #sorts loadings from largest to smallest, keeps top 5 variables 
  cat("TOP NEGATIVE:\n") #syntax for output: label for the variables with  largest negative loadings for the current PC
  print(head(sort(pca2$rotation[,pc], decreasing = FALSE), 5))} #sorts loadings from largest to smallest, keeps top 5 variables 


p1 + p2 + p3
