# 1.1 GLM -------------------------------------------------------------------


#generalized linear models 
m1 <- glm(mean_gi ~ mean_gonad_mass_g +
            total_biomass_g +
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
AIC.m2 <- stepAIC(m2,
                  direction = "both", 
                  trace = TRUE)

m3 <- glm(mean_gi ~ mean_gonad_mass_g + 
            total_biomass_g + 
            purple_urchin_densitym2 + 
            juveniles + 
            lamr + 
            n_macro_plants_m2
          , data = gi_predictors, family = gaussian)
summary(m3) #AIC: 371.21
r2(m3) #0.668
AIC.m3 <- stepAIC(m3,
                  direction = "both", 
                  trace = TRUE)

m4 <- glm(mean_gi ~ purple_urchin_densitym2 + #ecological drivers only 
            juveniles +
            lamr +
            n_macro_plants_m2 +
            macr +
            cov_mac_holdfast_live +
            cov_crustose_coralline,
          , data = gi_predictors, family = gaussian)
summary(m4)#AIC: 421.07
r2(m4) #0.436
stepAIC(m4, direction="both")

m5 <- glm(mean_gi ~ purple_urchin_densitym2 +
            lamr +
            n_macro_plants_m2 +
            macr +
            adults 
          , family = gaussian, data = gi_predictors)
summary(m5)#AIC: 413.67
r2(m5) #0.457
stepAIC(m5, direction="both")

# 1.2 BRM -----------------------------------------------------------------

prior_brm <- prior(horseshoe(), class = "b")

model <- brm(
  mean_gi ~ .,
  data = gi_numeric,
  family = gaussian(),
  prior = prior_brm,
  chains = 4, cores = 4
)


# 1.3 PCA -----------------------------------------------------------------

#principal component analysis 

#use only numeric columns 
gi_numeric <- gi_predictors %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit()

#ATTEMPT 1 
pca <- prcomp(gi_numeric %>% dplyr::select(-mean_gi), 
              scale. = TRUE) 

library(factoextra)
fviz_eig(pca) #scree plot 

pc_scores <- as.data.frame(pca$x[, 1:7]) #number based on scree plot

pc_data <- pc_scores %>%
  mutate(mean_gi = gi_numeric$mean_gi)  # add response variable back

priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(6.5, 2), class = "Intercept"),
  prior(student_t(3, 0, 2), class = "sigma"))

model_prior <- brm(mean_gi ~ .,
                   data = pc_data,
                   family = gaussian(),
                   prior = priors,
                   sample_prior = "only",
                   chains = 4, iter = 2000)

pp_check(model_prior)  # should look like plausible mean_gi values

#ATTEMPT 2
pca2 <- prcomp(gi_numeric, scale. = TRUE)

summary(pca2)

loadings <- pca2$rotation

round(loadings[, 1:4], 2)

strong_loadings <- loadings[,1:4]
strong_loadings[abs(strong_loadings) >= 0.4]
apply(loadings[,1:4], 2, function(x) names(x[abs(x) >= 0.25]))

cor_mat <- cor(gi_numeric)
summary(as.vector(cor_mat[upper.tri(cor_mat)]))


library(car)
vif(m2)



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

