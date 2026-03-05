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