require(librarian)
librarian::shelf(tidyverse, here, janitor, googlesheets4, lubridate, splitstackshape,
                 dplyr,ggplot2,pwr2,tidyr, broom, ggpubr, paletteer)

#don't run these two lines 
load(file.path("/Volumes/enhydra/data/students/sofia/zone_level_data.rda"))
saveRDS(quad_build3, file = "quad_build3.rds")

gonadindex_predictors <- readRDS("quad_build3.rds")


write.csv(data.frame(column_names = colnames(gonadindex_predictors)),
          "column_names.csv",
          row.names = FALSE)
