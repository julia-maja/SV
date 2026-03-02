# combine all datasets into one
library(tidyverse)

SV_data_avg <- read_csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
climate <- read_csv("/Users/juliamaja/Desktop/SV/SV_data_climate.csv")
spawning <- read_csv("/Users/juliamaja/Desktop/SV/SV_data_spawning.csv")
#salinity <- 
#latitude <- 
# salinity and latitude data are already in SV_data_avg

# > setdiff(SV_data_avg$Species,climate$Species)
# [1] "Corydoras acutus"       
# [2] "Trichogenes longipinnis"
# [3] "Anabas testudineus" 
# not sure why these species are not included in climate data
climate <- climate %>% select("Species", "climate")
source("/Users/juliamaja/Desktop/SV/get_climate.R")
get_climate("Corydoras acutus") # tropical
get_climate("Trichogenes longipinnis") # tropical
get_climate("Anabas testudineus") # tropical
new_rows <- data.frame(
  Species = c("Corydoras acutus", "Trichogenes longipinnis", "Anabas testudineus"),
  climate = c("tropical", "tropical", "tropical"))
climate <- rbind(climate, new_rows)

# get feeding type data
library(rfishbase)
specieslist <- as.vector(SV_data_avg[,5])
feed_type <- ecology(specieslist, c("Species", "FeedingType"))
# which(duplicated(feed_type$Species))
# feed_type[266,]
# Oreochromis niloticus is duplicated, so remove the NA entry
feed_type <- feed_type[-266,]

# diel data
diel_data <- readRDS("/Users/juliamaja/Desktop/SV/trait_data_fish.rds")
diel_data <- (diel_data 
              %>% select("species", "diel", "diel1", "diel2", "diel_continuous") 
              %>% rename("Species" = "species")
              %>% mutate(Species = gsub("_", " ", Species))
              )

SV_data_eco <- SV_data_avg %>% left_join(climate, by = "Species")
SV_data_eco <- SV_data_eco %>% left_join(spawning, by = "Species")
SV_data_eco <- SV_data_eco %>% left_join(feed_type, by = "Species")
SV_data_eco <- SV_data_eco %>% left_join(diel_data, by = "Species")

SV_data_eco <- SV_data_eco[, -c(1:4)] # removing the first 4 colums that are just numbering the rows
  
write.csv(SV_data_eco, "/Users/juliamaja/Desktop/SV/SV_data_eco.csv") 
