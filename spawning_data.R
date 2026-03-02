library(rfishbase)
library(tidyverse)

# this script gathers relevant eco/ morphological data on species in my SV dataset, from fishbase, and
# combines it into one neat data frame

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
specieslist <- as.vector(SV_data_avg[,5])

# reproduction() function
reproduction <- (reproduction(species_list = specieslist, 
                                   fields = c("Species", "Spawning", "RepGuild1", "RepGuild2", "AddInfos") ))

reproduction_missing <- as.data.frame(setdiff(SV_data_avg$Species, reproduction_data$Species))
colnames(reproduction_missing)[1] <- "Species"

# join the dataframe of reproduction data and missing species names to create a data frame
# of all my species including those with missing reproduction data

reproduction_data <- bind_rows(reproduction, reproduction_missing)
reproduction_data$Spawning <- tolower(reproduction_data$Spawning)

write_csv(reproduction_data, "/Users/juliamaja/Desktop/SV/SV_data_spawning.csv")
