# in this script I'm attempting to add climate data to my SV data

library(tidyverse)
library(rfishbase)

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

# some ingredients
# the function ecosystem()
ecosystem("Danio rerio")
# find the index of the column(s) you want
which(colnames(test_ecosystem) == "Climate")
# [1] 49
# test_ecosystem[49]
# 11 rows
# only keep ones where the column "Status" is "native"
# Status is 6th col
test_ecosystem[test_ecosystem$Status == "native", 49]
# now only 8 rows
# we need to keep only the most common climate, in this case "tropical" 
# make sure case insensitive!
climates <- test_ecosystem[test_ecosystem$Status == "native", 49]
climates$Climate <- tolower(climates$Climate)
freq_table <- table(climates)
max_freq <- max(freq_table)
most_common <- names(freq_table[freq_table == max_freq])

# creating the function:
get_climate <- function(species_name) {
  eco_data <- ecosystem(species_name)
  if (nrow(eco_data) == 0) {
    return(NA)
  }
  # filter for native  only
  native_data <- eco_data %>%
    filter(tolower(Status) != "introduced")
  if (nrow(native_data) == 0) {
    return(NA)
  }
  
  climates <- tolower(native_data$Climate)
  
  freq_table <- table(climates)
  max_freq <- max(freq_table)
  
  most_common <- names(freq_table[freq_table == max_freq])
  
  return(most_common)
}

# polar, boreal, temperate, subtropical, tropical

# now apply this function to my data frame to make a new column called climate

# test on a small subset 
SV_data_sub <- SV_data_avg[261:262,]
SV_data_sub$climate <- sapply(SV_data_sub$Species, get_climate)

# now apply it to the whole data frame
SV_data_avg$climate <- sapply(SV_data_avg$Species, get_climate)
# this is gonna take a hot minute 

write.csv(SV_data_avg, "/Users/juliamaja/Desktop/SV/SV_data_avg.csv") 

SV_data_climate <- SV_data_avg
SV_data_climate <- SV_data_climate %>% mutate(across(where(is.list), ~ sapply(., toString))) 
write.csv(as_data_frame(SV_data_climate), "/Users/juliamaja/Desktop/SV/SV_data_climate.csv") 

# perhaps save this function into my functions script and source it in data_cleaning
# and run it there? 

