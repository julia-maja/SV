library(rfishbase)
library(tidyverse)

### salinity ###

species_name <- "Danio rerio"

get_salinity <- function(species_name) {
  eco_data <- ecosystem(species_name)
  if (nrow(eco_data) == 0) {
    return(NA)
  }
  # filter out introduced species
  native_data <- eco_data %>%
    filter(tolower(Status) != "introduced")
  if (nrow(native_data) == 0) {
    return(NA)
  }
  
  salinities <- eco_data$Salinity
  
  freq_table <- table(salinities)
  max_freq <- max(freq_table)
  
  most_common <- names(freq_table[freq_table == max_freq])
  
  return(most_common)
}

#get_salinity("Danio rerio")

#SV_data_sub <- SV_data_avg[261:262,]
#SV_data_sub$salinity <- sapply(SV_data_sub$Species, get_salinity)

#SV_data_avg$salinity <- sapply(SV_data_avg$Species, get_latitude)

if (!"salinity" %in% names(SV_data_avg)) {
  SV_data_avg$salinity <- NA
}

# Loop through each row
for (i in seq_len(nrow(SV_data_avg))) {
  # Only compute if latitude is NA (i.e. not yet processed)
  if (is.na(SV_data_avg$salinity[i])) {
    species <- SV_data_avg$Species[i]
    
    # Try to get latitude, but handle potential errors
    result <- tryCatch(
      get_salinity(species),
      error = function(e) NA
    )
    
    # Assign result
    SV_data_avg$salinity[i] <- result
    
    # Optionally print progress
    cat("Processed row", i, "\n")
    
    # Save progress every few rows or every time
    if (i %% 10 == 0) {
      saveRDS(SV_data_avg, "SV_data_salinity_progress.rds")
    }
  }
}
SV_data_avg <- readRDS("SV_data_salinity_progress.rds")
write.csv(SV_data_avg, "/Users/juliamaja/Desktop/SV/SV_data_avg.csv") 

