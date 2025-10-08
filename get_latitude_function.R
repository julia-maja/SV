library(tidyverse)
library(rfishbase)

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
SV_data_latitude <- SV_data_avg

SV_data_avg <- SV_data_avg %>% left_join(SV_data_latitude[, c("tips", "latitude")], by = "tips")

# which(names(test_ecosystem) == "NorthernLat")

get_latitude <- function(species_name) {
  eco_data <- ecosystem(species_name)
  if (nrow(eco_data) == 0) {
    return(NA)
  }
  # filter out introduced or questionable species
  native_data <- eco_data %>%
    filter(tolower(Status) != "introduced")
  if (nrow(native_data) == 0) {
    return(NA)
  }
  native_data <- eco_data %>%
    filter(tolower(Status) != "questionable")
  if (nrow(native_data) == 0) {
    return(NA)
  }
  
  north_lat <- native_data$NorthernLat
  south_lat <- native_data$SouthernLat
  
  Nlat <- c(north_lat)
  Slat <- c(-south_lat)
  lat <- c(Nlat, Slat)
  med_lat <- median(lat, na.rm = TRUE)
  
  return(abs(med_lat))
}

#SV_data_sub <- SV_data_avg[261:262,]
#SV_data_sub$latitude <- sapply(SV_data_sub$Species, get_latitude)

SV_data_latitude$latitude <- sapply(SV_data_latitude$Species, get_latitude)
write.csv(SV_data_latitude, "/Users/juliamaja/Desktop/SV/SV_data_latitude.csv") 

# Check for an existing saved column or create a placeholder
if (!"latitude" %in% names(SV_data_avg)) {
  SV_data_avg$latitude <- NA
}

# Loop through each row
for (i in seq_len(nrow(SV_data_avg))) {
  # Only compute if latitude is NA (i.e. not yet processed)
  if (is.na(SV_data_avg$latitude[i])) {
    species <- SV_data_avg$Species[i]
    
    # Try to get latitude, but handle potential errors
    result <- tryCatch(
      get_latitude(species),
      error = function(e) NA
    )
    
    # Assign result
    SV_data_avg$latitude[i] <- result
    
    # Optionally print progress
    cat("Processed row", i, "\n")
    
    # Save progress every few rows or every time
    if (i %% 10 == 0) {
      saveRDS(SV_data_avg, "SV_data_latitude_progress.rds")
    }
  }
}
SV_data_avg <- readRDS("SV_data_latitude_progress.rds")
write.csv(SV_data_avg, "/Users/juliamaja/Desktop/SV/SV_data_avg.csv") 
