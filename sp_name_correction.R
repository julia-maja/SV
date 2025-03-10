library(rentrez)
library(dplyr)  
library(XML)

species_data <- read.csv("/Users/user2/Desktop/SV_sp_names_row.csv") 
#species_data <- c("Arapaima gigas","Heterotis niloticus","Osteoglossum bicirrhosum","Scleropages jardini","Pantodon buchholzi","Notopterus chitala")

#species_names <- species_data

# Function to get the correct scientific names for a (short) list of species names
get_correct_name <- function(species_names) {
  corrected_names <- vector("list", (nrow(species_names))/4 ) #
  for (i in 1:nrow(species_names)) {
    species_name <- species_names[i,]
    search_result <- entrez_search(db = "taxonomy", term = species_name)
    if (search_result$count > 0) {
      tax_id <- search_result$ids[1]
      tax_info <- entrez_fetch(db = "taxonomy", id = tax_id, rettype = "xml")
      tax_info_parsed <- xmlParse(tax_info)
      correct_name <- xpathSApply(tax_info_parsed, "//ScientificName", xmlValue)[1]
      corrected_names[[i]] <- correct_name
    } else {
      corrected_names[[i]] <- NA
    }
  }
  return(corrected_names)
}

corrected_species_names <- get_correct_name(species_data)

corrected_species_names <- as.data.frame(t(corrected_species_names)) 
rownames(corrected_species_names) <- NULL
updated_species_data <- bind_cols(species_data, corrected_species_names)

sample <- updated_species_data[c(1,2,4,31,33,35,45), ]

################################################################################
# get_correct_name <- function(species_name) {
#   corrected_names <- vector("list", length(species_name))
#   
#   for (i in 1:length(species_name)) {
#     species_name <- species_name[i]
#   
#   # If a match is found, fetch the correct scientific name
#   if (search_result$count > 0) {
#     # Get the first result and retrieve the scientific name
#     search_result <- entrez_search(db="taxonomy", term=species_name[i])
#     tax_id <- search_result$ids[1]
#     tax_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml")
#     tax_info_parsed <- xmlParse(tax_info)
#     correct_name <- xpathSApply(tax_info_parsed, "//ScientificName", xmlValue)[1]
#     corrected_names[[i]] <- correct_name
#     
#   } else {
#     # If no match is found, return NA (or the original name)
#     return(NA)
#   }
#   }
# }
# 
# 
# correct_names <- get_correct_name(species_data)
# 
# # View the corrected names
# view(species_data)

library(rentrez)
library(XML)

# Function to get the correct scientific names for a list of species names
get_correct_name <- function(species_names) {
  # Create an empty list to store the corrected names
  corrected_names <- vector("list", length(species_names))
  
  # Loop through each species name in the list
  for (i in 1:length(species_names)) {
    species_name <- species_names[i]
    
    # Query NCBI taxonomy database to find the species
    search_result <- entrez_search(db = "taxonomy", term = species_name)
    
    # If a match is found, fetch the correct scientific name
    if (search_result$count > 0) {
      # Get the first result and retrieve the scientific name
      tax_id <- search_result$ids[1]
      tax_info <- entrez_fetch(db = "taxonomy", id = tax_id, rettype = "xml")
      tax_info_parsed <- xmlParse(tax_info)
      correct_name <- xpathSApply(tax_info_parsed, "//ScientificName", xmlValue)[1]
      
      # Store the corrected name in the list
      corrected_names[[i]] <- correct_name
    } else {
      # If no match is found, store NA
      corrected_names[[i]] <- NA
    }
  }
  
  # Return the list of corrected names
  return(corrected_names)
}

# Example usage with a vector of species names
species_list <- c("Scleropages jardini", "Panthera leo", "Homo sapiens")
corrected_species_names <- get_correct_name(species_list)
print(corrected_species_names)


