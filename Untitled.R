library(rentrez)
library(dplyr)  
library(XML)

species_data <- read.csv("/Users/juliamaja/Downloads/SV_sp_names_row.csv")

get_correct_name <- function(species_name) {
  # Query NCBI taxonomy database to find the species
  search_result <- entrez_search(db="taxonomy", term=species_name)
  
  # If a match is found, fetch the correct scientific name
  if (search_result$count > 0) {
    # Get the first result and retrieve the scientific name
    search_result <- entrez_search(db="taxonomy", term="Scleropages jardini")
    tax_id <- search_result$ids[1]
    tax_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml")
    tax_info_parsed <- xmlParse(tax_info)
    correct_name <- xpathSApply(tax_info_parsed, "//ScientificName", xmlValue)[1]
    return(correct_name)
  } else {
    # If no match is found, return NA (or the original name)
    return(NA)
  }
}


correct_names <- get_correct_name(species_data)

# View the corrected names
view(species_data)




#troubleshooting: 

# Error in entrez_check(response) :
#HTTP failure 414, the request is too large. For large requests, try using web history as described in the rentrez tutorial
# the tutorial says I need to pist my list of tax IDs to the NCBI server first.

#shortening the function to only retrieve taxids:

get_taxid <- function(species_name) {
  # Query NCBI taxonomy database to find the species
  search_result <- entrez_search(db="taxonomy", term=species_name)
  
  # If a match is found, fetch the correct scientific name
  if (search_result$count > 0) {
    # Get the first result and retrieve the scientific name
    search_result <- entrez_search(db="taxonomy", term="Scleropages jardini")
    tax_id <- search_result$ids[1]
   
    return(taxid)
  } else {
    # If no match is found, return NA (or the original name)
    return(NA)
  }
}
taxids <- get_taxid(species_data)

# this gives the same error but how else do I convert my species names to taxids?
