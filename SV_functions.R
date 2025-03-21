
# function that uses ncbi taxonomy to retrieve the order given a species name
options(taxize_entrez_key = "ecfcbe684204381f87a6d00c9cb5c80a5a08")
get_order_ncbi <- function(species_name) {
  tax_data <- classification(species_name, db = "ncbi")
  tax_df <- tax_data[[1]]
  if ( any(is.na(tax_data[[1]])) ) {
    return(NA)  
  }
  order_row <- which(tax_df$rank == "order")
  if (length(order_row) > 0) {
    tax_order <- tax_df[order_row, 1]  
  } else {
    tax_order <- NA  
  }
  return(tax_order)
}

options(taxize_entrez_key = "ecfcbe684204381f87a6d00c9cb5c80a5a08")
get_family_ncbi <- function(species_name) {
  tax_data <- classification(species_name, db = "ncbi")
  tax_df <- tax_data[[1]]
  if ( any(is.na(tax_data[[1]])) ) {
    return(NA)  
  }
  order_row <- which(tax_df$rank == "family")
  if (length(order_row) > 0) {
    tax_order <- tax_df[order_row, 1]  
  } else {
    tax_order <- NA  
  }
  return(tax_order)
}

options(taxize_entrez_key = "ecfcbe684204381f87a6d00c9cb5c80a5a08")
get_genus_ncbi <- function(species_name) {
  tax_data <- classification(species_name, db = "ncbi")
  tax_df <- tax_data[[1]]
  if ( any(is.na(tax_data[[1]])) ) {
    return(NA)  
  }
  order_row <- which(tax_df$rank == "genus")
  if (length(order_row) > 0) {
    tax_order <- tax_df[order_row, 1]  
  } else {
    tax_order <- NA  
  }
  return(tax_order)
}

# function that finds the node number of the MCRA of each order
findMRCANode <- function(phylo, trait.data){
  nodes_list <- list()
  orders <- unique(trait.data$order)  
  
  for(current_order in orders){
    #ensure the species are in the tree you're working with
    trait.data.phylo <- trait.data[trait.data$tips %in% phylo$tip.label,]
    #remove any taxonomic levels with only one species (cannot find MRCA for one species)
    trait.data.phylo <- trait.data.phylo %>% group_by(order) %>% filter(n() > 1) %>% ungroup()
    #subset the trait data into species belonging to the same taxonomic group (current_order)
    trait.data.filtered <- trait.data.phylo %>% filter(order == current_order)
    #if no species in this order, skip this iteration
    if(nrow(trait.data.filtered) == 0){
      next
    }
    #take the vector of these species names (tips)
    tip_vector <- trait.data.filtered$tips  
    #find the node number of the MRCA of these species
    MRCA_node <- findMRCA(phylo, tip_vector)
    #create a dataframe for this order and its MRCA node
    loop_df <- data.frame(clade_name = current_order, node_number = MRCA_node)
    #add the result to the list
    nodes_list <- append(nodes_list, list(loop_df))
  }
  
  nodes_df <- do.call(rbind, nodes_list)
  return(nodes_df)
}