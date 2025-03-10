# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(ggplot2)
library(rphylopic)
library(RColorBrewer)
library(ggimage)
library(rfishbase)
library(rentrez)
library(taxize)

url <- 'https://docs.google.com/spreadsheets/d/1-5vhk_YOV4reklKyM98G4MRWEnNbcu6mnmDDJkO4DlM/edit#gid=0'
cetaceans_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

sp_data <- read.csv("/Users/user2/Desktop/SV_project/saccus_vasculosus.csv")
trait.data <- sp_data[,c("species", "SV", "order")]
trait.data$tips <- trait.data$species
trait.data$tips <- str_replace(trait.data$species, pattern = " ", replacement = "_")

otol_names <- tnrs_match_names(names = trait.data$species, context_name = "Vertebrates", do_approximate_matching = TRUE) 
resolved_names <- (otol_names[!(is.na(otol_names$unique_name)),] 
  %>%  rename( "species" = "search_string") 
  %>% mutate(species = str_replace(species, " ", "_"))
  %>% select(species, unique_name, "ott_id", "flags")
    )

# ensuring all instances of a species with genome assemblies have a "y"
# and shifting the SV classification scale
SV_data_1 <- read.csv("/Users/user2/Desktop/SV_project/saccus_vasculosus.csv")
SV_data_1 <- (SV_data_1 %>% filter(SV != "") 
  %>% mutate(species = str_replace(species, " ", "_"))
  %>% group_by(species) 
  %>% mutate(
    genome.assembly = if_else(
      genome.assembly == '' & any(genome.assembly == 'y'),  
      'y',  
      genome.assembly 
    )
  ) 
  %>% ungroup()
  %>% mutate(species = tolower(species))
  %>% left_join(resolved_names, by = "species")
  %>% select("species", "unique_name", "order", "SV", "genome.assembly", "ott_id", "flags")
  %>% mutate(tips = unique_name)
  %>% filter(SV != "no data")
  %>% mutate(SV = as.numeric(SV))
  %>% mutate(SV = case_when(
    SV == 2 ~ 1,
    SV == 0.5 ~ 1,
    SV == 2.5 ~ 2,
    SV == 3.5 ~ 4,
    TRUE ~ SV))
  %>% mutate(SV = as.character(SV))
)


# adding orders -----------------------------------------------------------
options(taxize_entrez_key = "ecfcbe684204381f87a6d00c9cb5c80a5a08")
get_order <- function(species_name) {
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

# #testing:
# SV_data <- SV_data[1:5,]
# SV_data_order <- SV_data %>% mutate(order = purrr::map_chr(unique_name, get_order)) 
# yippee! it works

SV_data <- SV_data_1 %>% mutate(order = purrr::map_chr(unique_name, get_order))

# problems:
# result <- get_order("Dicentrarchus labrax") # no "order" in rank column #solved
# result <- get_order("Heterocharax leptogrammus") # tax_data[[1]] is NA #solved

# some counts of the data:
all_species <- SV_data_1 %>% distinct(species) %>% count()
no_data <- SV_data_1 %>% filter(SV == "no data") %>% distinct(species) %>% count()
with_data <- SV_data_1 %>% filter(SV != "no data") %>% filter(SV != "") %>% distinct(species) %>% count()
searched <- no_data + with_data
with_genomes <- SV_data_1 %>% filter(genome.assembly == "y") %>% distinct(species) %>% count()
with_genomes_and_data <- SV_data_1 %>% filter(SV != "no data", genome.assembly == "y") %>% filter(SV != "", genome.assembly == "y") %>% distinct(species) %>% count()
with_genomes_no_data <- SV_data_1 %>% filter(SV == "no data", genome.assembly == "y") %>% distinct(species) %>% count()

summary_table <- data.frame(
  Category = c("all species", "no data", "with data", "searched", "with genomes", "with genomes and data", "with genomes and no data"),
  Count = c(all_species$n, no_data$n, with_data$n, searched$n, with_genomes$n, with_genomes_and_data$n, with_genomes_no_data$n)
)


# tree building:
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)


# Section 2: Creating the SV plots w OTL --------------------------------

display.brewer.all(type = "qual", colorblindFriendly = TRUE)
display.brewer.all(type = "seq")
display.brewer.all(type = "div")
display.brewer.all(type = "all")

sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette) 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot 

node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node),colour = "blue", hjust=-.2, size = 2) + geom_tiplab(size = 2, hjust = -0.1)
node_labels


# labelling families/ orders ----------------------------------------------
# find mrca of a set of species

# amelia's function:
findMRCANode <- function(phylo, trait.data, taxonomic_level_col){
  nodes_list <- list()
  for(i in unique(trait.data[,taxonomic_level_col])){
    #ensure the species are in the tree you're working with
    trait.data <- trait.data[trait.data$tips %in% phylo$tip.label,]
    #remove any taxonomic levels with only one species (cannot find MRCA for one species)
    trait.data <- trait.data %>% group_by_at(taxonomic_level_col) %>% filter(n()>1)
    #subset the trait data into species belonging to the same taxonomic group
    trait.data.filtered <- trait.data[trait.data[,taxonomic_level_col] == i, ]
    #take the vector of these species names
    tip_vector <- as.vector(trait.data.filtered[, "tips"])
    #find the node number of the MRCA of these species
    MRCA_node <- findMRCA(phylo, tip_vector$tips)
    #create a dataframe
    loop_df <- data.frame(clade_name = i, node_number = MRCA_node)
    #add to list, to extract later (out of for loop)
    nodes_list[[i]] <- loop_df}
  
  nodes_df <- do.call(rbind, nodes_list)
  return(nodes_df)
}

# my attempt at making it work:
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

nodes <- findMRCANode(phylo = tr, trait.data = SV_data)


all_nodes <- lapply(unique(SV_data$order), function(x) findMRCANode(phylo = tr, trait.data = SV_data))
nodes_df <- do.call(rbind, all_nodes)
nodes_df <- nodes_df %>% distinct()

#can now easily label all clades within that taxonomic level on the tree using the nodes_df
set.seed(30) 
color_palette <- #sample(colors(), 41)  
nodes_df$colour <- color_palette
order_tree <- sv.plot + geom_cladelab(node = nodes_df$node_number,
                                      label = nodes_df$clade_name,
                                      offset = 0.5,
                                      offset.text = 3,
                                      align = TRUE,
                                      fontsize = 2.5,
                                      barsize = 2,
                                      barcolour = nodes_df$colour
                                      )
order_tree








  
