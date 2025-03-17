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


url <- 'https://docs.google.com/spreadsheets/d/1Z2P6dAcoU0-Kh0UqVgIAeongX5IJohuMUKG1PyrNraU/edit?gid=84578656#gid=84578656'
sp_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
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
# SV_data <- (sp_data 
#   %>% mutate(species = str_replace(species, " ", "_"))
#   %>% group_by(species) 
#   %>% mutate(
#     genome.assembly = if_else(
#       genome.assembly == '' & any(genome.assembly == 'y'),  
#       'y',  
#       genome.assembly 
#     )
#   ) 
#   %>% ungroup()
#   %>% mutate(species = tolower(species))
#   %>% left_join(resolved_names, by = "species")
#   %>% select("species", "unique_name", "order", "SV", "genome.assembly", "presence", "ott_id", "flags")
#   %>% mutate(tips = unique_name)
#   %>% filter(source != "")
#   %>% mutate(presence = ifelse(presence == "" && SV == "0", "absent", "present"))
#   %>% mutate(SV = str_replace(SV, "2", "1"), str_replace(SV, "0.5", "1"), str_replace(SV, "2.5", "2"), str_replace(SV, "3.5", "4"))
#   %>% filter(unique_name !="")
#   %>% select("species", "unique_name", "order", "SV", "genome.assembly", "ott_id", "flags", "tips")
# )

# add existing orders from previous runs of this script:
SV_data <- SV_data %>% left_join(Species_order, by = c("unique_name" = "Species")) %>% mutate(order.x = ifelse(order.x == "", order.y, order.x)) %>% mutate(order.x = ifelse(is.na(order.x), order.y, order.x)) %>% rename("order" = "order.x") %>% select("species", "unique_name", "order", "SV", "genome.assembly", "ott_id", "flags", "tips") 

# taking the median across species with multiple entries 
# when a species has an entry with a 0 but other entries > 0, ignore the 0
# eventually implement data quality scores to get a more reliable average
SV_data_avg <- (SV_data %>% group_by(unique_name) 
            %>% mutate(unique_name = as.character(unique_name), SV = as.numeric(SV))
            %>% mutate(nonzero_exists = any(SV > 0)) 
            %>% filter(!(nonzero_exists & SV == 0)) 
            %>% summarise(SV = median(SV, na.rm = TRUE)) 
            %>% ungroup()
            %>% mutate(SV = as.character(SV))
            %>% left_join(SV_data, by = "unique_name")
            %>% select("unique_name", "SV.x", "order", "genome.assembly", "ott_id", "flags", "tips")
            #%>% rename("SV", "SV.x")
)

SV_data_avg$SV <- SV_data_avg$SV.x
SV_data_avg <- (SV_data_avg %>% select("unique_name", "SV", "order", "genome.assembly", "ott_id", "flags", "tips")
                %>% group_by(unique_name) %>% distinct() %>% ungroup()
                )
# SV_data_avg does not retain SV values of "present" because of as.numeric()

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

SV_data_orders <- SV_data_avg %>% mutate(order = ifelse(is.na(order), purrr::map_chr(unique_name, get_order), order))


# filling in missing orders using fishbase --------------------------------
fishbase_species <- rfishbase::load_taxa()
species <- SV_data_orders$unique_name
fishbase_species_sub <- fishbase_species[fishbase_species$Species %in% species,]
# missing 215 species! ncbi was only missing 25

# filling in those 25 with fishbase:
SV_data_orders_NA <- SV_data_orders[is.na(SV_data_orders$order), ]
species <- SV_data_orders_NA$unique_name
ncbi_missing_orders_in_fishbase <- fishbase_species[fishbase_species$Species %in% species,]
ncbi_missing_orders_in_fishbase <- ncbi_missing_orders_in_fishbase %>% select("Species", "Order")

SV_data_orders_complete <- SV_data_orders %>% rename("Species" = "unique_name") %>% left_join(ncbi_missing_orders_in_fishbase, by = "Species") %>% mutate(order = ifelse(is.na(order), Order, order))
# use equal instead of comma for rename(). 

# make a reference table for species and orders so you don't have to keep running th eget_order() function all over again each time you increase the dataset:
Species_order <- SV_data_orders_complete %>% group_by(Species) %>% distinct() %>% select("Species", "order")
write.csv(Species_order, "/Users/user2/Desktop/SV_project/Species_order.csv")



# change SV_data_avg so that it retains "present" values in SV



SV_data_presence <- (SV_data
              %>% mutate(presence = ifelse(SV != 0, "present", "absent"))
              %>% group_by(tips)
              %>% select("unique_name", "order", "tips", "presence", "ott_id", "flags", "genome.assembly")
              %>% distinct()
              %>% ungroup()
        )

# tree building for SV complexity:
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)

# tree bulding for SV presence/ absence:
tr <- tol_induced_subtree(ott_ids = SV_data_presence$ott_id[SV_data_presence$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data_presence$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)

# Section 2: Creating the SV plots w OTL --------------------------------

display.brewer.all(type = "qual", colorblindFriendly = TRUE)
display.brewer.all(type = "seq")
display.brewer.all(type = "div")
display.brewer.all(type = "all")

# SV complexity
sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot 

# SV presence/ absence
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data_presence[, c("tips", "presence", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette) 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot 

node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node),colour = "blue", hjust=-.2, size = 2) + geom_tiplab(size = 2, hjust = -0.1)
node_labels


# labelling families/ orders ----------------------------------------------
# find mrca of a set of species

# amelia's function:
# findMRCANode <- function(phylo, trait.data, taxonomic_level_col){
#   nodes_list <- list()
#   for(i in unique(trait.data[,taxonomic_level_col])){
#     #ensure the species are in the tree you're working with
#     trait.data <- trait.data[trait.data$tips %in% phylo$tip.label,]
#     #remove any taxonomic levels with only one species (cannot find MRCA for one species)
#     trait.data <- trait.data %>% group_by_at(taxonomic_level_col) %>% filter(n()>1)
#     #subset the trait data into species belonging to the same taxonomic group
#     trait.data.filtered <- trait.data[trait.data[,taxonomic_level_col] == i, ]
#     #take the vector of these species names
#     tip_vector <- as.vector(trait.data.filtered[, "tips"])
#     #find the node number of the MRCA of these species
#     MRCA_node <- findMRCA(phylo, tip_vector$tips)
#     #create a dataframe
#     loop_df <- data.frame(clade_name = i, node_number = MRCA_node)
#     #add to list, to extract later (out of for loop)
#     nodes_list[[i]] <- loop_df}
#   
#   nodes_df <- do.call(rbind, nodes_list)
#   return(nodes_df)
# }

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

nodes <- findMRCANode(phylo = tr, trait.data = SV_data_orders_complete)


all_nodes <- lapply(unique(SV_data_orders_complete$order), function(x) findMRCANode(phylo = tr, trait.data = SV_data_orders_complete))
nodes_df <- do.call(rbind, all_nodes)
nodes_df <- nodes_df %>% distinct()

#can now easily label all clades within that taxonomic level on the tree using the nodes_df
num_orders <- SV_data_orders_complete %>% group_by(order) %>% distinct() %>% count()
set.seed(30) 
color_palette <- sample(colors(), num_orders)  
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








  
