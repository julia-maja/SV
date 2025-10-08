library(rotl)
library(stringr)
library(ape)
library(geiger)
library(dplyr)


setwd("/Users/juliamaja/Desktop/SV")


SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
SV_data <- SV_data_avg #%>% filter(genome.assembly == "y") 
SV_data <- SV_data %>% filter(!is.na(presence)) %>% mutate( tips = str_replace(Species, "(species in domain Eukaryota)", "")) 
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id")
tr <- multi2di(tr)
# tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
# resolved_names <- readRDS(file = "resolved_names.rds") #from the data cleaning script 
resolved_names <- read.csv(file = "SV_data_avg.csv") #from the data cleaning script

# Make the reference file
# Ensure that the rownames and tip.labels in the target match the species names in the reference

resolved_names$ott_id <- paste("ott", resolved_names$ott_id, sep = "")
resolved_names <- resolved_names %>% rename("genus" = "Genus") %>% rename("family" = "Family") %>% rename("order" = "Order")
resolved_names <- resolved_names %>% rename("unique_name" = "Species")

reference.df <- resolved_names[resolved_names$ott_id %in% tr$tip.label,c("order", "family", "genus", "unique_name", "tips", "ott_id")]
  colnames(reference.df) <- c("order", "family", "genus", "unique_name", "tips_species", "tips")
rownames(reference.df) <- reference.df$tips

# some are duplicated, or have missing data, remove them
reference.df <- reference.df[!duplicated(reference.df$unique_name),]
reference.df <- reference.df[!is.na(reference.df$unique_name),]

# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points
# setwd("/scicore/home/schiera/gizevo30/projects/fish_sleep/")

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/homebrew/bin", sep = ":"))
# test treepl
system("treePL -h")

timetree_order <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_order.nwk")
timetree_family <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_family.nwk")
timetree_genus <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_genus.nwk")

# Use geiger to congruify the tree, works with treePL
# This seems to work up to genus, but not species (by replacing tip.labels with the same names)
geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy
tr.calibrated$tip.label <- resolved_names$tips[match(tr.calibrated$tip.label, resolved_names$ott_id)]

## Save out files

saveRDS(tr.calibrated, file = "tr_tree_calibrated_9_9.rds")

######################################################## 
###################### THIRD PART ######################
######################################################## 




# Below works if you modify the heights.phylo function
trace(geiger:::heights.phylo, edit = TRUE)
depth = max(xx[!(is.na(xx))])
# Also have to do it manually - which is ugh!
timetree_species <- ape::read.tree("timetree_data/actinopterygii_species.nwk")
timetree_species <- multi2di(timetree_species)


resolved_names <- read.csv(file = "SV_data_avg.csv") #from the data cleaning script 
resolved_names$ott_id <- paste("ott", resolved_names$ott_id, sep = "")
resolved_names <- resolved_names %>% rename("genus" = "Genus") %>% rename("family" = "Family") %>% rename("order" = "Order")
resolved_names <- resolved_names %>% rename("unique_name" = "Species")

resolved_names$ott_id <- paste("ott", resolved_names$ott_id, sep = "")
reference.df <- resolved_names[,c("order", "family", "genus", "unique_name", "tips", "ott_id")]
colnames(reference.df) <- c("order", "family", "genus", "unique_name", "tips_species", "tips")
rownames(reference.df) <- reference.df$tips


tr.calibrated <- readRDS("tr_tree_calibrated_9_9.rds")
tr.calibrated$tip.label <- reference.df$tips[match(tr.calibrated$tip.label, reference.df$tips_species)]


geiger.species <- congruify.phylo(reference = timetree_species, target = tr.calibrated, taxonomy = reference.df, tol = 0, scale = "treePL")
tr.calibrated <- geiger.species$phy

tr.calibrated$tip.label <- resolved_names$tips[match(tr.calibrated$tip.label, resolved_names$ott_id)]

# NA value is Diplodus sargus - can tell by plotting the tree
tr.calibrated$tip.label[320] <- "Diplodus_sargus"
tr.calibrated$tip.label <- gsub("\\(species in domain Eukaryota\\)", "", tr.calibrated$tip.label)
saveRDS(tr.calibrated, file = "tr_tree_calibrated_9_9.rds")



