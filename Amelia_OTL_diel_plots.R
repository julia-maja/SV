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
#install.packages("rphylopic")
library(rphylopic)
library(RColorBrewer)
#install.packages("ggimage")
library(ggimage)

setwd(here())


# Section 1: Load in and clean up cetacean data ---------------------------

#load in the cetacean diel behaviour dataframe
sp_data <- read.csv("/Users/user2/Desktop/SV.csv")

#in the open tree of life sperm whales are called Physeter catodon so change the name
#cetaceans_full$tips <- str_replace(cetaceans_full$tips, "Physeter_macrocephalus", "Physeter_catodon")


# Section 2: Subset the open tree of life ----------------------------------------

#create trait.data dataframe to use to build tree and add trait data (in geom tile)
#check that all the entries have a value for diel pattern (use diel pattern 3 instead of 1 because di/noc has NA for all cath sps)
#for trait data only needs the diel patterns, species names, tips (unless we want to add more data to the tree like confidence, method etc)
# trait.data <- cetaceans_full[cetaceans_full$Diel_Pattern_3 %in% c("diurnal", "nocturnal", "cathemeral", "crepuscular"),]
#for trait data only needs the diel patterns, species names, tips (unless we want to add more data to the tree like confidence etc)
# trait.data <- trait.data[, c("Species_name", "Diel_Pattern_1", "Diel_Pattern_2", "Diel_Pattern_3", "tips")]
#trait data has all the behavioural data for the 80 species with data (20 species have no behavioural data currently)

trait.data <- sp_data[,c("species", "SV", "order")]
trait.data$tips <- trait.data$species
trait.data$tips <- str_replace(trait.data$species, pattern = " ", replacement = "_")

#find the species names that are in the open tree of life
#looks for species in the same format as Species_names
resolved_names <- tnrs_match_names(names = trait.data$species, context_name = "Vertebrates", do_approximate_matching = TRUE)

## Remove any that don't have exact matches or are synonyms (this is really just for finding ancestral state, so missing a few species won't matter)
#Remove species without open tree of life entries. Removes Balaenoptera riceii (newly discovered, not in the tree)
resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
#Removes species that are synonyms for other species. Removes Inia humboldtiana (is a synonym/subspecies for Inia geoffrensis) and Neophocaena sunameri (subspecies of Neophocaena asiaeorientalis)
# resolved_names <- resolved_names[resolved_names$is_synonym == FALSE,]
#no approximate matches, so this removes no species currently
# resolved_names <- resolved_names[resolved_names$approximate_match == FALSE,]

#we've removed these three species, so now we have 77 cetaceans in our resolved_names dataframe to work with

# Remove excess information, clean up, and add tip label ids that will match the tree
resolved_names <- resolved_names[,c("search_string", "unique_name", "ott_id", "flags")]
resolved_names$tips <- str_replace(resolved_names$unique_name, " ", "_")
resolved_names <- resolved_names[!duplicated(resolved_names$tips),]

# Add data on activity patterns, starting with diel pattern 3 (crep, cath, di, noc)
# match() finds the position of its first argument in its second
resolved_names$SV <- trait.data$SV[match(resolved_names$search_string, tolower(trait.data$species))] # maybe trait.data$tips?
#check that the activity patterns look correct
table(resolved_names$SV)


## Fetch the tree
#subset the open tree of life to only include the species in resolved_names (find them by their ott_id)
tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") 
# I need to use the id option here, and then use that to map the tip labels from resolved_names (that way I don't run into the issue with the difference in formatting between the two tools)
## check flags and see which ones to remove

# Time calibrate it using geiger and timetree.org - the open tree of life is already time calibrated, so we don't need to do this step
# First resolve polytomies ~randomly using multi2dr

tr <- multi2di(tr)

#currently the tips are identified by their ott_id, can change to their species name (tip column)
tr$tip.label <- resolved_names$tips[match(tr$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]

#this is the tree of all the cetacean species we are working with in the open tree of life
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 2)

#subset the tree to make a tree with only the strictly diurnal and nocturnal species
#create drop tips with NA values leftover from the strictly diurnal/nocturnal column 

# trait.data$species <- str_replace[trait.data, " ", "Unknown"]
# to_drop <- trait.data[is.na(trait.data$SV), ]
#check that these are the species we want to drop (should be cathemeral species)
# to_drop$tips
#use functions drop tip to remove these species
# tr2 <- drop.tip(tr, to_drop$tips)

#this is the tree of all the non-cathemeral species
# ggtree(tr2, layout = "circular") + geom_tiplab(color = "black", size = 2)


# Section 2: Creating the diel plots w OTL --------------------------------

#pick colours! Using a custom colour palette for consistency between the plots
display.brewer.all(type = "qual", colorblindFriendly = TRUE)
display.brewer.all(type = "seq")
display.brewer.pal(8, "Set2")
custom.colours1 <- c("#FC8D62","#66C2A5")

sv.plot <- ggtree(tr, layout = "circular") %<+% resolved_names[, c("tips", "SV")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "YlOrRd") 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2)
sv.plot 

node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node),colour = "blue", hjust=-.2, size = 2) + geom_tiplab(size = 2, hjust = -0.1)
node_labels


# labelling families/ orders ----------------------------------------------
# find mrca of a set of species

trait.data <- trait.data %>% filter(order != "")
  
  findMRCANode2 <- function(phylo = tr, trait.data = trait.data, taxonomic_level_col = 3, taxonomic_level_name = "Araneae"){
    nodes_list <- list()
    trait.data.filtered <- trait.data[trait.data[,taxonomic_level_col] == taxonomic_level_name, ]
    tip_vector <- as.vector(trait.data.filtered[, "tips"])
    #find the node number of the MRCA of these species
    MRCA_node <- findMRCA(phylo, tip_vector$tips)
    node_df <- data.frame(clade_name = taxonomic_level_name, node_number = MRCA_node)
    return(node_df)
  }
  
  all_nodes <- lapply(unique(trait.data$order), function(x) findMRCANode2(phylo = tr, trait.data = trait.data, taxonomic_level_col  =3, taxonomic_level_name = x))
  nodes_df <- do.call(rbind, all_nodes)
  
  #can now easily label all clades within that taxonomic level on the tree using the nodes_df
  nodes_df$colour <- c("navy", "slateblue", "mediumpurple", "dodgerblue", "darkorchid1", "royalblue", "lightslateblue", "purple3", "steelblue", "red", "yellow", "blue", "green", "orange", "brown", "turquoise")
  order_tree <- sv.plot + geom_cladelab(node = nodes_df$node_number, label = nodes_df$clade_name, offset = 1.8, offset.text = 3, align = TRUE, fontsize = 2.5, barsize = 2.5, barcolour = nodes_df$colour)
  order_tree



# amelia's stuff ----------------------------------------------------------
# ##first we'll make a plot of the binary dataset (nocturnal, diurnal)
# #plot the subset tree (tr2) and add on the trait data for diel pattern 1
# diel.plot1 <- ggtree(tr2, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_1")]
# #create a layer to display the trait data (geom_tile)
# diel.plot1 <- diel.plot1 + geom_tile(data = diel.plot1$data[1:length(tr2$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_1), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = Set2)
# #add tip labels
# diel.plot1.labelled <- diel.plot1 + geom_tiplab(color = "black", size = 1)
# diel.plot1.labelled
# 
# #export plot as a png
# png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel1.png", width=11,height=8,units="cm",res=1200)
# print(diel.plot1.labelled)
# dev.off()
# 
# ##now we can repeat the process for diel pattern 2 (cath, di/crep, noc/crep, di, noc), using the non-subset tree
# custom.colours2 <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854")
# diel.plot2 <- ggtree(tr, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_2")]
# diel.plot2 <- diel.plot2 + geom_tile(data = diel.plot2$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours2)
# diel.plot2.labelled <- diel.plot2 + geom_tiplab(color = "black", size = 1)
# diel.plot2.labelled
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel2.png", width=14,height=10,units="cm",res=1200)
# print(diel.plot2.labelled)
# dev.off()
# 
# ## repeat process again for diel pattern 3 (cath, noc, di, crep)
# custom.colours3 <- c("#dd8ae7", "#d5cab2", "#FC8D62", "#66C2A5")
# custodiel.plot3 <- ggtree(tr, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_3")]
# diel.plot3 <- diel.plot3 + geom_tile(data = diel.plot3$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_3), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours3)
# diel.plot3.labelled <- diel.plot3 + geom_tiplab(color = "black", size = 1)
# diel.plot3.labelled
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel3.png", width=14,height=10,units="cm",res=1200)
# print(diel.plot3.labelled)
# dev.off()
# 
# 
# # Section 3: Additional diel plots ----------------------------------------
# 
# ##trinary dataset plot
# #create a dataframe of diel pattern dataset with three states (cathemeral, nocturnal, diurnal)
# trait.data3state <- cetaceans_full[!(is.na(cetaceans$Diel_Pattern_2)),]
# trait.data3state <- trait.data3state[,c("Species_name", "Diel_Pattern_2", "tips")]
# trait.data3state$Diel_Pattern_2 <- str_replace_all(trait.data3state$Diel_Pattern_2, "nocturnal/crepuscular", "nocturnal")
# trait.data3state$Diel_Pattern_2 <- str_replace_all(trait.data3state$Diel_Pattern_2, "diurnal/crepuscular", "nocturnal")
# 
# #plot out the three state tree
# custom.colours4 <- c("#dd8ae7", "#FC8D62", "#66C2A5")
# diel.plot3_state <- ggtree(tr, layout = "circular") %<+% trait.data3state[,c("tips", "Diel_Pattern_2")]
# diel.plot3_state <- diel.plot3_state + geom_tile(data = diel.plot3_state$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = Diel_Pattern_2), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = custom.colours4)
# diel.plot3_state <- diel.plot3_state + geom_tiplab(color = "black", size = 1.1)
# diel.plot3_state
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_diel_trinary.png", width=15,height=12,units="cm",res=1200)
# print(diel.plot3_state)
# dev.off()
# 
# ##diel plot with clade labels for each major family
# #find the nodes for each family to add clade labels
# node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node, colour = "red"), hjust=-.2, size = 1.8) + geom_tiplab(size = 1.8, hjust = -0.1)
# node_labels 
# 
# #create a dataframe with each cetacean family and its node
# family_names_right <- data.frame(names = c("Eschrichitiidae", "Neobalaenidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Platanistidae", "Kogiidae", "Physteridae"), nodes = c(72, 73, 151, 143, 88, 81, 139, 60))
# family_names_left <- data.frame(names = c("Ziphiidae", "Lipotidae", "Phocoenidae", "Iniidae", "Monodontidae"), nodes = c(128, 3, 123, 85, 127))
# #add the colour you want each family to be
# family_names_right$colour <- c("slateblue3", "royalblue4", "darkslateblue", "slateblue1", "cyan1", "dodgerblue3","royalblue3", "blueviolet" )
# family_names_left$colour <- c("royalblue1", "#1899f5", "#19e2ff","deepskyblue", "turquoise2")
# 
# diel_families <- diel.plot2 + geom_tiplab(color = "black", size = 1) +
# geom_cladelab(family_names_right, node = family_names_right$nodes, label = family_names_right$names, offset=11, align=FALSE, angle=12, offset.text=1, barsize=1, fontsize=2.5, textcolor=family_names_right$colour, barcolor=family_names_right$colour, alpha = 0.75) +
# geom_cladelab(family_names_left, node = family_names_left$nodes, label = family_names_left$names, offset=11, align=FALSE, angle=192, offset.text=1, barsize=1, fontsize=2.5, textcolor= family_names_left$colour, barcolor= family_names_left$colour, alpha = 0.75)
# #does this actually make the plot smaller?
# diel_families <- diel_families + theme(plot.margin = unit(c(14,8,14,8), "mm"))
# #geom_cladelab doesn't give a bar to clades with single sps so add them manually
# diel_families <- diel_families + geom_strip(60, 61, barsize=1, color='blueviolet', offset = 11)
# diel_families <- diel_families + geom_strip(73, 72, barsize=1, color='royalblue4', offset = 11)
# diel_families <- diel_families + geom_strip(63, 72, barsize=1, color='slateblue3', offset = 11)
# diel2_families <- diel_families + geom_strip(3, 6, barsize=1, color='#1899f5', offset = 11)
# diel2_families
#                           
# png("C:/Users/ameli/OneDrive/Documents/R_projects/diel2_with_families.png", width=15,height=12,units="cm",res=1200)
# print(diel2_families)
# dev.off()
# 
# #better family plot
# #create a dataframe with each cetacean family and its node
# family_names_right <- data.frame(names = c("Eschrichitiidae", "Neobalaenidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Platanistidae", "Kogiidae", "Physteridae"), nodes = c(72, 73, 151, 143, 88, 81, 139, 60))
# family_names_left <- data.frame(names = c("Ziphiidae", "Lipotidae", "Phocoenidae", "Iniidae", "Monodontidae"), nodes = c(128, 3, 123, 85, 127))
# #add the colour you want each family to be
# family_names_right$colour <- c("steelblue2", "steelblue3", "steelblue4", "steelblue1", "grey1", "grey50","grey60", "grey65" )
# family_names_left$colour <- c("grey45", "grey40", "grey15","grey30", "grey20")
# 
# #odontocetes in grey, mysticetes in blue
# family_names_right$parvorder <- "grey40"
# family_names_right[c(1,2,3,4), "parvorder"] <- "steelblue2"
# family_names_left$parvorder <- "grey40"
# 
# diel_families <- diel.plot3 +
#   geom_cladelab(family_names_right, node = family_names_right$nodes, label = family_names_right$names, offset=1, align=FALSE, angle=16, offset.text=1, barsize=1, fontsize=3, textcolor= family_names_right$colour, barcolor=family_names_right$parvorder, alpha = 0.75) +
#   geom_cladelab(family_names_left, node = family_names_left$nodes, label = family_names_left$names, offset=1, align=FALSE, angle=196, offset.text=1, barsize=1, fontsize=3, textcolor= family_names_left$colour, barcolor= family_names_left$parvorder, alpha = 0.75)
# 
# #geom_cladelab doesn't give a bar to clades with single sps so add them manually
# diel_families <- diel_families + geom_strip(60, 61, barsize=1, color='grey50', offset = 1, align=FALSE)
# diel_families <- diel_families + geom_strip(73, 72, barsize=1, color='steelblue1', offset = 1)
# diel_families <- diel_families + geom_strip(63, 72, barsize=1, color='steelblue1', offset = 1)
# diel3_families <- diel_families + geom_strip(3, 6, barsize=1, color='grey50', offset = 1)
# diel3_families
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/diel3_with_families.png", width=20,height=15,units="cm",res=1200)
# print(diel3_families)
# dev.off()
# 
# ##create a diel plot with phylopic images for major families 
# 
# #find images from phylopics collection, not necessary to rerun once you have the image info
# #keep commented out unless you want to find different phylopics ot use
# # img <- pick_phylopic(name = "Mysticeti", n = 19, view = 19)
# # img <- pick_phylopic(name = "Delphinida", n = 54, view = 54)
# # img <- pick_phylopic(name = "Platanistoidea", n = 3, view = 3)
# # img <- pick_phylopic(name = "pan-Physeteroidea", n = 10, view = 10)
# # img <- pick_phylopic(name = "Ziphioidea", n = 20, view = 20)
# # img <- pick_phylopic(name = "Lipotidae", n = 1, view = 1)
# # img <- pick_phylopic(name = "Phocoena", n = 4, view = 4)
# # img <- pick_phylopic(name = "Inia", n = 1, view = 1)
# # img <- pick_phylopic(name = "Monodon", n = 2, view = 2)
# 
# 
# #add image names into the dataframe
# family_names_right$images <- c("8b73f54f-15e8-41b8-8c9c-46c86a185104", "941169f7-bc86-4030-bb28-4419b78214de", "fdff3c1b-dd0a-44d5-9fd3-2ecda7939846", "012afb33-55c3-4fc6-9ae3-3a91fda32fd5", "3caf4fbd-ca3a-48b4-925a-50fbe9acd887",  "a55581b9-72c9-4ede-8fba-b908b08d94c9", "5bfb840e-071f-4a1a-b101-0a747a5453e7", "dc76cbdb-dba5-4d8f-8cf3-809515c30dbd")
# family_names_left$images <- c("2fad47bd-9a34-4bc1-b6cc-b7c4415b109d", "103dbb51-1bce-4df6-b1ac-ff0d38d188a3", "970f7e8d-a823-45b9-85b9-767704d0c13f", "f36c9daa-a102-42dd-88ac-a126753943d2", "bfe45de4-d8a8-423c-abf3-087a7c7d0d6c")
#   
# #map these onto the phylogeny
# diel3_fam_pics <- diel.plot3 + geom_cladelab(data = family_names_right, 
#                         mapping = aes(node = nodes, label = names, image = images), 
#                         geom = "phylopic", imagecolor = "grey30", 
#                         offset=0.5, offset.text=1, alpha = 1)
# 
# diel3_fam_pics <- diel3_fam_pics + geom_cladelab(data = family_names_left, 
#                                                  mapping = aes(node = nodes, label = names, image = images), 
#                                                  geom = "phylopic", imagecolor = "grey30", 
#                                                  offset=0.5, offset.text=1, alpha = 1)
# 
# diel3_fam_pics <- diel3_fam_pics +
#   geom_cladelab(family_names_right, node = family_names_right$nodes, label = family_names_right$names, offset=3, align=FALSE, angle=16, offset.text=1, barsize=0, fontsize=3, textcolor= "black", barcolor="white", alpha = 1) +
#   geom_cladelab(family_names_left, node = family_names_left$nodes, label = family_names_left$names, offset=3, align=FALSE, angle=196, offset.text=1, barsize=0, fontsize=3, textcolor= "black", barcolor= "white", alpha = 1)
# 
# diel3_fam_pics
# 
# png("C:/Users/ameli/OneDrive/Documents/R_projects/diel3_with_phylopics.png", width=20,height=15,units="cm",res=1200)
# print(diel3_fam_pics)
# dev.off()


