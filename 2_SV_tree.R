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
library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(xlsx)
library(geiger)
library(here)

source("/Users/juliamaja/Desktop/SV/SV_functions.R")
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

resolved_names <- SV_data_avg

# SV presence/ absence
tr <- tol_induced_subtree(ott_ids = resolved_names$ott_id[resolved_names$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
resolved_names$ott_id <- paste("ott", resolved_names$ott_id, sep = "")



# time calibrate the tree:

reference.df <- resolved_names[resolved_names$ott_id %in% tr$tip.label,c("Order", "Family", "Genus", "Species", "tips", "ott_id")]
colnames(reference.df) <- c("Order", "Family", "Genus", "Species", "tips_species", "tips")
rownames(reference.df) <- reference.df$tips

# # some are duplicated, or have missing data, remove them
# reference.df <- reference.df[!duplicated(reference.df$Species),]
# reference.df <- reference.df[!is.na(reference.df$Species),]


saveRDS(reference.df, "/Users/juliamaja/Desktop/SV/reference_df.rds")
saveRDS(tr, "/Users/juliamaja/Desktop/SV/tol_induced_tree.rds")


# Load the timetree tree (genus level data works, but not species)
# Have download timetree data for species, genus, family, and order
# Genus level data has the most calibration points

reference.df <- readRDS("/Users/juliamaja/Desktop/SV/reference_df.rds")
tr <- readRDS("/Users/juliamaja/Desktop/SV/tol_induced_tree.rds")

timetree_order <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_order.nwk")
timetree_family <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_family.nwk")
timetree_genus <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_genus.nwk")

# Use geiger to congruify the tree, works with treePL
# This seems to work up to genus, but not species (by replacing tip.labels with the same names)

geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")

geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")

tr.calibrated <- geiger.genus$phy

## Save out files

saveRDS(tr.calibrated, file = "calibrated_phylo.rds")

# Add in a stop here

print("this is the last message")
stop()
print("you should not see this")












# tree for sv absence/ presence:
tr$tip.label <- resolved_names$tips[match(tr$tip.label, paste("ott", resolved_names$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)

sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data_avg[, c("tips", "presence", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette) 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot 

# tree building for SV complexity:
SV_data <- SV_data_avg
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV3", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV3), inherit.aes = FALSE, color = "transparent") + scale_fill_gradient(low = "tomato", high = "palegreen", na.value = NA)
#sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_color_gradient(low = "white", high = "black")
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot



# add order labels to tree
node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node),colour = "blue", hjust=-.2, size = 2) + geom_tiplab(size = 2, hjust = -0.1)
node_labels
nodes <- findMRCANode(phylo = tr, trait.data = SV_data_orders)

all_nodes <- lapply(unique(SV_data_avg$order), function(x) findMRCANode(phylo = tr, trait.data = SV_data_avg))
nodes_df <- do.call(rbind, all_nodes)
nodes_df <- nodes_df %>% distinct()

#can now easily label all clades within that taxonomic level on the tree using the nodes_df
num_orders <- SV_data_avg %>% group_by(order) %>% distinct() %>% count()
set.seed(30) 
color_palette <- sample(colors(), 54)  
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