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

source("/Users/user2/Desktop/SV_project/SV_functions.R")
SV_data_avg <- read.csv("/Users/user2/Desktop/SV_project/SV_data_avg.csv")



# SV presence/ absence
tr <- tol_induced_subtree(ott_ids = SV_data_avg$ott_id[SV_data_avg$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data_avg$tips[match(tr$tip.label, paste("ott", SV_data_avg$ott_id, sep = ""))]

df <- data.frame(names = tr$tip.label, output = SV_data_avg$tips[match(tr$tip.label, paste("ott", SV_data_avg$ott_id, sep = ""))])

ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)

sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data_avg[, c("tips", "presence", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette) 
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