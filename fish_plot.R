library(ggtree)
library(dplyr)
library(stringr)
library(ape)
library(ggplot2)

fish_tree <- readRDS("/Users/user2/Desktop/tr_tree_calibrated_fish.rds")

sp_data <- read.csv("/Users/juliamaja/Downloads/saccus_vasculosus.csv") # 243 species
sp_data$species <- str_replace(sp_data$species, " ", "_")
sp_data <- sp_data[sp_data$species %in% fish_tree$tip.label, ] #95

sp_data <- (read.csv("/Users/user2/Desktop/saccus_vasculosus.csv") 
            %>% mutate(species = str_replace(sp_data$species, " ", "_"))
            %>% filter(species %in% fish_tree$tip.label)
)

fish_tree <- keep.tip(fish_tree, tip = sp_data$species)



fish_plot <- ggtree(fish_tree, layout="circular") %<+% sp_data[, c("species", "SV")] 
fish_plot <- fish_plot + geom_tile(data = fish_plot$data[1:length(fish_tree$tip.label), ], aes(x=x, y=y, fill=SV), inherit.aes=FALSE, color="transparent", width = 3 ) + geom_tiplab(size = 2)
fish_plot


fish_plot <- ggtree(fish_tree, layout="circular") %<+% sp_data[, c("species", "SV", "order")] 
fish_plot <- fish_plot + geom_tile(data = fish_plot$data[1:length(fish_tree$tip.label), ], aes(x=x, y=y, fill=SV), inherit.aes=FALSE, color="transparent", width = 3 ) + geom_tiplab(size = 2) 
fish_plot <- fish_plot + geom_tile(data = fish_plot$data[1:length(fish_tree$tip.label), ], aes(x=x+20, y=y, fill=order), inherit.aes=FALSE, color="transparent", width = 3 ) 
fish_plot



