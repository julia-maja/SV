# Section 0: Import and run packages -------------------------------------------------
# For retreiving data from the open tree of life
library(rotl)
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)
#add timescale to ggtree
library(deeptime)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)



# Section 1: Import data  -----------------------------------------------

#we want this script to be flexible enough to take any model input
#only need to change the model_results and file_name to create plots of new data

#currently I have max_clade_crep data for: artio max_crep, artio max_dinoc
#to do: cetacean max_crep, cetacean max_dinoc, artio w/out cetaceans max_crep, artio w/out cetaceans max_dinoc
all_model_results <- readRDS("SV_reconstruction_results.RDS")
#copy and paste first half of filename here (leave out the models)

#separate the results by the model types we want to use (ER, SYM, ARD, bridge_only)
#uncomment the model you want to plot

# model_results <- all_model_results[[1]]
# model_name <- "ER"

# model_results <- all_model_results[[2]]
# model_name <- "SYM"

model_results <- all_model_results[[3]]
model_name <- "ARD"

# model_results <- Loss_only_model
# model_name <- Loss_only_model

# Section 1: Plotting ancestral reconstruction from corHMM model  --------

lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
colnames(lik.anc) <- c("absent", "present")
phylo_tree <- model_results$phy
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + 
  aes(color = present) + geom_tippoint(aes(color = present), shape = 16, size = 1) + 
  scale_color_distiller(name = "SV presence", palette = "RdYlGn", direction = 1)  + 
  geom_tiplab(color = "black", size = 1.5, offset = 0.5) + 
  geom_tippoint(aes(color = present), shape = 16, size = 1.5)
#ancestral_plot <- ancestral_plot + geom_tiplab(hjust = -0.2, size = 1.5)
ancestral_plot



# # Adding orders and removing species names ------------------------------------

lik.anc$species <- rownames(lik.anc) 

ancestral_plot <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + 
  aes(color = present) + geom_tippoint(aes(color = present), shape = 16, size = 1) + 
  scale_color_distiller(name = "SV presence", palette = "RdYlGn", direction = 1)  + 
  geom_tiplab(color = "black", size = 1.5, offset = 0.5) + 
  geom_tippoint(aes(color = present), shape = 16, size = 1.5)
#ancestral_plot <- ancestral_plot + geom_tiplab(hjust = -0.2, size = 1.5)
ancestral_plot

tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_8_25.rds") # my tree


