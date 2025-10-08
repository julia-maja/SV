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

# Set the working directory and source the functions (not used yet)
setwd(here())

#with this script we will add cetaceans to the rest of artiodactyla and do ancestral trait reconstruction
#as well as modelling the evolution of diel patterns, to see if it's similar or different
#to the results seen with only cetaceans

source("scripts/fish_sleep_functions.R")
source("scripts/Amelia_functions.R")


# Section 1: Import data  -----------------------------------------------

#we want this script to be flexible enough to take any model input
#only need to change the model_results and file_name to create plots of new data

#currently I have max_clade_crep data for: artio max_crep, artio max_dinoc
#to do: cetacean max_crep, cetacean max_dinoc, artio w/out cetaceans max_crep, artio w/out cetaceans max_dinoc
all_model_results <- readRDS("SV_reconstruction_results.RDS")
#copy and paste first half of filename here (leave out the models)
file_name <- "SV_reconstruction_results_plot"

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

#from the model results file, tip states describes the trait states at the tips, states describes the trait states at the nodes
lik.anc <- as.data.frame(rbind(model_results$tip.states, model_results$states))
#for max_crep cath/crep makes more sense, for max_dinoc cathemeral makes more sense
colnames(lik.anc) <- c("absent", "present")
#colnames(lik.anc) <- c("cathemeral", "diurnal", "nocturnal")
phylo_tree <- model_results$phy
#associate each of these species and their trait states with its node
lik.anc$node <- c(1:length(phylo_tree$tip.label), (length(phylo_tree$tip.label) + 1):(phylo_tree$Nnode + length(phylo_tree$tip.label)))

#plot the ancestral reconstruction, displaying each of the three trait states (cathemeral, diurnal, nocturnal)
ancestral_plot <- ggtree(phylo_tree, layout = "circular") %<+% lik.anc + aes(color = present) + geom_tippoint(aes(color = present), shape = 16, size = 1) + scale_color_distiller(name = "SV presence", palette = "RdYlGn", direction = 1)  + geom_tiplab(color = "black", size = 1.5, offset = 0.5) + geom_tippoint(aes(color = present), shape = 16, size = 1.5)
#ancestral_plot <- ancestral_plot + geom_tiplab(hjust = -0.2, size = 1.5)
ancestral_plot


ggsave("/Users/juliamaja/Desktop/SV/ARD_reconstruction.svg", ancestral_plot, device = "svg")

#create the name of the file by pasting together ancestral recon, the diel state and the file_name 
pdf(paste("/Users/juliamaja/Desktop/SV/plots", "ARD_SV", "_", model_name, ".pdf", sep = ""))
ancestral_plot
dev.off()



# # Pie chart ancestral reconstruction ------------------------------------

#load in ARD model data
all_model_results <- readRDS("SV_reconstruction_results.RDS")
model_results <- all_model_results[[3]]
file_name <- "SV_reconstruction_results_pie_plot"

phylo_tree <- model_results$phy

#rename column names for consistency in the next steps
colnames(model_results$data) <- c("tips", "presence")

#to make more clear we can colour the tips separately using geom_tipppoint 
#may have to adjust what trait data column is called in each
base_tree <- ggtree(phylo_tree, layout = "rectangular") + geom_tiplab(size = 2, hjust = -0.1)
base_tree <- base_tree %<+% model_results$data[, c("tips", "presence")]
base_tree <- base_tree + geom_tippoint(aes(color = presence), size = 3) 
base_tree

#make the dataframe of likelihoods at the internal nodes without the tips
lik.anc <- as.data.frame(model_results$states)

#for cetaceans we have to add 72 because we are skipping the tips (nodes 1-72)
#the internal nodes start at 73 and end at node 143
#for artiodactyla we add 300 because we are skipping the tips (nodes 1-300)
#the internal nodes start at 301 and end at node 599
lik.anc$node <- c(1:nrow(lik.anc)) + nrow(model_results$data)

#get the pie charts from this database using nodepie
#the number of columns changes depending on how many trait states
pie <- nodepie(lik.anc, 1:(length(lik.anc)-1))

pie_tree <- base_tree + geom_inset(pie, width = .01, height = .01) 
#this adds a the timescale for the entire tree
pie_tree <- pie_tree + theme_tree2()
#reverses the timescale so it starts at 0mya at the tips and extends back to 50mya at ancestor
pie_tree <- revts(pie_tree)
pie_tree

#save out
png(paste("C:/Users/ameli/OneDrive/Documents/R_projects/New_ancestral_recon/pie_chart/", "pie_chart_anc_recon_", file_name, ".png", sep = ""), width=40,height=20, units="cm",res=600)
pie_tree
dev.off()

