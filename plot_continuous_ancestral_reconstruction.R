
library(phytools)
library(ape)
library(plotrix)
library(dplyr)

#ancestral reconstruction for continuous trait data: SV complexity

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

# filter for genome-only species:
# SV_data_avg <- SV_data_avg %>% filter(genome.assembly == "y")

# filter out species with SV type = NA
SV_data_avg <- SV_data_avg %>% filter(!is.na(SV3))


# trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/fish_time_tree.rds")
trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") # my tree

# removing jawless fish
trpy_n <- drop.tip(trpy_n, c("Eptatretus_atami", "Myxine_glutinosa"))

# filter species with possible false absents
#SV_data_avg <- SV_data_avg %>% filter(is.na(SV.not.in.figure))
#trpy_n <- drop.tip(tr, setdiff(tr$tip.label, SV_data_avg$tips))

# estimate ancestral states
# re-order data to match the tree
SV_data_avg <- SV_data_avg[match(trpy_n$tip.label, SV_data_avg$tips), ]
SV_data_avg <- SV_data_avg %>% filter(!is.na(SV3))

# Create a named vector of SV3 values
SV3 <- SV_data_avg$SV3

# prune the tree to remove species not in my dataset
tr <- drop.tip(trpy_n, setdiff(trpy_n$tip.label, SV_data_avg$tips))

anc_SV3 <- fastAnc(tr, SV3, vars = TRUE, CI = TRUE)


#### plot
# apparently you don't need the results of the ancestral reconstruction 
# to make a plot because contMap does it already

# read in data and make the row names be the tip values
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv", row.names=14)
SV_data_avg <- SV_data_avg %>% mutate(Species = str_replace(Species, " ", "_"))

# re-order data to match the tree
SV_data_avg <- SV_data_avg[match(trpy_n$tip.label, SV_data_avg$Species), ]
#SV_data_avg <- SV_data_avg %>% filter(!is.na(SV3))
# prune the tree to remove species not in my dataset

# filter out species with SV type = NA
SV_data_avg <- SV_data_avg %>% filter(!is.na(SV3))

# trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/fish_time_tree.rds")
trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") # my tree
# removing jawless fish
trpy_n <- drop.tip(trpy_n, c("Eptatretus_atami", "Myxine_glutinosa"))
tr <- drop.tip(trpy_n, setdiff(trpy_n$tip.label, SV_data_avg$Species))

SV_data_avg <- SV_data_avg %>% select(SV3)

svl <- as.matrix(SV_data_avg)[,1]

obj <- contMap(tr, svl, plot = FALSE)

plot(obj, type="fan", lwd=4)











