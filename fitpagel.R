
library(phytools)

# load in dataset and tree
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") 
# note: tr_tree_calibrated_9_8 is built without 'false absent' species but SV_data_avg has them
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
# tr_tree_calibrated_9_9.rds has all species.

# model function requires binary data (0/1)
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
SV_data_avg <- SV_data_avg %>% filter(salinity == "freshwater" | salinity == "saltwater")
SV_data_avg$salinity_num <- ifelse(SV_data_avg$salinity == "saltwater", 1, 0)

# it also requires the row names to be species names so it can match the data to the tree
rownames(SV_data_avg) <- SV_data_avg$tips

data <- SV_data_avg
trait1 <- as.integer(data$presence_num)
trait2 <- as.integer(data$salinity_num)

tr <- drop.tip(tr, setdiff(tr$tip.label, SV_data_avg$tips))
SV_data_avg <- SV_data_avg %>% filter(tips %in% tr$tip.label)

fit <- fitPagel(tr, trait1, trait2, model = "ARD")  # or "ER", "SYM"
print(fit)

# Error in x[tree$tip.label, ] : subscript out of bounds
# matching the tree and data species does not fix this

# why:
# > setdiff(tr$tip.label, SV_data_avg$tips)
# character(0)
# > setdiff(SV_data_avg$tips, tr$tip.label)
# character(0)
