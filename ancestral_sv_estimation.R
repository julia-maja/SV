

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")


# ace ---------------------------------------------------------------------
trait.vector_n <- SV_data_avg$presence

names(trait.vector_n) <- SV_data_avg$tips
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("present", "absent")]
trait.vector_n <- trait.vector_n[match(tr$tip.label, names(trait.vector_n))]
#trait.vector_n <- na.omit(trait.vector_n) 
trpy_n <- keep.tip(tr, tip = names(trait.vector_n))
trait.vector_n <- trait.vector_n[match(trpy_n$tip.label, names(trait.vector_n))]
 # Removes NA values
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
trait.data_n <- trait.data[trait.data$species %in% trpy_n$tip.label,]

standard_tests <- list()
# Equal rates, symmetric rates (same as ER), and All rates different
standard_tests[[1]] <- ace(trait.vector_n, trpy_n, model = "ER", type = "discrete")
standard_tests[[2]] <- ace(trait.vector_n, trpy_n, model = "SYM", type = "discrete")
standard_tests[[3]] <- ace(trait.vector_n, trpy_n, model = "ARD", type = "discrete")



# corHMM ------------------------------------------------------------------
library("corHMM")

# #430 species
# SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# #SV_data_avg <- SV_data_avg[1:10, ]
# 
# tr <- readRDS("/Users/juliamaja/Desktop/SV/julia_fish_tree.rds")
# 
# 
# #242 species in the calibrated tree
# SV_data_avg <- SV_data_avg[SV_data_avg$tips %in% tr$tip.label, ]
# trpy_n <- keep.tip(tr, tip = SV_data_avg$tips)
# trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
# #trpy_n <- drop.clade(trpy_n, NA)
# 
# ggtree(trpy_n, layout = "circular") + geom_tiplab(size = 1.5)

# #check which should actually be removed
# row.names(SV_data_avg) <- SV_data_avg$tips
# SV_data_avg <-SV_data_avg[!is.na(SV_data_avg$presence), ]
# #242 species
  
trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/fish_time_tree.rds")

ER_model <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 1, model = "ER", node.states = "marginal")
SYM_model <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 1, model = "SYM", node.states = "marginal")
ARD_model <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 1, model = "ARD", node.states = "marginal")
# ARD_model_2HR <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 2, model = "ARD", node.states = "marginal")

#"If analyzing a binary or multistate character, the order of root.p is the same order as the traits – e.g., for states 1, 2, 3, a root.p=c(0,1,0) would fix the root to be in state 2"
#since 1 = echo and 2 = no echo we can use the following to set the state at the root to no echo
ER_set_root_present <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 1, model = "ER", node.states = "marginal", root.p = c(0, 1)) # (absent, present)
ER_set_root_absent <- corHMM(phy = trpy_n, data = SV_data_avg[, c("tips", "presence")], rate.cat = 1, model = "ER", node.states = "marginal", root.p = c(1, 0)) # (absent, present)



model_results_list <- list(ER_model, SYM_model, ARD_model)
saveRDS(model_results_list, "SV_reconstruction_results.RDS")

presence <- SV_data_avg$presence
names(presence) <- SV_data_avg$tips
fitER <- ace(presence, trpy_n,model="ER",type="discrete")


# "total garbage" test:
n = length(tr$tip.label) # size of tree
n_0 = count(SV_data_avg %>% filter(presence == "present")) # num of "presents"
p = n_0/n
lnL_garb = n_0*log(p)+(n-n_0)*log(1-p) # in R, log is the natural logarithm. log base 10 is log10
print(lnL_garb)



