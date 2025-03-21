library("corHMM")

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

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
#SV_data_avg <- SV_data_avg[1:10, ]

SV_data_avg <- SV_data_avg[SV_data_avg$tips %in% tr$tip.label, ]
trpy_n <- keep.tip(tr, tip = SV_data_avg$tips)
trpy_n$edge.length[trpy_n$edge.length == 0] <- 0.001 
#trpy_n <- drop.clade(trpy_n, NA)

#check which should actually be removed

row.names(SV_data_avg) <- SV_data_avg$tips

SV_ER_reconstruction <- corHMM(phy = trpy_n, data = SV_data_avg[!is.na(trpy_n$tip.label), !is.na(c("presence", "tips"))], rate.cat = 1, model = "ER", node.states = "marginal")

presence <- SV_data_avg$presence
names(presence) <- SV_data_avg$tips
fitER <- ace(presence, trpy_n,model="ER",type="discrete")

