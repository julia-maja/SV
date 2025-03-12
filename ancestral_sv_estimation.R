
trait.vector_n <- SV_data_presence$presence


names(trait.vector_n) <- SV_data_presence$tips
trait.vector_n <- trait.vector_n[trait.vector_n %in% c("present", "absent")]
trait.vector_n <- trait.vector_n[match(tr$tip.label, names(trait.vector_n))]
trait.vector_n <- na.omit(trait.vector_n) 
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



