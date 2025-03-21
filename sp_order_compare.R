sp_order <- read.csv("/Users/juliamaja/Desktop/SV/Species_order.csv")

sp_order_compare <- sp_order %>% left_join(fishbase_species, by = "Species") %>% select("Species", "order", "Order") %>% filter(order != Order)
