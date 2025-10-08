library(tidyverse)
library(ggtree) 
library(RColorBrewer)


teleost_proteomes <- read.csv("/Users/juliamaja/Downloads/proteomes_taxonomy.csv")
teleost_proteomes$Organism <- str_replace(teleost_proteomes$Organism, " ", "_")
teleost_proteomes$Organism <- str_replace_all(teleost_proteomes$Organism, "\\s*\\(.*?\\)", "")

SV_data_avg <- readRDS("/Users/juliamaja/Desktop/SV/SV_data_avg.rds")

setdiff(teleost_proteomes$Organism, SV_data_avg$tips)
intersection <- as_data_frame(intersect(SV_data_avg$tips, teleost_proteomes$Organism))

sp_w_proteomes <- SV_data_avg[SV_data_avg$tips %in% intersection$value, ]


# tree building for SV complexity + orders:
sv_palette <- c("black", "transparent")
SV_data <- sp_w_proteomes %>% filter(!is.na(SV3)) %>% mutate( Species = str_replace(Species, "(species in domain Eukaryota)", "")) %>% mutate(Species = str_replace(Species, "_", "")) %>% mutate(tips = Species)
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV3", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV3), inherit.aes = FALSE, color = "transparent") + scale_fill_gradient(name = "SV complexiity", low = "tomato", high = "palegreen", na.value = NA)
sv.plot <- sv.plot + new_scale_fill() + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_discrete("okabe")
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 2)
sv.plot

# tree plot but time calibrated:
trpy_n <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_8_25.rds")# my tree
trpy_n <- drop.tip(trpy_n, setdiff(trpy_n$tip.label, SV_data$tips))
SV_data <- sp_w_proteomes %>% filter(!is.na(SV3)) %>% mutate( Species = str_replace(Species, "(species in domain Eukaryota)", "")) %>% mutate(Species = str_replace(Species, "_", "")) %>% mutate(tips = Species)
ggtree(trpy_n, layout = "circular") + geom_tiplab(color = "black", size = 2) 
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV3", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(trpy_n$tip.label),], aes(y=y, x=x, fill = SV3), inherit.aes = FALSE, color = "transparent") + scale_fill_gradient(name = "SV complexiity", low = "tomato", high = "palegreen", na.value = NA)
sv.plot <- sv.plot + new_scale_fill() + geom_tile(data = sv.plot$data[1:length(trpy_n$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_discrete("okabe")
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 2)
sv.plot
