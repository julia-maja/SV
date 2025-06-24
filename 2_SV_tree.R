# Section 0: Packages -----------------------------------------------------
library(ape) 
library(corHMM)
library(phangorn)
library(stringr)
library(here)
library(rotl)
library(ggtree) 
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(ggplot2)
library(rphylopic)
library(RColorBrewer)
library(ggimage)
library(rfishbase)
library(rentrez)
library(taxize)
library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(xlsx) # problem
library(geiger)
library(here)
library(ggtree)
library(randomcoloR)
library(ggnewscale)
library(svglite)
library(Polychrome)



SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

#
# time calibrating --------------------------------------------------------


# # time calibrate the tree:
# 
# reference.df <- resolved_names[resolved_names$ott_id %in% tr$tip.label,c("Order", "Family", "Genus", "Species", "tips", "ott_id")]
# colnames(reference.df) <- c("Order", "Family", "Genus", "Species", "tips_species", "tips")
# rownames(reference.df) <- reference.df$tips
# 
# # # some are duplicated, or have missing data, remove them
# # reference.df <- reference.df[!duplicated(reference.df$Species),]
# # reference.df <- reference.df[!is.na(reference.df$Species),]
# 
# 
# saveRDS(reference.df, "/Users/juliamaja/Desktop/SV/reference_df.rds")
# saveRDS(tr, "/Users/juliamaja/Desktop/SV/tol_induced_tree.rds")
# 
# 
# # Load the timetree tree (genus level data works, but not species)
# # Have download timetree data for species, genus, family, and order
# # Genus level data has the most calibration points
# 
# reference.df <- readRDS("/Users/juliamaja/Desktop/SV/reference_df.rds")
# tr <- readRDS("/Users/juliamaja/Desktop/SV/tol_induced_tree.rds")
# 
# timetree_order <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_order.nwk")
# timetree_family <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_family.nwk")
# timetree_genus <- ape::read.tree("/Users/juliamaja/Desktop/SV/timetree_data/actinopterygii_genus.nwk")
# 
# # Use geiger to congruify the tree, works with treePL
# # This seems to work up to genus, but not species (by replacing tip.labels with the same names)
# 
# geiger.order <- congruify.phylo(reference = timetree_order, target = tr, taxonomy = reference.df, tol = 0, scale = "treePL")
# 
# geiger.family <- congruify.phylo(reference = timetree_family, target = geiger.order$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
# 
# geiger.genus <- congruify.phylo(reference = timetree_genus, target = geiger.family$phy, taxonomy = reference.df, tol = 0, scale = "treePL")
# 
# tr.calibrated <- geiger.genus$phy
# 
# ## Save out files
# 
# saveRDS(tr.calibrated, file = "calibrated_phylo.rds")
# 
# # Add in a stop here
# 
# print("this is the last message")
# stop()
# print("you should not see this")
# 
# cannot do this without treePL


# time-calibration with max's fish tree -----------------------------------



# tr <- readRDS("/Users/juliamaja/Downloads/tr_tree_calibrated_fish.rds")
tr <- readRDS("/Users/juliamaja/Desktop/SV/julia_fish_tree.rds")


SV_data_avg <- SV_data_avg[SV_data_avg$tips %in% tr$tip.label,] #%>% mutate( Species = str_replace(Species, "(species in domain Eukaryota)", "")) %>% mutate(Species = str_replace(Species, "_", "")) %>% mutate(tips = Species)
tr <- keep.tip(tr, tip = SV_data_avg$tips)
ggtree(tr, layout = "circular") + geom_tiplab(size = 1.5)


saveRDS(tr, "/Users/juliamaja/Desktop/SV/fish_time_tree.rds")




# plotting ----------------------------------------------------------------

source("/Users/juliamaja/Desktop/SV/SV_functions.R")
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

resolved_names <- SV_data_avg

# color palettes
rainbow_pal <- colorRampPalette(rainbow(7))(100)
set.seed(12)
random_colors <- distinctColorPalette(99)
P50 = c("#483d8b", "#3cb371", "#bc8f8f", "#bdb76b", "#008b8b", "#4682b4", "#d2691e",
        "#F8C471", "#cd5c5c", "#00008b", "#32cd32", "#daa520", "#8fbc8f", "#8b008b",
        "#9932cc", "#ff4500", "#ff8c00", "#ffd700", "#F9E79F", "#244fa4", "#BCBD22",
        "#00fa9a", "#dc143c", "#00ffff", "#00bfff", "#0000ff", "#a020f0", "#adff2f",
        "#ff6347", "#da70d6", "#d8bfd8", "#ff00ff", "#1e90ff", "#db7093", "#dda0dd",
        "#add8e6", "#ff1493", "#7b68ee", "#ffa07a", "#98fb98", "#7fffd4", "#FCCDE5",
        "#ff69b4", "#2f4f4f", "#556b2f", "#a0522d", "#006400",
        "#708090", "#8b0000", "#A6D854"
)
P60 = c( "#000", "#a9a9a9", "#d3d3d3", "#2f4f4f", "#556b2f", "#6b8e23", "#a0522d",
         "#2e8b57", "#228b22", "#800000", "#191970", "#006400", "#808000", "#483d8b",
         "#b22222", "#5f9ea0", "#778899", "#3cb371", "#bc8f8f", "#663399", "#008080",
         "#bdb76b", "#4682b4", "#d2691e", "#9acd32", "#20b2aa", "#cd5c5c", "#00008b",
         "#4b0082", "#32cd32", "#daa520", "#7f007f", "#8fbc8f", "#b03060", "#d2b48c",
         "#66cdaa", "#9932cc", "#ff0000", "#ff8c00", "#ffa500", "#ffd700", "#ffff00",
         "#c71585", "#0000cd", "#40e0d0", "#7fff00", "#00ff00", "#ba55d3", "#00fa9a",
         "#00ff7f", "#4169e1", "#dc143c", "#00ffff", "#00bfff", "#f4a460", "#9370db",
         "#0000ff", "#a020f0", "#adff2f", "#ff6347", "#d8bfd8", "#b0c4de", "#ff00ff",
         "#1e90ff", "#db7093", "#f0e68c", "#fa8072", "#ffff54", "#6495ed", "#dda0dd",
         "#b0e0e6", "#90ee90", "#ff1493", "#7b68ee", "#ffa07a", "#f5deb3", "#ee82ee",
         "#87cefa", "#7fffd4", "#ff69b4"
)
pres_abs <- c("red", "green")

# SV presence/ absence
SV_data <- SV_data_avg
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
sv_palette <- c("oldlace", "#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "slateblue")
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "presence", "genome.assembly")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_brewer(palette = "RdPu") 
#sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x+15, fill = genome.assembly), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette) 
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot  

# SV presence/ absence + order
SV_data <- SV_data_avg
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
rainbow_pal <- colorRampPalette(rainbow(7))(100)
set.seed(12)
random_colors <- distinctColorPalette(99)
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = random_colors) 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = random_colors)
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)

# Compute Order label positions (insert after all tiles and before labeling)
tip_data <- sv.plot$data[1:length(tr$tip.label), ]
order_labels <- tip_data %>% group_by(Order) %>% summarize(x = mean(x + 20), y = median(y), .groups = "drop")
sv.plot <- sv.plot + geom_text(data = order_labels, aes(x = x, y = y, label = Order), size = 3, hjust = 0)
sv.plot <- sv.plot + theme(legend.position = "none")
sv.plot 


# SV presence/ absence + order  (genome assembly subset)
SV_data <- SV_data_avg %>% filter(genome.assembly == "y") 
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[!is.na(sv.plot$data$presence) & seq_len(nrow(sv.plot$data)) <= length(tr$tip.label), ], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = pres_abs) 
sv.plot <- sv.plot + new_scale_fill() + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = P50)
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot

# presence/ absence of whole set with and without genome assembly
SV_data <- SV_data_avg 
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 0.8)
SV_data <- SV_data_avg
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[!is.na(sv.plot$data$presence) & seq_len(nrow(sv.plot$data)) <= length(tr$tip.label), ], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = pres_abs) 
sv.plot <- sv.plot + new_scale_fill() + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = P60)
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1)
sv.plot

# Compute Order label positions (insert after all tiles and before labeling)
tip_data <- sv.plot$data[1:length(tr$tip.label), ]
order_labels <- tip_data %>% group_by(Order) %>% summarize(x = mean(x + 20), y = median(y), .groups = "drop")
sv.plot <- sv.plot + geom_text(data = order_labels, aes(x = x, y = y, label = Order), size = 3, hjust = 0)
sv.plot <- sv.plot + theme(legend.position = "none")
sv.plot 





# tree building for SV complexity:
sv_palette <- c("black", "transparent")
SV_data <- SV_data_avg %>% filter(!is.na(SV3)) %>% mutate( Species = str_replace(Species, "(species in domain Eukaryota)", "")) %>% mutate(Species = str_replace(Species, "_", "")) %>% mutate(tips = Species)
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 1.5)
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "SV3", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = SV3), inherit.aes = FALSE, color = "transparent") + scale_fill_gradient(name = "SV complexiity", low = "tomato", high = "palegreen", na.value = NA)
#new_scale_color() + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = sv_palette)
sv.plot <- sv.plot + geom_tiplab(hjust = -0.2, size = 1.5)
sv.plot
file_name <- "/Users/juliamaja/Desktop/SV/plots/sv_complexity.svg"
ggsave(file_name, sv.plot, device = "svg")



num_species_total <- sp_data %>% group_by(species) %>% count()
num_species_searched <- sp_data %>% group_by(species) %>% filter(SV != "" | presence != "") %>% count()
num_species_w_data <- sp_data %>% group_by(species) %>% filter((SV != "" & SV != "no data") | presence != "") %>% count()
num_species_w_data_and_genome <- sp_data %>% group_by(species) %>% filter(((SV != "" & SV != "no data") | presence != "") & genome.assembly == "y") %>% count()
num_sp_w_genome <- sp_data %>% group_by(species) %>% filter(genome.assembly == "y")

# add order labels to tree
node_labels <- ggtree(tr, layout = "circular") + geom_text(aes(label=node),colour = "blue", hjust=-.2, size = 2) + geom_tiplab(size = 2, hjust = -0.1)
node_labels
nodes <- findMRCANode(phylo = tr, trait.data = SV_data_avg)

all_nodes <- lapply(unique(SV_data_avg$order), function(x) findMRCANode(phylo = tr, trait.data = SV_data_avg))
nodes_df <- do.call(rbind, all_nodes)
nodes_df <- nodes_df %>% distinct()

#can now easily label all clades within that taxonomic level on the tree using the nodes_df
num_orders <- SV_data_avg %>% group_by(Order) %>% distinct() %>% count()
set.seed(30) 
color_palette <- sample(colors(), 54)  
nodes_df$colour <- color_palette
order_tree <- sv.plot + geom_cladelab(node = nodes_df$node_number,
                                      label = nodes_df$clade_name,
                                      offset = 0.5,
                                      offset.text = 3,
                                      align = TRUE,
                                      fontsize = 2.5,
                                      barsize = 2,
                                      barcolour = nodes_df$colour
)
order_tree


SV_data_avg <- SV_data_avg %>% filter(genome.assembly == "y")
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")

# highlighting my species on the fish tree of life (11637 species)
ref_tree <- read.tree("/Users/juliamaja/Downloads/actinopt_12k_treePL.tre")

sample_species <- gsub(" ", "_", SV_data_avg$Species)
tip_df <- data.frame(label = ref_tree$tip.label) %>%
  mutate(sampled = ifelse(label %in% sample_species, "Sampled", "Other"))

ref.plot <- ggtree(ref_tree, layout = "circular") %<+% tip_df + 
  geom_tippoint(aes(color = sampled), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Sampled" = "red", "Other" = "gray")) +
  theme(legend.position = "right") +
  ggtitle("My dataset highlighted on fish tree of life")
ref.plot

# add orders
orders_df <- read.csv("/Users/juliamaja/Downloads/PFC_taxonomy.csv")
orders_df <- orders_df %>% rename("Species" = "genus.species")
orders_df$Species <- gsub(" ", "_", orders_df$Species)

ref.plot <- ref.plot + geom_tile(data = ref.plot$data[1:length(ref_tree$tip.label),], aes(y=y, x=x + 15, fill = order), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = random_colors)


library(ggtreeExtra)
library(ggnewscale)
ref.plot <- ggtree(ref_tree, layout = "circular") %<+% tip_df +
  geom_tippoint(aes(color = sampled), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Sampled" = "red", "Other" = "gray")) +
  new_scale_fill() +  # allow second fill scale
  
  geom_fruit(
    data = tip_df,
    geom = geom_tile,
    mapping = aes(y = label, fill = order),
    width = 0.04,
    offset = 0.015,
    color = NA
  ) +
  theme(legend.position = "right") +
  ggtitle("Fish Tree with Sampled Tips and Taxonomic Order Ring")

ref.plot


# Taxonomic breadth:
# comparing order breadth with fishbase (35731 species)
# can comment out or in log-scale

ref_counts <- fishbase_species %>%
  count(Order, name = "Ref_Species")
sample_counts <- SV_data_avg %>%
  count(Order, name = "Sample_Species")
merged_counts <- full_join(ref_counts, sample_counts, by = "Order") %>%
  replace_na(list(Ref_Species = 0, Sample_Species = 0)) %>%
  filter(Ref_Species > 0) %>%  # avoid divide by zero
  mutate(Percent = round((Sample_Species / Ref_Species) * 100, 1))
merged_counts <- merged_counts %>%
  arrange(desc(Ref_Species)) %>%
  mutate(Order = factor(Order, levels = unique(Order)))

ggplot(merged_counts, aes(x = Order)) +
  geom_bar(aes(y = Ref_Species), stat = "identity", fill = "gray70") +
  geom_bar(aes(y = Sample_Species), stat = "identity", fill = "red", alpha = 0.8) +
  geom_text(
    aes(y = Sample_Species, label = paste0(Percent, "%")),
    vjust = -0.2, size = 3, color = "black"
  ) +
  #scale_y_log10() +
  labs(
    title = "Species per Order (Fishbase)",
    x = "Order",
    y = "Number of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# my sample size compared to fishbase: 0.01656825
592 / 35731 


# statistical power:
library(ape)
library(geiger)

# calculate lambda and effective sample size:
# Convert "present"/"absent" to 1/0
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)

# Fix species names to match tree tip labels (replace spaces with underscores)
SV_data_avg$Species <- gsub(" ", "_", SV_data_avg$Species)

# Create named numeric vector of trait
trait_vec <- SV_data_avg$presence_num
names(trait_vec) <- SV_data_avg$Species

# Keep only species present in tree
trait_vec <- trait_vec[names(trait_vec) %in% tree$tip.label]

# Reorder trait_vec to match tree tip labels
trait_vec <- trait_vec[tree$tip.label]
# Your tree
tree <- trpy_n

# Estimate phylogenetic signal (Pagel's λ)
lambda_est <- phylosig(tree, x = trait_vec, method = "lambda")

# Estimate effective sample size (ESS)
# High λ = more phylogenetic signal = lower ESS
# Formula: ESS ≈ n / (1 + (n - 1) * λ)
n <- length(tree$tip.label)
lambda <- lambda_est$lambda

ESS <- n / (1 + (n - 1) * lambda)
ESS

