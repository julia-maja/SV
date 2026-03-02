
library(dplyr)
library(ggplot2)
library(RColorBrewer)

rom_cichlids <- readxl::read_excel("/Users/juliamaja/Downloads/ROMI Cichlidae-Nov 2025-to Julia Maja.xls")

cichlids <- rom_cichlids %>%
  mutate(gen_sp = paste(Genus, SpecificEpithet, sep = "_"))

uniq_cichlids <- (cichlids 
                  %>% select(gen_sp, `No of Preserved Specimens`, Country) 
                  %>% group_by(gen_sp) 
                  #%>% summarise("No of Preserved Specimens" = sum(as.numeric(`No of Preserved Specimens`), na.rm = TRUE)) 
                  %>% rename(Species = gen_sp)
)

cichlids_waterbody <- cichlids %>% select(gen_sp, Country, Waterbody) %>% unique()



P50 = c("#483d8b", "#3cb371", "#bc8f8f", "#bdb76b", "#008b8b", "#4682b4", "#d2691e",
        "#F8C471", "#cd5c5c", "#00008b", "#32cd32", "#daa520", "#8fbc8f", "#8b008b",
        "#9932cc", "#ff4500", "#ff8c00", "#ffd700", "#F9E79F", "#244fa4", "#BCBD22",
        "#00fa9a", "#dc143c", "#00ffff", "#00bfff", "#0000ff", "#a020f0", "#adff2f",
        "#ff6347", "#da70d6", "#d8bfd8", "#ff00ff", "#1e90ff", "#db7093", "#dda0dd",
        "#add8e6", "#ff1493", "#7b68ee", "#ffa07a", "#98fb98", "#7fffd4", "#FCCDE5",
        "#ff69b4", "#2f4f4f", "#556b2f", "#a0522d", "#006400",
        "#708090", "#8b0000", "#A6D854"
        )


df <- cichlids_waterbody %>%
  count(Country) %>%
  mutate(
    percent = round(n / sum(n) * 100, 1),
    legend_label = paste0(Country, " — ", n, " (", percent, "%)")
  )

ggplot(df, aes(x = "", y = n, fill = Country)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Occurrences per Country", fill = "Country") +
  scale_fill_manual(
    values = P50,                       
    labels = df$legend_label               
  )

## species

top_species <- uniq_cichlids %>%
  filter(`No of Preserved Specimens` != max(`No of Preserved Specimens`)) %>% 
  arrange(desc(`No of Preserved Specimens`)) %>%
  slice(1:50) 

ggplot(top_species, aes(x = "", y = `No of Preserved Specimens`, fill = gen_sp)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = P50)
  labs(title = "Top 15 Species by Preserved Specimens", fill = "Species")

  uniq_cichlids<-uniq_cichlids %>%
    filter(`No of Preserved Specimens` != max(`No of Preserved Specimens`))

  ggplot(uniq_cichlids, aes(x = `No of Preserved Specimens`)) +
    geom_histogram(binwidth = 50, fill = "steelblue", color = "black") +
    labs(title = "Frequency of Species by Number of Preserved Specimens",
         x = "Number of Preserved Specimens",
         y = "Number of Species") +
    theme_minimal() 
 
  
  library(rfishbase)
   
  SV_data_avg <- read_csv("Desktop/SV/SV_data_avg.csv")
  my_cichlids <- (SV_data_avg %>% filter(Family == "Cichlidae") 
    #%>% filter(genome.assembly == "y") 
    %>% select(Species, SV3, presence) 
    %>% mutate(Species = str_replace(Species, " ", "_")) 
  )
  
  # which ROM species do I not have in my dataset

new_cichlids <- uniq_cichlids %>%
  anti_join(my_cichlids, by = "Species") %>%    # keep only species not in my_cichlids
  group_by(Species) %>%
  summarise(
    Country = first(Country),                   # take the first country
    `No of Preserved Specimens` = sum(as.numeric(`No of Preserved Specimens`), na.rm = TRUE)
  ) %>%
  ungroup() 
# new_cichlids <- as_data_frame(setdiff(uniq_cichlids$gen_sp, my_cichlids$Species))
  
  # which ROM species have genome assemblies
  library(rentrez)
   # Create species_list from uniq_cichlids$gen_sp
  species_list <- uniq_cichlids$gen_sp
  
  # Replace "-" with a space
  species_list <- gsub("-", " ", species_list)
  
  check_assembly <- function(species) {
    # search NCBI Assembly database
    res <- entrez_search(db = "assembly", term = paste0(species, "[ORGN]"))
    # return number of assemblies found
    data.frame(
      species = species,
      n_assemblies = res$count
    )
  }
  
  results <- lapply(species_list, check_assembly)
  results_df <- do.call(rbind, results)
  
  species_with_genomes <- results_df %>%
    filter(n_assemblies > 0) %>%
    pull(species)  # extract the species names as a vector
  
  # Subset your original data frame to just species and country
  result_df <- cichlids %>%
    filter(gen_sp %in% species_with_genomes) %>%
    select(gen_sp, Country, Waterbody) %>% 
    group_by(gen_sp) %>% unique()
  
  ###########
  
  library(stringr)
  
  my_sv_cichlids <- (SV_data_avg  
    %>% filter(Family == "Cichlidae")
    %>% filter(genome.assembly == "y")
    %>% mutate(Species = str_replace(Species, " ", "_"))
    #%>% pull(Species)
  )
   
rom_cichlids_2 <- cichlids %>% 
  filter(gen_sp %in% my_sv_cichlids) %>% 
  select(gen_sp) %>% group_by(gen_sp) %>% unique() 
  %>% mutate(gen_sp = str_replace(gen_sp, " ", "_"))


# presence/ absence of whole set with and without genome assembly
SV_data <- SV_data_avg 
SV_data <- SV_data %>% filter(!is.na(presence)) %>% mutate( tips = str_replace(Species, "(species in domain Eukaryota)", "")) 
tr <- tol_induced_subtree(ott_ids = SV_data$ott_id[SV_data$flags %in% c("sibling_higher", "")], label_format = "id") 
tr <- multi2di(tr)
tr$tip.label <- SV_data$tips[match(tr$tip.label, paste("ott", SV_data$ott_id, sep = ""))]
ggtree(tr, layout = "circular") + geom_tiplab(color = "black", size = 0.8)
sv.plot <- ggtree(tr, layout = "circular") %<+% SV_data[, c("tips", "presence", "genome.assembly", "Order")] 
sv.plot <- sv.plot + geom_tile(data = sv.plot$data[!is.na(sv.plot$data$presence) & seq_len(nrow(sv.plot$data)) <= length(tr$tip.label), ], aes(y=y, x=x, fill = presence), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = pres_abs) 
sv.plot <- sv.plot + new_scale_fill() + geom_tile(data = sv.plot$data[1:length(tr$tip.label),], aes(y=y, x=x + 15, fill = Order), inherit.aes = FALSE, color = "transparent") + scale_fill_manual(values = P60)
sv.plot <- sv.plot + geom_tiplab(aes(color = ifelse(genome.assembly == "y", "highlight", "default")), hjust = -0.2, size = 1) + scale_color_manual(values = c("highlight" = "blue", "default" = "black"))
sv.plot
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  