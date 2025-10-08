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

# load google sheet
url <- 'https://docs.google.com/spreadsheets/d/1Z2P6dAcoU0-Kh0UqVgIAeongX5IJohuMUKG1PyrNraU/edit?gid=84578656#gid=84578656'
sp_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) # as of april 17: 667 sp with data, 1183 left to check. the rest are "no data".
trait.data <- sp_data[,c("species", "SV", "order")]
trait.data$tips <- trait.data$species
#trait.data$tips <- str_replace(trait.data$species, pattern = " ", replacement = "_") 
trait.data <- trait.data %>% mutate(species = str_replace(species, " ", "_"))

# resolve spcies names with open tree of life
otol_names <- tnrs_match_names(names = trait.data$species, context_name = "Vertebrates", do_approximate_matching = TRUE) 
resolved_names <- (otol_names[!(is.na(otol_names$unique_name)),] 
                   %>%  rename( "species" = "search_string") 
                   %>% mutate(species = str_replace(species, " ", "_"))
                   %>% dplyr::select(species, unique_name, "ott_id", "flags")
)
#resolved_names <- distinct(resolved_names, ott_id, .keep_all = TRUE)
saveRDS(resolved_names, "resolved_names.rds")
resolved_names <- readRDS("resolved_names.rds")

SV_data <- sp_data %>% mutate(species = str_replace(species, " ", "_")) # as of april 17: 597 obs. as of june 6: 773 obs as of june 19, 2456 obs
SV_data <- SV_data %>% mutate(species = str_replace(species, " ", "_"))
SV_data <- SV_data %>% mutate(species = tolower(species))
SV_data <- SV_data %>% left_join(resolved_names, by = "species")  
SV_data <- SV_data %>% mutate(tips = unique_name)

SV_data <- SV_data %>% group_by(unique_name) 
SV_data <- SV_data %>% mutate(
  genome.assembly = if_else(
    genome.assembly == '' & any(genome.assembly == 'y'),  
    'y',  
    genome.assembly 
  )
) 
SV_data <- SV_data %>% ungroup()

SV_data <- SV_data %>% rename("Source" = "source")
#SV_data <- SV_data %>% filter(Source != "") 
SV_data <- SV_data %>% filter(SV != "no data") 
SV_data <- SV_data %>% filter(SV != "m") 
SV_data <- SV_data %>% mutate(SV = str_replace(SV, "2", "1"), str_replace(SV, "0.5", "1"), str_replace(SV, "2.5", "2"), str_replace(SV, "3.5", "4")) 

SV_data_3 <- SV_data %>% filter(is.na(unique_name)) %>% filter(!grepl("sp.", species)) %>% filter(species !="")   # to see which species names are not in OTOL # 11 cases as of july 7th
SV_data_2 <- SV_data %>% filter(unique_name !="")
SV_data_4 <- anti_join(SV_data, SV_data_2)

SV_data <- SV_data %>% filter(unique_name !="")
SV_data <- SV_data %>% dplyr::select("species", "unique_name", "order", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure", "image.files", "description", "Source")

SV_data <- SV_data %>% 
  filter(!(is.na(SV) | SV == "") | !(is.na(presence) | presence == ""))
#SV_data <- SV_data %>% filter(SV != "" & presence != "")
SV_data <- SV_data %>% mutate(presence = ifelse(SV == "0", "absent", "present"))
SV_data <- SV_data %>% mutate(SV = ifelse(SV == "present", NA, SV))
SV_data <- SV_data %>% mutate(SV = ifelse(SV == "present ", NA, SV))

SV_data <- SV_data %>% rename("Species" = "unique_name")
SV_data <- SV_data %>% dplyr::select("Species", "order", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure", "image.files", "description", "Source")  

# add order, family, genus to the data
# SV_data_orders_ncbi <- SV_data %>% mutate(order = ifelse(is.na(order), purrr::map_chr(unique_name, get_order), order))

# fishbase_species <- rfishbase::load_taxa()
# SV_data <- SV_data %>% left_join(fishbase_species, by = "Species")
# SV_data <- SV_data %>% select("Species", "Order", "Family", "Genus", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure", "image.files", "description", "Source")

get_rank <- function(tax_info, rank_name) {
  lineage <- tax_lineage(tax_info)[[1]]
  values <- lineage$name[lineage$rank == rank_name]
  if (length(values) == 0) return(NA_character_) else return(values[1])
}

# fill in taxonomy
df <- SV_data %>%
  rowwise() %>%
  mutate(
    tax_info = list(taxonomy_taxon_info(ott_id, include_lineage = TRUE)),
    order = get_rank(tax_info, "order"),
    family = get_rank(tax_info, "family"),
    genus = get_rank(tax_info, "genus")
  ) %>%
  ungroup() %>%
  dplyr::select(-tax_info)
SV_data <- df %>% rename("Order" = "order") %>% rename("Genus" = "genus") %>% rename("Family" = "family")
SV_data <- SV_data %>% dplyr::select("Species", "Order", "Family", "Genus", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure", "image.files", "description", "Source")

missing_tax <- SV_data %>% filter(is.na(Genus) | is.na(Family) | is.na(Order)) # check where there are problems # 19 
# fill in the gaps with FishBase
fishbase_species <- rfishbase::load_taxa()
missing_tax <- missing_tax %>% left_join(fishbase_species, by = "Species")
missing_tax <- missing_tax %>% rename("Order" = "Order.y") %>% rename("Genus" = "Genus.y") %>% rename("Family" = "Family.y")
missing_tax <- missing_tax %>% dplyr::select("Species", "Order", "Family", "Genus", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure", "image.files", "description", "Source")

SV_data <- SV_data %>% filter(!is.na(Genus)) %>% filter(!is.na(Family)) %>% filter(!is.na(Order)) #filter out rows with missing taxonomy values 
SV_data <- SV_data %>% bind_rows(missing_tax) # add the FishBase correced rows
SV_data <- SV_data %>% filter(!is.na(Genus)) %>% filter(!is.na(Family)) %>% filter(!is.na(Order)) #remove any rows where FishBase couldn't solve the msising value

write.csv(SV_data, "/Users/juliamaja/Desktop/SV/SV_data.csv") 
SV_data <- read.csv("/Users/juliamaja/Desktop/SV/SV_data.csv")

SV_data_avg <- SV_data %>% group_by(Species)
SV_data_avg <- SV_data_avg %>% mutate(Species = as.character(Species), SV = as.numeric(SV))
SV_data_avg <- SV_data_avg %>% mutate(SV = ifelse(SV == 0, NA, SV))


#SV_data_avg <- SV_data_avg %>% mutate(test = !(any(SV > 0) & SV ==0))
#SV_data_avg <- SV_data_avg %>% filter(!(any(SV > 0) & SV == 0)) 

SV_data_avg <- SV_data_avg %>% mutate(SV2 = median(SV, na.rm = TRUE))
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(SV3 = ifelse(all(is.na(SV)) & presence == "absent",0,SV2))
#SV_avg = 
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(SV3 = ifelse( presence == "absent" & any(presence == "present"), NA, SV3)) 
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(presence = ifelse( presence == "absent" & any(presence == "present"), "present", presence ))
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(SV3 = ifelse(is.na(SV3) & any(!is.na(SV3)), first(na.omit(SV3)), SV3))

# SV_data_avg <- SV_data_avg %>% mutate(SV = ifelse(present == "present" & is.na(SV), "present", SV))
#SV_data_avg <- SV_data_avg %>% rename("SV" = "SV3")

SV_data_avg <- SV_data_avg %>% dplyr::select("Species", "Order", "Family", "Genus", "SV3", "presence", "genome.assembly", "ott_id", "flags", "tips", "SV.not.in.figure") 
#dealing with cases where SV is not seen in figure:
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(SV.not.in.figure = if (any(SV.not.in.figure != "y")) ifelse(SV.not.in.figure == "y", "", SV.not.in.figure) else SV.not.in.figure) %>% ungroup()
SV_data_avg <- SV_data_avg %>% distinct() 
SV_data_avg <- SV_data_avg %>% ungroup()
SV_data_avg <- SV_data_avg %>% mutate(tips = str_replace(tips, " ", "_")) # as of april 17: 478 species. as of june 6: 592


# duplicates because of diff orders for same species: 
# undefined_orders <- SV_data_avg %>% group_by(Species) %>% filter(n_distinct(order) > 1)
# ohh this is beacause initially some species have orders in one entry but not in others! perhaps remove all those and redo orders with fishbase
# solved

# check for other duplicates in species names:

#SV_data_avg <- SV_data_avg %>% filter(!duplicated(SV_data_avg$Species))
SV_duplicates <- SV_data_avg[duplicated(SV_data_avg$Species),]
View(SV_duplicates)

# SV_data_avg <- SV_data_avg %>% filter(SV.not.in.figure != "y")

write.csv(SV_data_avg, "/Users/juliamaja/Desktop/SV/SV_data_avg.csv") 

###

# SV_data_avg_1 <- readRDS("/Users/juliamaja/Desktop/SV_data_avg.RDS") 
# SV_data_avg_1 <-SV_data_avg_1 %>% dplyr::select(-climate)
# 
# SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# new_tips <- setdiff(SV_data_avg_1$tips, SV_data_avg$tips)
# new_rows <- SV_data_avg_1 %>% filter(tips %in% new_tips)
# SV_data_avg_full <- bind_rows(SV_data_avg, new_rows)
# write.csv(SV_data_avg_full, "/Users/juliamaja/Desktop/SV/SV_data_avg.csv") 

# extra stats and stuff ---------------------------------------------------


# some stats:
num_w_genome <- SV_data_avg %>% filter(genome.assembly == "y") %>% count() #171 #205 as of july 18 # 171 as of july 23 ? but I didn't change anything
absent_w_genome <- SV_data_avg %>% filter(genome.assembly == "y" & presence == "absent") %>% count() #42 #50
false_absent_w_genome <- SV_data_avg %>% filter(genome.assembly == "y" & presence == "absent" & SV.not.in.figure == "y") %>%  count() #13 #11

library(tidyr)
# Count species per order in your dataset
your_order_counts <- SV_data_avg %>%
  count(Order, name = "n_sampled")

# Count species per order in FishBase
fishbase_order_counts <- fishbase_species %>%
  count(Order, name = "n_total")

# Join the tables and calculate proportion sampled
order_coverage <- left_join(fishbase_order_counts, your_order_counts, by = "Order") %>%
  mutate(
    n_sampled = replace_na(n_sampled, 0),
    prop_sampled = n_sampled / n_total
  )
library(ggplot2)

ggplot(order_coverage, aes(x = n_total, y = n_sampled, label = Order)) +
  geom_point() +
  geom_text(check_overlap = TRUE, hjust = 1.1, vjust = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    x = "Total species in FishBase (per Order)",
    y = "Sampled species in your dataset",
    title = "Order-Level Sampling Coverage"
  ) +
  theme_minimal()

sp_w_gen <- SV_data_avg %>% filter(genome.assembly == "y")
sp_w_gen <- sp_w_gen %>% dplyr::select(Species)
sp_w_gen <- sp_w_gen %>%
  mutate(Species = str_replace(Species, " ", "_"))
write.csv(sp_w_gen, "/Users/juliamaja/Desktop/species_to_keep.csv")

sp_w_genomes_205 <- read.csv("/Users/juliamaja/Desktop/species_to_keep.csv")

#mystery_sp <- anti_join(sp_w_genomes_205, sp_w_gen)

