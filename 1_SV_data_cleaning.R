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
sp_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
trait.data <- sp_data[,c("species", "SV", "order")]
trait.data$tips <- trait.data$species
trait.data$tips <- str_replace(trait.data$species, pattern = " ", replacement = "_")

# resolve spcies names with open tree of life
otol_names <- tnrs_match_names(names = trait.data$species, context_name = "Vertebrates", do_approximate_matching = TRUE) 
resolved_names <- (otol_names[!(is.na(otol_names$unique_name)),] 
                   %>%  rename( "species" = "search_string") 
                   %>% mutate(species = str_replace(species, " ", "_"))
                   %>% select(species, unique_name, "ott_id", "flags")
)
resolved_names <- distinct(resolved_names, ott_id, .keep_all = TRUE)


SV_data <- sp_data %>% mutate(species = str_replace(species, " ", "_"))
SV_data <- SV_data %>% mutate(species = str_replace(species, " ", "_"))
SV_data <- SV_data %>% mutate(species = tolower(species))
SV_data <- SV_data %>% group_by(species) 
SV_data <- SV_data %>% mutate(
  genome.assembly = if_else(
    genome.assembly == '' & any(genome.assembly == 'y'),  
    'y',  
    genome.assembly 
  )
) 
SV_data <- SV_data %>% ungroup()
SV_data <- SV_data %>% left_join(resolved_names, by = "species")  
SV_data <- SV_data %>% mutate(tips = unique_name)
SV_data <- SV_data %>% rename("Source" = "source")
SV_data <- SV_data %>% filter(Source != "")
SV_data <- SV_data %>% mutate(SV = str_replace(SV, "2", "1"), str_replace(SV, "0.5", "1"), str_replace(SV, "2.5", "2"), str_replace(SV, "3.5", "4"))
SV_data <- SV_data %>% filter(unique_name !="")
SV_data <- SV_data %>% select("species", "unique_name", "order", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips")
SV_data <- SV_data %>% filter(SV != "")
SV_data <- SV_data %>% mutate(presence = ifelse(SV == "0", "absent", "present"))
SV_data <- SV_data %>% mutate(SV = ifelse(SV == "present", NA, SV))
SV_data <- SV_data %>% rename("Species" = "unique_name")
SV_data <- SV_data %>% select("Species", "order", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips")

# add order, family, genus to the data
# SV_data_orders_ncbi <- SV_data %>% mutate(order = ifelse(is.na(order), purrr::map_chr(unique_name, get_order), order))

fishbase_species <- rfishbase::load_taxa()
SV_data <- SV_data %>% left_join(fishbase_species, by = "Species")
SV_data <- SV_data_orders %>% select("Species", "Order", "Family", "Genus", "SV", "presence", "genome.assembly", "ott_id", "flags", "tips")


SV_data_avg <- SV_data %>% group_by(Species) 
SV_data_avg <- SV_data_avg %>% mutate(Species = as.character(Species), SV = as.numeric(SV))
SV_data_avg <- SV_data_avg %>% mutate(SV = ifelse(SV == 0, NA, SV))
#SV_data_avg <- SV_data_avg %>% mutate(test = !(any(SV > 0) & SV ==0))
#SV_data_avg <- SV_data_avg %>% filter(!(any(SV > 0) & SV == 0)) 

SV_data_avg <- SV_data_avg %>% mutate(SV2 = median(SV, na.rm = TRUE))
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(SV3 = ifelse(all(is.na(SV)) & presence == "absent",0,SV2))
SV_data_avg <- SV_data_avg %>% group_by(Species) %>% mutate(presence = ifelse( presence == "absent" & any(presence == "present"), "present", presence ))

SV_data_avg <- SV_data_avg %>% select("Species", "Order", "Family", "Genus", "SV3", "presence", "genome.assembly", "ott_id", "flags", "tips") 
SV_data_avg <- SV_data_avg %>% distinct() 
SV_data_avg <- SV_data_avg %>% ungroup()
SV_data_avg <- SV_data_avg %>% mutate(tips = str_replace(tips, " ", "_"))

# duplicates because of diff orders for same species: 
# undefined_orders <- SV_data_avg %>% group_by(Species) %>% filter(n_distinct(order) > 1)
# ohh this is beacause initially some species have orders in one entry but not in others! perhaps remove all those and redo orders with fishbase
# solved

# check for other duplicates in species names:
SV_duplicates <- SV_data_avg[duplicated(SV_data_avg$Species),]
View(SV_duplicates)

write.csv(SV_data_avg, "/Users/user2/Desktop/SV_project/SV_data_avg.csv") 




