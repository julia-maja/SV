library(tidyverse)
library(ggplot2)
library(virdis)
library(ggsankey)

SV_data <- read.csv("/Users/juliamaja/Desktop/SV/SV_data.csv")
SV_data <- SV_data %>%
  mutate(data_type = case_when(
    image.files != "" & description != "" ~ "photo and description",
    image.files != "" ~ "photo",
    description != "" ~ "description",
    TRUE ~ NA_character_ 
  ))

SV_data <- SV_data %>%
  mutate(
    photo_type = paste(
      ifelse(str_detect(image.files, regex("ventral", ignore_case = TRUE)), "ventral", ""),
      ifelse(str_detect(image.files, regex("lateral", ignore_case = TRUE)), "lateral", ""),
      ifelse(str_detect(image.files, regex("sagittal", ignore_case = TRUE)), "sagittal section", ""),
      ifelse(str_detect(image.files, regex("cross", ignore_case = TRUE)), "cross section", ""),
      sep = ", "
    ),
   photo_type = str_replace_all(photo_type, "^,\\s*|,\\s*$", ""), # remove leading/trailing commas
   photo_type = str_replace_all(photo_type, "^,\\s*|", ""), 
   photo_type = str_replace_all(photo_type, "^,\\s*|", ""),
   photo_type = str_replace_all(photo_type, ",\\s*$", ""),
   photo_type = str_replace_all(photo_type, ",\\s*$", ""),
   photo_type = str_replace_all(photo_type, ",\\s*,", ","),
   photo_type = str_replace_all(photo_type, ",\\s*,", ",")  # clean up double commas
  ) %>%
  mutate(photo_type = ifelse(photo_type == "", NA, photo_type))  # set blank results to NA

#creating a confusion plot
SV_sp_concordance <- SV_data %>% select(.,  c("Species", "presence", "data_type", "photo_type"))

SV_sp_concordance1 <- SV_sp_concordance %>% filter(data_type == "photo and description") %>% separate(., data_type, into = c("data_type1", "data_type2"), sep = " and ")
SV_sp_concordance1 <- pivot_longer(data = SV_sp_concordance1, cols = c(data_type1, data_type2), names_to = "column", values_to = "data_type")
SV_sp_concordance1 <-SV_sp_concordance1 %>% select("Species", "data_type", "photo_type", "presence")
SV_sp_concordance <- bind_rows(SV_sp_concordance1, SV_sp_concordance)
  
SV_sp_concordance2 <- SV_sp_concordance[!is.na(SV_sp_concordance$photo_type), ] %>% separate(., photo_type, into = c("photo_type1", "photo_type2", "photo_type3", "photo_type4"), sep = ", ")
SV_sp_concordance2[SV_sp_concordance2 == ""] <- NA
SV_sp_concordance2 <- pivot_longer(data = SV_sp_concordance2, cols = c(photo_type1, photo_type2, photo_type3, photo_type4), names_to = "column", values_to = "photo_type")
SV_sp_concordance2 <- SV_sp_concordance2[!is.na(SV_sp_concordance2$photo_type),]
SV_sp_concordance2 <-SV_sp_concordance2 %>% select("Species", "photo_type", "presence")
colnames(SV_sp_concordance2) <- c("Species", "data_type", "presence")

SV_sp_concordance <-SV_sp_concordance %>% select("Species", "data_type", "presence")
SV_sp_concordance <- SV_sp_concordance %>% filter(data_type == "description")

SV_sp_concordance <- bind_rows(SV_sp_concordance2, SV_sp_concordance)

multiple_sources <- SV_sp_concordance %>% count(Species) %>% filter(n>1)
SV_sp_concordance <- SV_sp_concordance[SV_sp_concordance$Species %in% multiple_sources$Species, ]

write.csv(SV_sp_concordance, "/Users/juliamaja/Desktop/SV/SV_sp_concordance.csv", row.names = FALSE)

species_list <- table(SV_sp_concordance$Species)
species_list <- names(species_list[species_list >1])

#function Max wrote for comparing entries
compTwo <- function(comp1 = "comp1", comp2 = "comp2") {
  
  if(any(is.na(c(comp1, comp2)))) {
    return(NA)
  } else {
    #then compares if any of the components match
    if(any(comp1 %in% comp2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
}

#species <- "Latimeria chalumnae"

#apply this function across all species with multiple entries
output <- lapply(species_list, function(species) {
  
  #filter for one species at a time
  df <- SV_sp_concordance[SV_sp_concordance$Species == species,]
  #rename the column names to be unique for every entry for this species (ie for multiple column 2 entries column 2.1, 2.2 etc)
  df$data_type <- make.unique(df$data_type)
  
  #converts the dataframe so it compares every entry with each other (ie for A,B,C A-A, A-B, A-C, B-A, B-B, B-C, etc)
  df_lists_comb <- expand(df, nesting(var = data_type, vector = presence), nesting(var2 = data_type, vector2 = presence), .name_repair = "universal")
  
  #??? idk
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  #evaluates the activity patterns for each of these sources and returns if they agree or not (TRUE or FALSE)
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  #manipulate the strings for both variable names to revert them back to the original name (back to column 2 from col 2.1)
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 5)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 5)
  
  #create a column returning the comparison being made (ie col2-col2, col1-col2, etc)
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "_")
  
  #return just the comparison result column (TRUE or FALSE match) and the comparison being made (ie col1 vs col1)
  return(comparisons[,c("comp","var_final")])
})

#combine this list of results 
output <- Reduce(rbind, output)

table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "_"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "_"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table
plot_freq <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + 
  geom_tile() + geom_text() #+ scale_fill_viridis(limits = c(0,1))
plot_freq
plot_count <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) +
  geom_tile() + geom_text() #+ scale_fill_viridis(limits = c(0,1))
plot_count


# sankey ------------------------------------------------------------------
SV_sp_concordance <- read.csv("/Users/juliamaja/Desktop/SV/SV_sp_concordance.csv")

species_list <- table(SV_sp_concordance$Species)
species_list <- names(species_list[species_list >1])

output <- lapply(species_list, function(species){
  df <- SV_sp_concordance[SV_sp_concordance$Species == species, ]
  df$data_type <- make.unique(df$data_type)
  return(df)
})

SV_sp_concordance <- Reduce(rbind, output)

SV_wide <- pivot_wider(SV_sp_concordance, names_from = "data_type", values_from = "presence")

#all entires comparison
# SV_wide <- SV_wide %>% make_long(2:19)
#first instance of ventral, lateral, sagital or description source data
# SV_wide <- SV_wide %>% make_long(c(2,3,4,7))
#how well do all ventral agree
# SV_wide <- SV_wide %>% make_long(c("ventral", "ventral.1", "ventral.2", "ventral.3"))
SV_wide <- SV_wide %>% make_long(c("cross section", "description"))

SV_wide <- SV_wide %>% filter_at(vars(node, next_node), any_vars(!is.na(.)))

ggplot(SV_wide, aes(x=x, next_x = next_x, node = node,  next_node = next_node, fill = factor(node))) + geom_sankey() + theme_sankey(base_size = 16)
