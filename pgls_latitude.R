library(tidyverse)
library(phytools)
library(ggplot2)



SV_data_climate <- readRDS("/Users/juliamaja/Desktop/SV/SV_data_climate_progress.RDS")


tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") 

SV_data_climate <- SV_data_climate %>% filter(climate != "")

# check for NA values
any(is.na(SV_data_climate$climate))

SV_data_climate <- SV_data_climate %>% select(-"X.1", -"X")

# distribution of climates
climate_counts <- SV_data_climate %>%
  group_by(climate) %>%
  summarise(count = n(), .groups = "drop")
pie(climate_counts$count, climate_counts$climate, main = "distribution of climate types in my data")

# presence/ absence by climate
climate_counts <- SV_data_climate %>%
  group_by(presence, climate) %>%
  summarise(count = n(), .groups = "drop")

climate_counts_present <- climate_counts %>% filter(presence == "present")
plot <- pie(climate_counts_present$count, climate_counts_present$climate, main = "distribution of climate types for SV present")

climate_counts_absent <- climate_counts %>% filter(presence == "absent")
plot <- pie(climate_counts_absent$count, climate_counts_absent$climate, main = "distribution of climate types for SV absent")

# SV complexity by climate

SV_data_climate <- SV_data_climate %>%
  mutate(SV3_group = case_when(
    SV3 == 0 ~ "0",
    SV3 > 0 & SV3 < 3 ~ "0.5–2.99",
    SV3 >= 3 & SV3 <= 4 ~ "3–4"
  ))

climate_counts <- SV_data_climate %>%
  group_by(SV3_group, climate) %>%
  summarise(count = n(), .groups = "drop")

as.character(climate_counts$SV3_group)

climate_counts_0 <- climate_counts %>% filter(SV3_group == "0")
plot <- pie(climate_counts_0$count, climate_counts_0$climate, main = "distribution of climate types for SV = 0")

climate_counts_1_3 <- climate_counts %>% filter(SV3_group == "0.5–2.99")
plot <- pie(climate_counts_1_3$count, climate_counts_1_3$climate, main = "distribution of climate types for SV between 1 and 3")

climate_counts_3_4 <- climate_counts %>% filter(SV3_group == "3–4")
plot <- pie(climate_counts_3_4$count, climate_counts_3_4$climate, main = "distribution of climate types for SV between 3 and 4")

####
# pgls for trait category vs climate

library(caper)

SV_data_latitude <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_latitude.csv")
SV_data_latitude$latitude <- as.factor(SV_data_latitude$latitude)
SV_data_latitude <- SV_data_latitude[!is.na(SV_data_latitude$SV3), ]

tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") # my tree

SV_data_climate <- SV_data_latitude
# Set species names as row names
rownames(SV_data_climate) <- SV_data_climate$tips
SV_data_climate <- SV_data_climate %>% filter(!is.na(latitude))
common_species <- intersect(tr$tip.label, rownames(SV_data_climate))
tr <- drop.tip(tr, setdiff(tr$tip.label, common_species))
SV_data_climate <- SV_data_climate[common_species, ]

comp_data <- comparative.data(phy = tr, data = SV_data_climate,
                              names.col = "tips", vcv = TRUE, warn.dropped = TRUE)

#model <- pgls(SV3 ~ latitude, data = comp_data, lambda = 0.0997) # fixed lambda
model <- pgls(SV3 ~ latitude, data = comp_data, lambda = "ML") # letting the model estimate lambda
summary(model)

phylosig(tr, as.vector(SV_data_climate$SV3), method="K", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
# lambda = 0.0997 but I think this is incorrect - i think that because the order of  
# trait values in the vector don't match the order of the tip labels in the tree
# the phylogenetic signal seems close to zero 

SV_data_climate <- SV_data_climate[match(tr$tip.label, SV_data_climate$tips), ]


###

# phylogenetic logistic regresseion for presence/ absence
# logistic regression is used for binary variables (linear regression is for continuous)


library(phylolm)

# model function requires binary data (0/1)
SV_data_climate$presence_num <- ifelse(SV_data_climate$presence == "present", 1, 0)
SV_data_climate <- SV_data_climate %>% filter(!is.na(latitude))

model <- phyloglm(presence_num ~ latitude, phy = tr, data = SV_data_climate,
                  method = "logistic_MPLE")

# alpha: 0.001463066 i think alpha is phylogenetic signal and 
# this is the same problem as before where the order of 
# trait values don't match the order of the tip labels in the tree
# so the model estimates low phylogenetic signal
#
# but when I do SV_data_climate$tips and tr$tip.label the species are in the same order

summary(model)

# output:
# Coefficients:
#   Estimate     StdErr z.value  p.value   
# (Intercept) -0.1880351  0.7586217 -0.2479 0.804240   
# latitude     0.0153949  0.0055526  2.7725 0.005562 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#  latitude 0.0154 means that for every 1° increase in latitude
#  the log-odds of SV presence increases by 0.0154
# 
#  p-value = 0.0056 = statistically significant (**).
# 
#  strong evidence of a positive relationship between latitude and SV presence

# plotting

library(ggplot2)

coef_vals <- coef(model)
intercept <- coef_vals[1]
slope <- coef_vals[2]

SV_data_avg$predicted_logit <- intercept + slope * SV_data_avg$latitude

SV_data_avg$predicted_prob <- 1 / (1 + exp(-SV_data_avg$predicted_logit))

ggplot(SV_data_avg, aes(x = latitude, y = presence_num)) +
  #geom_jitter(height = 0.05, alpha = 0.3) +
  geom_line(aes(y = predicted_prob), color = "blue", size = 1) +
  labs(y = "Predicted Probability of Presence", x = "Latitude") +
  theme_minimal()
plot(presence_num ~ latitude, data = SV_data_avg,
     pch = 16, col = rgb(0, 0, 0, 0.3),
     ylab = "Trait Presence (0/1)", xlab = "Latitude")


ggplot(SV_data_avg, aes(x = presence, y = latitude)) + geom_boxplot()

### 

# phylogenetic ANOVA

SV_data_climate <- SV_data_climate[match(tr$tip.label, SV_data_climate$tips), ]
tr <- keep.tip(tr, tip = SV_data_climate$tips)
x = SV_data_climate$presence
y = SV_data_climate$latitude

phylANOVA(tr, x, y, nsim=1000, posthoc=TRUE, p.adj="holm")

boxplot(y ~ x, 
        ylab = "Latitude", 
        xlab = "Presence", 
        main = "Trait differences between groups",
        )
pval <- 0.016  
text(x = 1, y = max(y), labels = paste("p =", pval))





