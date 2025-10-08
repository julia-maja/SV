library(phylolm)


# SV presence/ absence vs latitude ----------------------------------------
# load in dataset and tree
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") 
# note: tr_tree_calibrated_9_8 is built without 'false absent' species but SV_data_avg has them
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
# tr_tree_calibrated_9_9.rds has all species.

# model function requires binary data (0/1)
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
# it also requires the row names to be species names so it can match the data to the tree
rownames(SV_data_avg) <- SV_data_avg$tips

#### test SV presence/ absence vs latitude ####

# filter species with possible false absents
SV_data_avg <- SV_data_avg %>% filter(is.na(SV.not.in.figure))
tr <- drop.tip(tr, setdiff(tr$tip.label, SV_data_avg$tips))

# filter rows where latitude is NA
SV_data_avg <- SV_data_avg %>% filter(!is.na(latitude))

# apply the logistic regression
model <- phyloglm(presence_num ~ latitude, phy = tr, data = SV_data_avg,
                  method = "logistic_MPLE", boot = 1000, save = TRUE)

summary(model)

hist(model[["bootstrap"]][, "latitude"], breaks = 100)

# output:
# Parameter estimate(s):
#   alpha: 0.001363882 
# 
# Coefficients:
#   Estimate    StdErr z.value  p.value   
# (Intercept) 0.4332899 0.7353995  0.5892 0.555734   
# latitude    0.0093320 0.0035148  2.6551 0.007928 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# still significant even with false absents removed

# Coefficients:
#   Estimate    StdErr   z.value lowerbootCI upperbootCI p.value  
# (Intercept) 0.4432736 0.7827276 0.5663192   0.4434409      0.4447 0.57118  
# latitude    0.0078289 0.0033194 2.3585071   0.0035235      0.0193 0.01835 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### plotting



# SV presence/ absence vs salinity ----------------------------------------


SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
rownames(SV_data_avg) <- SV_data_avg$tips
# filter species with possible false absents
SV_data_avg <- SV_data_avg %>% filter(is.na(SV.not.in.figure))
tr <- drop.tip(tr, setdiff(tr$tip.label, SV_data_avg$tips))

# filter rows where salinity is NA
SV_data_avg <- SV_data_avg %>% filter(!is.na(salinity))
# filter the one species that is "brackish"
SV_data_avg <- SV_data_avg %>% filter(salinity != "brackish")

# apply the logistic regression
model <- phyloglm(presence_num ~ salinity, phy = tr, data = SV_data_avg,
                  method = "logistic_MPLE", boot = 1000, full.matrix = TRUE, save = TRUE)

summary(model)

hist(model[["bootstrap"]][, "salinitysaltwater"], breaks = 100)
# output:
# alpha: 0.004655865 
# 
# Coefficients:
#   Estimate  StdErr z.value  p.value   
# (Intercept)        0.39265 0.49157  0.7988 0.424418   
# salinitysaltwater  0.84130 0.28085  2.9956 0.002739 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# salinity is a significant predictor of SV presence - high salinity = high likelihood of SV presence 
# this agrees with the 1992 paper which fond that most freshwater fish have no or poorly developed SV


### plotting


# SV occurrence vs diel pattern -------------------------------------------
# load in dataset and tree
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") 
# note: tr_tree_calibrated_9_8 is built without 'false absent' species but SV_data_avg has them
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
# tr_tree_calibrated_9_9.rds has all species.

# model function requires binary data (0/1)
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
# it also requires the row names to be species names so it can match the data to the tree
rownames(SV_data_avg) <- SV_data_avg$tips

diel_data <- readRDS("/Users/juliamaja/Desktop/SV/trait_data_fish.rds")

noct_diurn <- diel_data %>% filter(diel == "diurnal" | diel == "nocturnal")

SV_data_noct_diurn <- left_join(SV_data_avg, noct_diurn, by = "tips")
SV_data_noct_diurn <- SV_data_noct_diurn %>% filter(!is.na(diel))
rownames(SV_data_noct_diurn) <- SV_data_noct_diurn$tips

model <- phyloglm(presence_num ~ diel, phy = tr, data = SV_data_noct_diurn,
                  method = "logistic_MPLE", boot = 1000, save = TRUE)

summary(model)
hist(model[["bootstrap"]][, "dielnocturnal"], breaks = 100)

# output:
# Parameter estimate(s):
#   alpha: 0.00637774 
# 
# Coefficients:
#   Estimate  StdErr z.value  p.value   
# (Intercept)    0.23401 0.44538  0.5254 0.599299   
# dielnocturnal  0.67187 0.25863  2.5978 0.009382 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### plotting

#################################################
#################################################
############## crepuscular vs not ###############
#################################################
#################################################

# load in dataset and tree
SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
# tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_8.rds") 
# note: tr_tree_calibrated_9_8 is built without 'false absent' species but SV_data_avg has them
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
# tr_tree_calibrated_9_9.rds has all species.

# model function requires binary data (0/1)
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
# it also requires the row names to be species names so it can match the data to the tree
rownames(SV_data_avg) <- SV_data_avg$tips

diel_data <- readRDS("/Users/juliamaja/Desktop/SV/trait_data_fish.rds")

crep <- diel_data %>% filter(!is.na(crepuscular))

SV_data_crep <- left_join(SV_data_avg, crep, by = "tips")
rownames(SV_data_crep) <- SV_data_crep$tips
SV_data_crep <- SV_data_crep %>% filter(!is.na(crepuscular))

model <- phyloglm(presence_num ~ crepuscular, phy = tr, data = SV_data_crep,
                  method = "logistic_MPLE", boot = 1000, full.matrix = TRUE, save = TRUE)

summary(model)

hist(model[["bootstrap"]][, "crepuscularnon_crepuscular"], breaks = 100)

# output:
# Parameter estimate(s):
#   alpha: 0.004655612 
# 
# Coefficients:
#   Estimate  StdErr z.value p.value  
# (Intercept)                 0.42898 0.67977  0.6311 0.52800  
# crepuscularnon_crepuscular  0.60286 0.35851  1.6816 0.09265 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### plotting



# bootstrapping -----------------------------------------------------------
### control ###

SV_data_avg <- read.csv("/Users/juliamaja/Desktop/SV/SV_data_avg.csv")
tr <- readRDS("/Users/juliamaja/Desktop/SV/tr_tree_calibrated_9_9.rds") 
SV_data_avg$presence_num <- ifelse(SV_data_avg$presence == "present", 1, 0)
rownames(SV_data_avg) <- SV_data_avg$tips
# filter species with possible false absents
SV_data_avg <- SV_data_avg %>% filter(is.na(SV.not.in.figure))
tr <- drop.tip(tr, setdiff(tr$tip.label, SV_data_avg$tips))

set.seed(12345)
SV_data_avg$variable <- sample(c(0, 1), size = nrow(SV_data_avg), replace = TRUE)


# apply the logistic regression
model <- phyloglm(presence_num ~ variable, phy = tr, data = SV_data_avg,
                  method = "logistic_MPLE")

summary(model)

# seed: 123
# output: this random variable is significantly correlated :(
# AIC     logLik Pen.logLik 
# 420.6     -207.3     -204.5 
# 
# Method: logistic_MPLE
# Mean tip height: 688.9458
# Parameter estimate(s):
#   alpha: 0.004907834 
# 
# Coefficients:
#   Estimate   StdErr z.value  p.value   
# (Intercept)  0.83662  0.49737  1.6821 0.092550 . 
# variable    -0.41067  0.14037 -2.9255 0.003439 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Note: Wald-type p-values for coefficients, conditional on alpha=0.004907834

# seed: 1678: output:
# Coefficients:
#   Estimate  StdErr z.value   p.value    
# (Intercept)  1.24653 0.33511  3.7198 0.0001994 ***     # no idea what it means for this to be significant
#   variable     0.22285 0.16478  1.3524 0.1762509    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# random variable simulation test
set.seed(123)
n_sim <- 1000
p_vals <- numeric(n_sim)

for (i in 1:n_sim) {
  SV_data_avg$random_var <- rbinom(nrow(SV_data_avg), 1, 0.5)
  model_rand <- phyloglm(presence_num ~ random_var, phy = tr, data = SV_data_avg,
                         method = "logistic_MPLE")
  p_vals[i] <- summary(model_rand)$coefficients[2, 4]  # get p-value for random_var
}
mean(p_vals < 0.05) # this gives the false positive rate i.e. the proportion of p's less than 0.05 out of 100 permutations

# if p > 5%, then the model is overfitting or assumptions are violated
# p = 0.35 so 35 out of 100 random simulations are significant. yikes!

# random variable simulation test
set.seed(123)
n_sim <- 100
p_vals <- numeric(n_sim)

for (i in 1:n_sim) {
  SV_data_avg$random_var <- sample(-100:100, size = nrow(SV_data_avg), replace = TRUE)
  model_rand <- phyloglm(presence_num ~ random_var, phy = tr, data = SV_data_avg,
                         method = "logistic_MPLE")
  p_vals[i] <- summary(model_rand)$coefficients[2, 4]  # get p-value for random_var
}
mean(p_vals < 0.05)
# p = 0.94 

# permutation test for salinity
set.seed(456)
n_perm <- 1000
p_vals_perm <- numeric(n_perm)

for (i in 1:n_perm) {
  SV_data_avg$sal_perm <- sample(SV_data_avg$salinity)
  model_perm <- phyloglm(presence_num ~ sal_perm, phy = tr, data = SV_data_avg,
                         method = "logistic_MPLE")
  p_vals_perm[i] <- summary(model_perm)$coefficients[2, 4]
}
mean(p_vals_perm < 0.05) 
hist(p_vals_perm)
# p = 0.03

# permutation test for latitude
set.seed(456)
n_perm <- 100
p_vals_perm <- numeric(n_perm)

for (i in 1:n_perm) {
  SV_data_avg$lat_perm <- sample(SV_data_avg$latitude)
  model_perm <- phyloglm(presence_num ~ lat_perm, phy = tr, data = SV_data_avg,
                         method = "logistic_MPLE")
  p_vals_perm[i] <- summary(model_perm)$coefficients[2, 4]
}
mean(p_vals_perm < 0.05)
# p = 0.76


# try the bootstrapping method that is implemented in the package.



