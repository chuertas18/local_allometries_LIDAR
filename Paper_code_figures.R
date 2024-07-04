# Clear the workspace
rm(list = ls())

# Load required libraries
library(data.table)   # For fast data manipulation
library(ggplot2)      # For data visualization
library(brms)         # For Bayesian regression modeling
library(ggpubr)       # For publication ready plots
library(ggExtra)      # For adding marginal histograms to ggplot2, and more
library(tidyverse)    # Includes ggplot2, dplyr, tidyr, readr, etc.

# Working directory
setwd("./Data/")

########################################################################################
# Figure 5 - Boxplot of individual tree height estimation error in meters per plot using 
# a pantropical equation (Chave et al. 2014), the reference model 4 - MM and the canopy 
# constructor predictions (mean of 5 runs). The dotted line represents the 0 error. 
#######################################################################################
###################################################################################################
# 1/ Databases
###################################################################################################

# Read allometric database and filter it
data_allo <- fread("Dataset1") 
data_allo <- data_allo[outliers == 0 & duplicates_als == 0, ]
data_allo$src <- as.factor(data_allo$src) 
data_allo <- data_allo[src == "ALS", ]
summary(data_allo)

# Subset the data
newdata <- data_allo[operator == "CIRAD", .(idtree, sel_allom, plot, dbh_allo, dbh_year, H_year, trait, species_allo, HC, h95)]

# Function to calculate expected height from pantropical model
pantropical_Chave14 <- function(DBH) {
  E <- -0.1074988
  H <- exp(0.893 - E + 0.760 * log(DBH) - 0.0340 * (log(DBH))^2) * exp(0.243^2 / 2)
  return(H)
}

# Load height prediction data from different versions
v6 <- fread("P1to16_0006_trees_final.txt")
v7 <- fread("P1to16_0007_trees_final.txt")
v8 <- fread("P1to16_0008_trees_final.txt")
v9 <- fread("P1to16_0009_trees_final.txt")
v10 <- fread("P1to16_0010_trees_final.txt")

# Merge height prediction data
cc <- v6[, .(id, dbh, height)]
setnames(cc, "height", "H6")
cc <- merge(cc, v7[, .(id, height)], by = "id")
setnames(cc, "height", "H7")
cc <- merge(cc, v8[, .(id, height)], by = "id")
setnames(cc, "height", "H8")
cc <- merge(cc, v9[, .(id, height)], by = "id")
setnames(cc, "height", "H9")
cc <- merge(cc, v10[, .(id, height)], by = "id")
setnames(cc, "height", "H10")

# Calculate correlation and predicted height
cor(cc[, 3:7])
cc[, predCC := rowMeans(.SD), .SDcols = 3:7]

# Merge new data with height prediction data
newdata <- merge(newdata, cc[, .(id, dbh, predCC)], by.x = "idtree", by.y = "id")

# Load CIRAD data and species list
# Data from CIRAD plots (Dataset2), available for consultation with permission:
#   Derroire, G., Hérault, B., Rossi, V., Blanc, L., Gourlet- Fleury, S., & Schmitt, L. (2023). Paracou permanent plots. CIRAD Dataverse. https://doi.org/10.18167/DVN1/8G8AHY
data_CIRAD <- fread("Dataset2.csv")
not_sing <- names(table(data_allo$species_allo))[which(table(data_allo$species_allo) > 1)]

# Update species names to match unified nomenclature
data_CIRAD$species_allo <- paste(data_CIRAD$genus, data_CIRAD$species, sep = ' ')
data_CIRAD$species_allo[data_CIRAD$species_allo == "Agonandra silvatica"] <- "Agonandra sylvatica"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Alexa wachenheimii"] <- "Alexa wachenheimi"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Eschweilera grandiflora_form2"] <- "Eschweilera grandiflora"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Inga capitata_form2"] <- "Inga capitata"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Lecythis persistens subsp. aurantiaca"] <- "Lecythis persistens"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Licaria cannella"] <- "Licaria canella"
data_CIRAD$species_allo[data_CIRAD$species_allo == "Parinari rodolphii"] <- "Parinari rodolphi"

# Mark singletons
data_CIRAD$singleton <- !data_CIRAD$species_allo %in% not_sing
data_CIRAD$species_noSing <- data_CIRAD$species_allo
data_CIRAD$species_noSing[data_CIRAD$singleton] <- "bin"

# Select and rename columns for model prediction
sel_data_CIRAD <- data_CIRAD[data_CIRAD$idtree %in% newdata$idtree, .(idtree, dbh, species_noSing, CHM2015_25ppm_fo_comp)]
setnames(sel_data_CIRAD, c("idtree", "dbh_allo", "species_noSing", "HC"))

# Load models and predict heights
model4_MM_nsing <- readRDS("model4_MM_nsing_field.rds")
sel_data_CIRAD$predMM4 <- fitted(model4_MM_nsing, newdata = sel_data_CIRAD)[, 1]

model1_MM <- readRDS("model1_MM_field.rds")
sel_data_CIRAD$predMM1 <- fitted(model1_MM, newdata = sel_data_CIRAD)[, 1]

# Merge predicted heights with new data
newdata <- merge(newdata, sel_data_CIRAD, by.x = "idtree", by.y = "id")
newdata <- newdata[, .(idtree, plot, dbh_allo.x, trait, species_noSing, h95, predCC, predMM4, predMM1)]
setnames(newdata, c("idtree", "plot", "dbh", "trait", "species_noSing", "h_OBS", "h_CC", "h_MM4", "h_MM1"))

# Calculate pantropical predicted height
newdata$h_PAN <- pantropical_Chave14(newdata$dbh)

# Create scatter plot matrix
pairs(newdata[, .(h_OBS, h_CC, h_MM4, h_MM1, h_PAN)], pch = 19)

# Create boxplot data frames
df1 <- data.frame(plot = newdata$plot, model = "CC", error = newdata$h_CC - newdata$h_OBS)
df2 <- data.frame(plot = newdata$plot, model = "Pan", error = newdata$h_PAN - newdata$h_OBS)
df3 <- data.frame(plot = newdata$plot, model = "MM4", error = newdata$h_MM4 - newdata$h_OBS)
df4 <- data.frame(plot = newdata$plot, model = "MM1", error = newdata$h_MM1 - newdata$h_OBS)

df <- rbind(df1, df2, df3, df4)

# Create boxplot
height_boxplot <- ggplot(df, aes(x = as.factor(plot), y = error, fill = model)) +
  geom_boxplot(outlier.colour = NULL, outlier.shape = 16, outlier.size = 1, notch = TRUE) +
  labs(x = "Plot number", y = "Predicted height - Observed height (m)") +
  geom_hline(aes(yintercept = 0), color = "black", linetype = "twodash") +
  scale_fill_discrete(labels = c("Canopy Constructor", "Bayesian model", "Pantropical model")) +
  theme_light() +
  theme(text = element_text(size = 15), legend.position = "bottom")

print(height_boxplot)

# Calculate mean and standard deviation of errors
mean(df1$error)
mean(df2$error)
mean(df3$error)
mean(df4$error)
sd(df1$error)
sd(df2$error)
sd(df3$error)
sd(df4$error)

# Create scatter plots
ggplot(newdata, aes(x = h_CC, y = h_OBS)) + 
  geom_point() +
  geom_smooth(method = lm) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1)

ggplot(newdata, aes(x = h_MM4, y = h_OBS)) + 
  geom_point() +
  geom_smooth(method = lm) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1)

plot(h_OBS ~ h_CC, data = newdata)

# Create density scatter plots for high density data
smoothScatter(newdata$h_OBS ~ newdata$dbh)
x <- newdata[, .(h_OBS, dbh)]
plot(x, col = densCols(x), pch = 20)

x <- cc[, .(dbh, predCC)]
plot(x, col = densCols(x), pch = 20)
smoothScatter(x$predCC ~ x$dbh)

######################################################################################
## Figure 6 - Per species H-DBH relationships of trees with DBH >20 cm belonging to 
# the 6 most abundant species in dataset 1 (Hobsv, green line). The relationship is also 
# plotted after replacing observed heights with CC-derived heights (H ̂CC, extracted from 
#  dataset 2, purple lines). Lines are based on a Generalized Additive Model 
# (geom_smooth method=loess in ggplot2)
######################################################################################
# Loading allometric data and performing initial filtering
data_allo <- fread("Dataset1")
data_allo <- data_allo[outliers == 0 & duplicates_als == 0, ]
data_allo$src <- as.factor(data_allo$src) # Converting 'src' to a factor
data_allo <- data_allo[src == "ALS", ]

# Loading and combining inventory data from multiple sources
inventory_CC_files <- list(
  "P1to16_0006_trees_final.txt",
  "P1to16_0007_trees_final.txt",
  "P1to16_0008_trees_final.txt",
  "P1to16_0009_trees_final.txt",
  "P1to16_0010_trees_final.txt"
)
inventory_CC_all_data <- rbindlist(lapply(inventory_CC_files, fread))

# Aggregating data by 'id' and calculating mean for specified columns
inventory_CC <- inventory_CC_all_data[, .(mean_dbh = mean(dbh, na.rm = TRUE), 
                                          mean_height = mean(height, na.rm = TRUE)), by = id]

# Merging and further processing of allometric and inventory data
data_allo_cc <- merge(data_allo[, .(idtree, dbh_allo, dbh, H, trait, habitat, species_allo, parcela)], inventory_CC, by = "idtree")
data_allo_cc <- data_allo_cc[dbh >= 20, ]

# Generating summary table of species data
species_summary <- data_allo_cc %>%
  group_by(species_allo) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

print(species_summary) # Print the summary table

# Plotting
fig6 <- ggplot(data_allo_cc, aes(x = dbh, y = H)) +
  geom_point(aes(color = species_allo)) +
  geom_smooth(method = "loess", se = FALSE, aes(color = species_allo)) +
  labs(x = "DBH (cm)", y = "Height (m)", color = "Species") +
  theme_minimal()
print(fig6)

######################################################################################
## Figure 7 - Box plot of height residuals across habitat units (upper row) and 
## forest types (lower row) for different models (1-MM to 4-MM) applied to CC predictions.
######################################################################################

# Load BRM Models
load("model1_MM_dataCCmean.Rdata")
model1_MM_nsing_CC$formula

load("model2_MM_dataCCmean.Rdata")
model2_MM_nsing_CC$formula

load("model3_MM_dataCCmean.Rdata")
model3_MM_nsing_CC$formula

load("model4_MM_dataCCmean.Rdata")
fit2_MM4_5s_2023$formula

# Prediction and Evaluation Function
evaluate_model <- function(model, model_name, data_path) {
  preds <- predict(model)
  lm_model <- lm(model$data$height ~ preds[, 1])
  rmse <- round(sqrt(mean(lm_model$residuals^2)), 2)
  write.csv(preds, data_path)
  return(list(preds = preds, lm_model = lm_model, rmse = rmse))
}

# Evaluate all models
eval_MM1 <- evaluate_model(model1_MM_nsing_CC, "Model 1", "PredMM1.csv")
eval_MM2 <- evaluate_model(model2_MM_nsing_CC, "Model 2", "PredMM2.csv")
eval_MM3 <- evaluate_model(model3_MM_nsing_CC, "Model 3", "PredMM3.csv")
eval_MM4 <- evaluate_model(fit2_MM4_5s_2023, "Model 4", "PredMM4.csv")

# Load training data and add residuals
# Data from CIRAD plots (Dataset2), available for consultation with permission:
#   Derroire, G., Hérault, B., Rossi, V., Blanc, L., Gourlet- Fleury, S., & Schmitt, L. (2023). Paracou permanent plots. CIRAD Dataverse. https://doi.org/10.18167/DVN1/8G8AHY
train <- fread("Dataset2_train29944.csv")
train[, `:=`(
  res1 = eval_MM1$lm_model$residuals,
  res2 = eval_MM2$lm_model$residuals,
  res3 = eval_MM3$lm_model$residuals,
  res4 = eval_MM4$lm_model$residuals
)]
train$trait[train$trait == "T0"] <- "SUF"
train$trait[train$trait == "T4"] <- "TUF"
train$trait <- factor(train$trait, levels = c("TUF", "SUF", "T1", "T2", "T3"))

# Plotting function for residuals
plot_residuals <- function(train, res_col, subtitle, rmse, y_lim = c(-10, 10)) {
  ggplot(train, aes(x = trait, y = -get(res_col))) +
    ylab("residuals(Hpred-Href)") +
    ylim(y_lim) +
    xlab("forest type") +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', size = 0.5) +
    labs(subtitle = subtitle) +
    annotate("text", x = 1.2, y = y_lim[2], label = paste0("rmse=", rmse), size = 3) +
    theme_minimal()
}

# Create residual plots by forest type
MM1_trait <- plot_residuals(train, "res1", "Model 1", eval_MM1$rmse)
MM2_trait <- plot_residuals(train, "res2", "Model 2 (LHC)", eval_MM2$rmse)
MM3_trait <- plot_residuals(train, "res3", "Model 3 (species)", eval_MM3$rmse)
MM4_trait <- plot_residuals(train, "res4", "Model 4 (LHC+species)", eval_MM4$rmse)

# Create residual plots by habitat
plot_residuals_hab <- function(train, res_col, subtitle, rmse, y_lim = c(-10, 10)) {
  ggplot(train[!is.na(habitat)], aes(x = habitat, y = -get(res_col))) +
    ylab("residuals(Hpred-Href)") +
    ylim(y_lim) +
    xlab("habitat") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', size = 0.5) +
    labs(subtitle = subtitle) +
    annotate("text", x = 1.2, y = y_lim[2], label = paste0("rmse=", rmse), size = 3) +
    theme_minimal()
}

MM1_hab <- plot_residuals_hab(train, "res1", "Model 1", eval_MM1$rmse)
MM2_hab <- plot_residuals_hab(train, "res2", "Model 2 (LHC)", eval_MM2$rmse)
MM3_hab <- plot_residuals_hab(train, "res3", "Model 3 (species)", eval_MM3$rmse)
MM4_hab <- plot_residuals_hab(train, "res4", "Model 4 (LHC+species)", eval_MM4$rmse)

# Arrange and save the figure
fig7 <- ggarrange(MM1_hab, MM2_hab, MM3_hab, MM4_hab, MM1_trait, MM2_trait, MM3_trait, MM4_trait, ncol = 4, nrow = 2, common.legend = TRUE)

# Compute Basal Area (BA) and Quadratic Mean Diameter (QMD)
BA <- function(vec_dbh) {
  sum((vec_dbh / 2)^2 * pi / 10000) / 1.56
}

# Create data frame for BA and QMD
df <- data.frame(
  ID = as.numeric(tapply(train$square_125, train$square_125, unique)),
  trait = as.numeric(tapply(train$trait, train$square_125, unique)),
  dens = as.numeric(tapply(train$dbh.y, train$square_125, length)) / 1.56,
  BA = as.numeric(tapply(train$dbh.y, train$square_125, BA)),
  res1 = as.numeric(tapply(train$res1, train$square_125, mean)),
  res2 = as.numeric(tapply(train$res2, train$square_125, mean)),
  res3 = as.numeric(tapply(train$res3, train$square_125, mean)),
  res4 = as.numeric(tapply(train$res4, train$square_125, mean))
)

df$QMD <- sqrt((df$BA / df$dens) / pi) * 2
df$trait <- factor(df$trait, levels = c("TUF", "SUF", "T1", "T2", "T3"))

# Plot Residuals vs. BA and QMD
plot_residuals_BA <- function(df, res_col, subtitle) {
  ggplot(df, aes(x = BA, y = -get(res_col))) +
    geom_point(aes(color = trait)) +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(label.x = 18) +
    labs(x = expression(paste("Basal area (", m^2, ".", ha^-1, ")")), y = "residuals(Hpred-Href) in m") +
    ylim(-3, 3) +
    xlim(17, 30) +
    labs(subtitle = subtitle) +
    theme_minimal()
}

ResMM1_BA <- plot_residuals_BA(df, "res1", "Model 1")
ResMM2_BA <- plot_residuals_BA(df, "res2", "Model 2")
ResMM3_BA <- plot_residuals_BA(df, "res3", "Model 3")
ResMM4_BA <- plot_residuals_BA(df, "res4", "Model 4")


plot_residuals_QMD <- function(df, res_col, subtitle) {
  ggplot(df, aes(x = QMD, y = -get(res_col))) +
    geom_point(aes(color = trait)) +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(label.x = 0.30) +
    labs(x = "Quadratic Mean Diameter, in m", y = "residuals(Hpred-Href) in m") +
    ylim(-3, 3) +
    xlim(0.29, 0.43) +
    labs(subtitle = subtitle) +
    theme_minimal()
}


ResMM1_QMD <- plot_residuals_QMD(df, "res1", "Model 1")
ResMM2_QMD <- plot_residuals_QMD(df, "res2", "Model 2")
ResMM3_QMD <- plot_residuals_QMD(df, "res3", "Model 3")
ResMM4_QMD <- plot_residuals_QMD(df, "res4", "Model 4")

# Arrange and save the figure
fig8 <- ggarrange(ResMM1_BA, ResMM2_BA, ResMM3_BA, ResMM4_BA, ResMM1_QMD, ResMM2_QMD, ResMM3_QMD, ResMM4_QMD, ncol = 4, nrow = 2, common.legend = TRUE)



###################################################################################
## AGB - Figure 8 - Estimates of above-ground biomass in Mg.ha-1 per plot of 1.56 
## ha. The colors indicate the different treatments. 
###################################################################################
# Processes allometric and inventory data from various sources, corrects species 
# names, merges datasets, applies allometric models, calculates above-ground biomass (AGB), 
# and visualizes the results through plots

# Load inventory data
# Data from CIRAD plots (Dataset2), available for consultation with permission:
#   Derroire, G., Hérault, B., Rossi, V., Blanc, L., Gourlet- Fleury, S., & Schmitt, L. (2023). Paracou permanent plots. CIRAD Dataverse. https://doi.org/10.18167/DVN1/8G8AHY
data_CIRAD <- fread("Dataset2.csv")

# Import list of species with species-specific estimates into CIRAD database
not_sing <- names(table(data_allo$species_allo))[which(table(data_allo$species_allo) > 1)]

data_CIRAD$singleton <- TRUE
data_CIRAD$species_allo <- paste(data_CIRAD$genus, data_CIRAD$species, sep = ' ')

# Correct inconsistent naming of species in databases (Revised nomenclature in lcvp)
corrections <- list(
  "Agonandra silvatica" = "Agonandra sylvatica",
  "Alexa wachenheimii" = "Alexa wachenheimi",
  "Eschweilera grandiflora_form2" = "Eschweilera grandiflora",
  "Inga capitata_form2" = "Inga capitata",
  "Lecythis persistens subsp. aurantiaca" = "Lecythis persistens",
  "Licaria cannella" = "Licaria canella",
  "Parinari rodolphii" = "Parinari rodolphi"
)
data_CIRAD$species_allo <- as.factor(data_CIRAD$species_allo)
data_CIRAD$species_allo <- ifelse(data_CIRAD$species_allo %in% names(corrections), corrections[data_CIRAD$species_allo], data_CIRAD$species_allo)

# Update singleton status
data_CIRAD$singleton[data_CIRAD$species_allo %in% not_sing] <- FALSE
data_CIRAD$species_noSing <- data_CIRAD$species_allo
data_CIRAD$species_noSing[data_CIRAD$singleton == TRUE] <- "bin"  # Singletons are classified in the same class

# Import BRMS adjusted models
model4_MM_nsing <- readRDS("model4_MM_nsing_field.rds")
model1_MM <- readRDS("model1_MM_field.rds")

# Create dataset to apply allometric models
newdat <- data_CIRAD[family != "Arecaceae", .(idtree, dbh, species_noSing, CHM2015_25ppm_fo_comp, trait, wd, square_125)]
names(newdat) <- c("id", "dbh_allo", "species_noSing", "HC", "trait", "wd", "square_125")
newdat$pred1 <- fitted(model1_MM, newdata = newdat)[, 1]
newdat$pred4 <- fitted(model4_MM_nsing, newdata = newdat)[, 1]

# Import CC predictions of height per tree
files <- list.files(path = "P1to16_00006_0010", pattern = "trees_final.txt", full.names = TRUE)
cc <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
cc <- dcast(cc, id ~ ., value.var = "height", fun.aggregate = mean, na.rm = TRUE)
names(cc) <- c("id", paste0("H", 6:10))
cc[, predCC := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0("H", 6:10)]
newdat <- merge(newdat, cc[, .(id, predCC)], by = "id", all.x = TRUE)

# AGB calculation functions Chave et al. 2014
agb_eq_H <- function(wd, dbh, height) {
  agb <- (0.0673 * (wd * (dbh^2) * height)^0.976) * 0.001
  return(agb)
}

pantropical_Chave14 <- function(DBH, TS, CWD, PS) {
  E <- -0.1074988
  H <- exp(0.893 - E + 0.760 * log(DBH) - 0.0340 * (log(DBH))^2)
  return(H)
}

# E oleracea AGB calculation
AGB_pinot <- function(dbh) {
  AGB <- exp(-3.863 + 2.987 * log(dbh)) * 0.001
  return(AGB)
}

# Other palms AGB calculation
AGB_palm <- function(dbh) {
  AGB <- exp(-3.3488 + 2.7483 * log(dbh)) * exp(0.588^2 / 2) * 0.001
  return(AGB)
}

# AGB calculation for palms and non-palms
palm <- data_CIRAD[family == "Arecaceae", .(idtree, dbh, species_allo, CHM2015_25ppm_fo_comp, trait, wd, square_125)]
names(palm) <- c("id", "dbh_allo", "species_noSing", "HC", "trait", "wd", "square_125")
palm[, `:=` (pred1 = -99, pred4 = -99, predCC = -99, AGB0 = 0)]
palm[species_noSing != "Euterpe oleracea", AGB0 := AGB_palm(dbh_allo)]
palm[species_noSing == "Euterpe oleracea", AGB0 := AGB_pinot(dbh_allo)]
palm[, `:=` (AGB1 = AGB0, AGB4 = AGB0, AGBcc = AGB0)]

newdat[, AGB0 := agb_eq_H(wd, dbh_allo, pantropical_Chave14(dbh_allo, TS = 542, CWD = -103, PS = 41))]
newdat[, `:=` (AGB1 = agb_eq_H(wd, dbh_allo, pred1), AGB4 = agb_eq_H(wd, dbh_allo, pred4), AGBcc = agb_eq_H(wd, dbh_allo, predCC))]

fulldat <- rbind(newdat, palm)
fulldat <- na.omit(fulldat)
fulldat <- na.omit(fulldat)

# AGB summaries by trait
tapply(fulldat$AGB0, fulldat$trait, sum, na.rm = TRUE) / c(37.5, 18.75, 18.75, 18.75, 25)
tapply(fulldat$AGB1, fulldat$trait, sum, na.rm = TRUE) / c(37.5, 18.75, 18.75, 18.75, 25)
tapply(fulldat$AGB4, fulldat$trait, sum, na.rm = TRUE) / c(37.5, 18.75, 18.75, 18.75, 25)
tapply(fulldat$AGBcc, fulldat$trait, sum, na.rm = TRUE) / c(37.5, 18.75, 18.75, 18.75, 25)

# Filter out trees missing square_125
fulldat <- fulldat[!is.na(square_125), ]
fulldat$treatment <- factor(ifelse(fulldat$trait == "T4", "TUF", ifelse(fulldat$trait == "T0", "SUF", fulldat$trait)), levels = c("TUF", "SUF", "T1", "T2", "T3"))

# Plot
paleta<-c("#6a3d9a","#33a02c","#1f78b4","#ff7f00","#e31a1c")
figa<-ggplot(byplot, aes(x=AGBcc, y=AGB0, color=treatment)) + 
  geom_point( ) +
  geom_abline(intercept = 0, slope = 1,  color = "black")+
  ylim(250, 550) +
  xlim(250, 550) +
  ylab("Pantropical model") + 
  xlab("Canopy constructor") +
  scale_color_manual(values = paleta,aesthetics = c("colour", "fill"))+
  theme_minimal()

figb<-ggplot(byplot, aes(x=AGBcc, y=AGB1, color=treatment)) + 
  geom_point( ) +
  geom_abline(intercept = 0, slope = 1,  
              color = "black")+
  ylim(250, 550) +
  xlim(250, 550) +
  ylab("Model 1 (local)") + 
  xlab("Canopy constructor") +
  scale_color_manual(values = paleta,aesthetics = c("colour", "fill"))+
  theme_minimal()

figc<-ggplot(byplot, aes(x=AGBcc, y=AGB4, color=treatment)) + 
  geom_point( ) +
  geom_abline(intercept = 0, slope = 1,  
              color = "black")+
  ylim(250, 550) +
  xlim(250, 550) +
  ylab("Model 4 (species + LCH") + 
  xlab("Canopy constructor") +
  scale_color_manual(values = paleta,aesthetics = c("colour", "fill"))+
  theme_minimal()

figd<-ggplot(byplot, aes(x=AGB4, y=AGB0, color=treatment)) + 
  geom_point( ) +
  geom_abline(intercept = 0, slope = 1,  
              color = "black") +
  ylim(250, 550) +
  xlim(250, 550) +
  ylab("Model 4 (species + LCH") + 
  xlab("Pantropical model") +
  scale_color_manual(values = paleta,aesthetics = c("colour", "fill"))+
  theme_minimal()

fig8<-ggarrange(figa, figb, figc, figd, ncol = 2, nrow = 2, common.legend = T)

                      