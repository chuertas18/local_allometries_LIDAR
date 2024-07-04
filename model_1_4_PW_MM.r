rm(list=ls())
##################################################################################################################
# Paracou local allometry model process
# Created by Claudia Huertas
# Adapted by previous code: Mélaine Aubry-kientz and Fabian Fischer
# This script allows the local adjustment of H:DBH allometries, using bayesian method with a power function.
##################################################################################################################

##################################################################################################################
# 1/ Libraries
##################################################################################################################
# list.of.packages <- c("brms")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

# library(data.table) # Enables more efficient and faster dataframe reading and writing
library(brms) # Connection to STAN for Bayesian modeling
library(Rcpp)


##################################################################################################################
# 2/ Databases
##################################################################################################################
# Data from CIRAD plots (Dataset2), available for consultation with permission:
#   Derroire, G., Hérault, B., Rossi, V., Blanc, L., Gourlet- Fleury, S., & Schmitt, L. (2023). Paracou permanent plots. CIRAD Dataverse. https://doi.org/10.18167/DVN1/8G8AHY
bd_paracou_HC<-read.csv("Dataset2.csv")
colnames(bd_paracou_HC)[which(names(bd_paracou_HC) == "CHM2015_25ppm_fo_comp")]  <- 'HC'
bd_paracou_HC$species_allo<-paste(bd_paracou_HC$genus,bd_paracou_HC$species)

my_read_delim <- function(path){
  # readr::read_delim(path, "\t", escape_double = FALSE, trim_ws = TRUE)
  read.table(path, sep="\t", header=TRUE)
}
trees.final.PL = do.call(rbind,lapply(list.files(path = "/Data/P1to16_00006_0010/P1to16.MM",pattern = "final.txt", full.names = T), my_read_delim))


# trees.final.PL[, CAparam := ifelse(paramID <= 5,"original","convHull")]
# Check if you are reducing to the mean value.
trees.final.PL = trees.final.PL[which(trees.final.PL$paramID > 5 & trees.final.PL$dbh>=0.20),]


# Union of the two databases
data_cc_cirad<-merge(trees.final.PL,bd_paracou_HC[,c("idtree","dbh", "wd",  "HC","species_allo",
                                                     "plot","habitat","trait","square_125")],
                     by.x = 'id', by.y = 'idtree')

# Determination of singleton species
singletons=names(table(data_cc_cirad$species_allo))[which(table(data_cc_cirad$species_allo)==1)]
data_cc_cirad$singleton=F
data_cc_cirad$singleton[which(data_cc_cirad$species_allo %in% singletons)]<-T
data_cc_cirad$species_noSing<-data_cc_cirad$species_allo
data_cc_cirad$species_noSing[data_cc_cirad$singleton==T]<-"bin" # Singletons are classified in the same class

# Filter with diameter greater than or equal to 20 cm
data_cc_cirad<-data_cc_cirad[which(data_cc_cirad$dbh.y>=20),]

##################################################################################################################
# 3/ Run Bayesian model
##################################################################################################################
###################################################################################################
# Model 1/ Single allometry for the entire site 
###################################################################################################
# Model 1 - Power law
model1_pw_nsing_CC = brm(
  bf(height~alpha*dbh.y^beta, ## Power law
     alpha ~ 1,
     beta ~ 1,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 1000, 
  warmup = 500, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)

# Model 1 - Michaelis Menten
model1_MM_nsing_CC = brm(
  bf(height~(alpha*dbh.y)/(beta+dbh.y), ## 
     alpha ~ 1,
     beta ~ 1,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 1000, 
  warmup = 500, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# ###################################################################################################
## 4/ Second model (model2) includes canopy height in a local neighborhood (HC)
# ###################################################################################################
# Model 2 - Power law
model2_pw_nsing_CC = brm(
  bf(height~alpha*dbh.y^beta, ## Power law
     alpha ~ 1 + HC,
     beta ~ 1,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 2000, 
  warmup = 1000, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


# Model 2 - Michaelis Menten
model2_MM_nsing_CC = brm(
  bf(height~(alpha*dbh.y)/(beta+dbh.y), 
     alpha ~ 1 + HC,
     beta ~ 1,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha), 
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  # parameters of the MCMC algorithm
  iter = 2000, 
  warmup = 1000, 
  chains = 4, 
  cores = 7,
  control = list(adapt_delta = 0.8), 
  seed = 25,
  silent = FALSE
)


###################################################################################################
# 5/ model 3 Species
###################################################################################################
# Model 3- Power law
model3_pw_nsing_CC = brm(
  bf(height~alpha*dbh.y^beta,
     alpha ~ 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  # control = list(adapt_delta = 0.95,
  #                max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
  # backend = "cmdstanr",
  # threads = threading(4)
)



# Model 3 Michaelis-Menten
model3_MM_nsing_CC = brm(
  bf(height~(alpha*dbh.y)/(beta+dbh.y),
     alpha ~ 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
            prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  cores = 7,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE,
  backend = "cmdstanr",
  threads = threading(4)
)


###################################################################################################
# 4/ model 4 Local canopy height + species
###################################################################################################
# Model 4- Power law
model4_pw_nsing_CC= brm(
  bf(height~alpha*dbh.y^beta,
     alpha ~ HC + 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  # in case there are problems with fitting the model
  # control = list(adapt_delta = 0.95),
  seed = 25,
  silent = FALSE
  # backend = "cmdstanr",
  # threads = threading(4)
)


# Model 4 Michaelis-Menten
model4_MM_nsing_CC= brm(
  bf(height~(alpha*dbh.y)/(beta+dbh.y),
     alpha ~ HC + 1|species_noSing,
     beta ~ 1|species_noSing,
     nl = TRUE),
  data = data_cc_cirad,
  family = "gaussian",
  prior = c(
    prior(gamma(1.5, 0.01), lb = 0, nlpar = alpha),
    prior(gamma(1.5, 0.01), lb = 0, nlpar = beta)
  ),
  iter = 3000,
  warmup = 1500,
  chains = 4,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12), # adapt_delta can be set to 0.99 and max_treedepth = 12,
  seed = 25,
  silent = FALSE
)
