library(tidyverse)
library(rjags)

# Models 
## Currently working with the model built specifically for the year 2000 
#bef.combined.model = read_file("Working scripts/Models/bef_combined.txt")
#zero.inflated.model = read_file("Working scripts/Models/zero_inflated.txt")


model.2000 = read_file("Data/Models/cscu_model2")

# options for running jags
parms <- list(); parms$n.iter=3*10^5; parms$n.burnin=1*10^4; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)
# Loading in all Jags Dat from other scripts -------

#### Load in Isotope Data ----------
load("Data/VaryingIsotopesData/EarlyData/IsotopeArray_early.RData")
# Number of observations
n.obs <- dim(arr)[2]
# Number of isotopes
n.iso <- dim(arr)[3]

#### Load in Mass Data ----------

load("Data/VaryingIsotopesData/PreData/n_obs_mass_pre.RData")
load("Data/VaryingIsotopesData/PreData/observed_lengths_pre.RData")

#### Load in Effort Data --------

## Load in effort data
load("Data/VaryingIsotopesData/PreData/pre_effortdata.RData")


#### Load in the Catch Data --------------
load(file = "Data/VaryingIsotopesData/PreData/mean_catch_array_pre.RData")


## CMax for innits is the maximum number for any site/species combo
Cmax = apply(mean_catch_array, c(1,4), max,na.rm=T) ## Set up innits


# Number of sites, years, and species




## Set up variables for JAGS data

## Set up variables for JAGS data
sites =  dim(mean_catch_array)[1]
years = dim(mean_catch_array)[3] ## I messed around with labelling here. The years are actually the temporal replicates within the year 2000, and the reps are just 1. This makes it easier when adjustin the model to work across time periods...
replicates = (mean_catch_array%>% dim())[2]
species = (mean_catch_array %>% dim())[4]



#### Load in covariates -------------------
load("Data/VaryingIsotopesData/PreData/habs.RData")
load("Data/VaryingIsotopesData/PreData/ice_off_pre.RData")


load("Data/VaryingIsotopesData/LateData/shoreline_length.RData")
load("Data/VaryingIsotopesData/LateData/hab_cov.RData")


# Jags Data
alpha <- 1
# Jags Data
beta <- 1



# Specify parameters to monitor ------------------------------------------------




priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 ## Prior setup




## Prior setup
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 



# Data for the JAGS model
jags_data = list(sites = sites,
                 
                 replicates = replicates, 
                 species = species ,
                 year = years,
                 C = mean_catch_array, 
                 shoreline_length =(shoreline_length$shoreline_length %>% as.numeric()),
                 mean_length = mean((shoreline_length$shoreline_length)),
                 hab = hab$sub,
                 num_hab = 2, 
                 wood = hab$wood, 
                 ice_off = rep(ice_off[1], length = years),
                 ## mass components
                 mass_observed_lake = observed_lengths.pre,
                 n_obs_mass = n_obs_mass.pre,
                 # Isotope data
                 "n.obs" = n.obs,
                 "n.iso" = n.iso,
                 "Y" = arr,
                 # Priors
                 "R" = priors$R,
                 "k" = priors$k, 
                 "tau.mu" = priors$tau.mu
)


## Innits
#With size and abundance

inits = function() list(N = Cmax, 
                        avg_mass  = rnorm(species),
                        mu =array(data = rnorm(species * n.iso),
                                  dim = c(species, n.iso)))


## Size data -------------------------------------------------
# Specify parameters to monitor ------------------------------------------------
jags_parameters = c("lambda", "p", "N_total","avg_mass", "mu", "Sigma2")
# Compile and run the model ----------------------------------------------------

jags_model <- jags.model(textConnection(model.2000), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)


jags_output.pre <- coda.samples(jags_model, 
                                variable.names = jags_parameters, 
                                n.iter = parms$n.iter, 
                                thin = parms$n.thin,
                                inits = inits)

save(jags_output.pre, file = "Data/VaryingIsotopesData/PreData/PreJagsOutput.RData")

