library(tidyverse)
library(rjags)

# Models 
## Currently working with the full.mix model
#bef.combined.model = read_file("Working scripts/Models/bef_combined.txt")
#zero.inflated.model = read_file("Working scripts/Models/zero_inflated.txt")
full_model = read_file("Data/Models/updated_full.txt")
model = read_file("Data/Models/cscu_model2") ## Model two gets rid of one of the varience parameters for detection probability based on site

# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=2*10^5; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)



# Loading in all Jags Data from other scripts -------

#### Load in Isotope Data ----------
load("Data/VaryingIsotopesData/LateData/IsotopeArray_late.RData")
# Number of observations
n.obs <- dim(arr)[2]
# Number of isotopes
n.iso <- dim(arr)[3]

#### Load in Mass Data ----------

load("Data/VaryingIsotopesData/LateData/n_obs_mass_late.RData")
load("Data/VaryingIsotopesData/LateData/observed_lengths_late.RData")

#### Load in Effort Data --------
load("Data/VaryingIsotopesData/LateData/late_effortdata.RData")


#### Load in the Catch Data --------------
load("Data/VaryingIsotopesData/LateData/mean_catch_array_late.RData") ## New model

# Number of sites, years, and species


## CMax for innits is the maximum number for any site/species combo

Cmax = apply(mean_catch_array, c(1,4), max,na.rm=T) ## Set up innits


## Set up variables for JAGS data
sites =  dim(mean_catch_array)[1]
year = dim(mean_catch_array)[3]
replicates = dim(mean_catch_array)[2]
species = dim(mean_catch_array)[4]

#### Load in covariates -------------------


## New model covariates
load("Data/VaryingIsotopesData/LateData/ice_off.RData") # Ice off new model
load("Data/VaryingIsotopesData/LateData/hab_cov.RData") # Habitat new model
load("Data/VaryingIsotopesData/LateData/shoreline_length.RData") # Length new model
load("Data/VaryingIsotopesData/site_dist.RData") # Site dissimirity matrix based on adjacency and habitat

# Jags Data
alpha <- 1
# Jags Data
beta <- 1



# Specify parameters to monitor ------------------------------------------------





### RJAGS




## Prior setup
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 



jags_data = list(sites = sites,
                 year = year,

                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 shoreline_length = as.numeric(shoreline_length$shoreline_std),
                 hab = scale(as.numeric(hab$sub)) %>% as.numeric,
                 num_hab = 2, 
                 wood = as.numeric(hab$wood), 
                 ice_off = ice_off,
                 site_dist = dist_vec,
                 ## mass components
                 mass_observed_lake = observed_lengths.late,
                 n_obs_mass = n_obs_mass.late,
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
n_species = species


inits = function() list(N = Cmax, 
                        avg_mass  = rnorm(n_species),
                        mu =array(data = rnorm(n_species * n.iso),
                                  dim = c(n_species, n.iso)))



## Size data -------------------------------------------------
# Specify parameters to monitor ------------------------------------------------
jags_parameters = c("lambda", "p", "N_total","avg_mass", "mu", "Sigma2")
# Compile and run the model ----------------------------------------------------


jags_model <- jags.model(textConnection(model), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)


jags_output.late <- coda.samples(jags_model, 
                                 variable.names = jags_parameters, 
                                 n.iter = parms$n.iter, 
                                 thin = parms$n.thin,
                                 inits = inits)



save(file = "Data/VaryingIsotopesData/LateData/LateJagsOutput.RData", jags_output.late)
load("Data/VaryingIsotopesData/LateData/LateJagsOutput.RData")
