 library(tidyverse)
library(rjags)

# Models 
## Currently working with the model specifically for the 2000
#bef.combined.model = read_file("Working scripts/Models/bef_combined.txt")
#zero.inflated.model = read_file("Working scripts/Models/zero_inflated.txt")
full_model = read_file("Data/Models/Model_2000.txt")

# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)








# Loading in all Jags Dat from other scripts -------

#### Load in Isotope Data ----------
load("Data/VaryingIsotopesData/LateData/IsotopeArray_late.RData") ## Use this for the unvaried isotopes
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
load("Data/VaryingIsotopesData/PreData/mean_catch_array_pre.RData") ## New model

## Set up variables for JAGS data

## CMax for innits is the maximum number for any site/species combo
Cmax = apply(mean_catch_array, c(1,3), max,na.rm=T) ## Set up innits


# Number of sites, years, and species

## Set up variables for JAGS data
sites =  dim(mean_catch_array[,,1])[1]

replicates = (mean_catch_array%>% dim())[2]
species = (mean_catch_array %>% dim())[3]


#### Load in covariates -------------------



load("Data/VaryingIsotopesData/PreData/ice_off_pre.RData")

load("Data/VaryingIsotopesData/PreData/habs.RData")


# Jags Data
alpha <- 1
# Jags Data
beta <- 1



# Specify parameters to monitor ------------------------------------------------
## Prior setup
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 

jags_data = list(sites = sites,
                 
                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 
                 hab = as.numeric(habs[,3]),
                 num_hab = 2, 
                 wood = as.numeric(habs[,4]), 
                 ice_off = ice_off[1],
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
n_species = species


inits = function() list(N = Cmax, 
                        avg_mass  = rnorm(n_species),
                        mu =array(data = rnorm(n_species * n.iso),
                                  dim = c(n_species, n.iso)))


## Size data -------------------------------------------------
# Specify parameters to monitor ------------------------------------------------
jags_parameters = c("lambda", "p", "N","avg_mass", "mu", "Sigma2")
# Compile and run the model ----------------------------------------------------


jags_model <- jags.model(textConnection(full_model), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)


jags_output.pre.un <- coda.samples(jags_model, 
                                     variable.names = jags_parameters, 
                                     n.iter = parms$n.iter, 
                                     thin = parms$n.thin,
                                     inits = inits)
save(file = "Data/UnvariedIsotopeData/PreData/PreJagsOutput.RData", jags_output.pre.un)
