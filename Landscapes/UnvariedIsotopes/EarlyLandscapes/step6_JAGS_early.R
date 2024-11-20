library(tidyverse)
library(rjags)

# Models 
## Currently working with the updated_full model
#bef.combined.model = read_file("Working scripts/Models/bef_combined.txt")
#zero.inflated.model = read_file("Working scripts/Models/zero_inflated.txt")
full_model = read_file("Data/Models/updated_full.txt")


# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)








# Loading in all Jags Dat from other scripts -------

#### Load in Isotope Data ----------
load("Data/VaryingIsotopesData/LateData/IsotopeArray_late.RData") ## This is the late isotope data that im using across unvaried landscapes
# Number of observations
n.obs <- dim(arr)[2]
# Number of isotopes
n.iso <- dim(arr)[3]

#### Load in Mass Data ----------

load("Data/VaryingIsotopesData/EarlyData/n_obs_mass_early.RData")
load("Data/VaryingIsotopesData/EarlyData/observed_lengths_early.RData")

#### Load in Effort Data --------

## Load in effort data
load("Data/VaryingIsotopesData/EarlyData/early_effortdata.RData")


#### Load in the Catch Data --------------
load("Data/VaryingIsotopesData/EarlyData/mean_catch_array_early.RData") ## New model

## Set up variables for JAGS data

Cmax = apply(mean_catch_array, c(1,4), max,na.rm=T) ## Set up innits

sites =  dim(mean_catch_array[,,1,1])[1]
year = (mean_catch_array[,,,1] %>% dim())[3]
replicates = (mean_catch_array[,,,1] %>% dim())[2]
species = (mean_catch_array %>% dim())[4]


#### Load in covariates -------------------
load("Data/VaryingIsotopesData/EarlyData/habs_cov.RData")
load("Data/VaryingIsotopesData/EarlyData/covariate_temp.RData")
load("Data/VaryingIsotopesData/EarlyData/z_covariate.RData")


load("Data/VaryingIsotopesData/EarlyData/ice_off_early.RData")
load("Data/VaryingIsotopesData/EarlyData/shoreline_length_early.RData")
load("Data/VaryingIsotopesData/EarlyData/hab_cov.RData")



# Jags Data
alpha <- 1
# Jags Data
beta <- 1



# Specify parameters to monitor ------------------------------------------------
# Specify parameters to monitor ------------------------------------------------





### RJAGS




## Prior setup
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 

jags_data = list(sites = sites,
                 year = year,
                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 site_proportion = shoreline_length[,1:2] %>% 
                   mutate(`2` = replace_na(`2`, 1)),
                 hab = as.numeric(hab$sub),
                 num_hab = 2, 
                 wood = as.numeric(hab$mean_w), 
                 ice_off = ice_off,
                 ## mass components
                 mass_observed_lake = observed_lengths.early,
                 n_obs_mass = n_obs_mass.early,
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


jags_output.early.un <- coda.samples(jags_model, 
                                 variable.names = jags_parameters, 
                                 n.iter = parms$n.iter, 
                                 thin = parms$n.thin,
                                 inits = inits)



save(file = "Data/UnvariedIsotopeData/EarlyData/EarlyJagsOutput.RData", jags_output.early.un)
