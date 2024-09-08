library(tidyverse)
library(rjags)

# Models 
## Currently working with the zero.aggregate model
#bef.combined.model = read_file("Working scripts/Models/bef_combined.txt")
#zero.inflated.model = read_file("Working scripts/Models/zero_inflated.txt")

zero.inflated.model = read_file("Data/Models/zero_inflated.txt")


# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)








# Loading in all Jags Dat from other scripts -------

#### Load in Isotope Data ----------
load("Data/EarlyData/IsotopeArray_early.RData")
# Number of observations
n.obs <- dim(arr)[2]
# Number of isotopes
n.iso <- dim(arr)[3]

#### Load in Mass Data ----------

load("Data/EarlyData/n_obs_mass_early.RData")
load("Data/EarlyData/observed_lengths_early.RData")

#### Load in Effort Data --------

## Load in effort data
load("Data/EarlyData/early_effortdata.RData")


#### Load in the Catch Data --------------
load("Data/EarlyData/Early_arraydata.RData")

# Number of sites, years, and species
n_sites <- dim(species_cpue.data)[2] -2
n_years <- dim(effort)[2]
n_species <- dim(array_data)[3]

#### Load in covariates -------------------
load("Data/EarlyData/habs_cov.RData")
load("Data/EarlyData/covariate_temp.RData")
load("Data/EarlyData/z_covariate.RData")
# Jags Data
alpha <- 1
# Jags Data
beta <- 1



# Specify parameters to monitor ------------------------------------------------

jags_parameters = c("avg_mass","N", "mu", "Sigma2")
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 ## Prior setup




## Prior setup
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 



# Data for the JAGS model
data_list <- list(
  # Values for for loops
  n_sites = n_sites,
  n_years = n_years,
  n_species = n_species, 
  
  # Catch and Effort
  y = array_data, ## observed data in the practice one, i've changed it here
  effort = effort,
  
  # Covariates
  x = covariate_T, # temp for early model
  z_covariate =  matrix(z_covariate, nrow = n_sites, ncol = n_years),
  habitat = covariate_hab,
  
  
  # Mass data
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

inits <- replicate(parms$n.chains,
                   list(
                     mu =array(data = rnorm(n_species * n.iso),
                               dim = c(n_species, n.iso)),## Isotope model innit
                     tau = array(data = array(0, 
                                              dim = c(n_species, 
                                                      n.iso, n.iso)),
                                 dim = c(n_species, n.iso, n.iso)), ## Isotope model innit
                     avg_mass = rnorm(n_species), ## Mass innit
                     beta0 = rep(0, n_species),
                     beta1 = rep(0, n_species),
                     alpha0 = rep(0, n_species), 
                     alpha1 = rep(0, n_species), 
                     sigma_u = runif(1,0,1),
                     sigma_v = runif(1,0,1)),
                   simplify = FALSE)

## Size data -------------------------------------------------
# Specify parameters to monitor ------------------------------------------------
jags_parameters = c("avg_mass","N", "mu", "Sigma2")
# Compile and run the model ----------------------------------------------------
jags_model <- jags.model(textConnection(zero.inflated.model), 
                         data = data_list, 
                         n.chains = parms$n.chains)
jags_output.early <- coda.samples(jags_model, 
                            variable.names = jags_parameters, 
                            n.iter = parms$n.iter, 
                            thin = parms$n.thin,
                            inits = inits)


save(file = "Data/EarlyData/EarlyJagsOutput.RData", jags_output.early)
