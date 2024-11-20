library(rjags)
library(tidyverse)

setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape/")
source("Landscapes/VaryingIsotopes/PreLandscapes/step1_LML_source_pre.R")

## I'm a little worried that the multiple observation issue is going to catch up again here. But, for now, I'm simplifying this by just taking average catch and average effort for each site in each year. Not putting in each separately.


sample %>% filter(YEAR == 2000, MONTH < 7, WATER == "LML", GEAR == "BEF") %>%
  group_by(SITE_N, GEAR_CODE) %>%
  summarize(count = n()) %>%
  print(n = 100)



sample %>% filter(YEAR == 2000, MONTH < 7, WATER == "LML", GEAR == "BEF") %>%
  ggplot(aes(x = DAY_N, y = SITE_N)) + 
  geom_point()

step_a = LML.data %>%
  filter(YEAR >= year_min & YEAR != 2002 & YEAR <= year_max, ## here it details the min and max
         SITE != "BEF.LML.NA") %>%
  separate(SPECIES, into = c("SPECIES", "AGE")) %>%
  select(YSAMP_N, DAY_N, YEAR, SEASON, WATER, SITE, SPECIES,
         FISH_N, WEIGHT, LENGTH,  EFFORT) %>%
  group_by(WATER, DAY_N, YEAR, SITE, SPECIES, EFFORT) %>%
  count() %>% ## Abundance per year, site, species
  mutate(CPUE_seconds = n / EFFORT) %>%
  ungroup() %>%
  #mutate(CPUE_seconds = normalize(n/EFFORT)) %>%
  complete(., nesting(WATER, YEAR, DAY_N, SITE, EFFORT), SPECIES) %>%
  replace_na(list(CPUE_seconds = 0, n = 0)) %>% 
  #group_by(YEAR, SPECIES, SITE) %>%
  #summarize(mean_catch = mean(n), mean_effort = mean(EFFORT)) %>%
  ungroup() %>%
  filter(SPECIES %in% species) %>%
  rename("mean_catch" = "n")






#save(file = "Data/PreData/Pre_arraydata.RData", array_data)


### New model data setup --------------------------------------------


## Load in the data that says which sites to pair
rep_group = read.csv("Data/rep_groups.csv") ## Round values for modifying p

## Calculate what proportion of a shoreline is sampled during each replicate
shoreline_length = rep_group %>%
  select(shoreline_length, rep_group_simple) %>%
  na.omit() %>% group_by(rep_group_simple) %>%
  mutate(count = c(1:length(shoreline_length))) %>% ## index for pivoting wider
  pivot_wider(names_from = count, values_from = shoreline_length) %>% ## Pivot it so reps are columns
  ungroup() %>%
  select(-rep_group_simple) %>% ## remove columns not needed in RJAGS
  mutate(total_length = rowSums(.)) %>% ## Total length of shoreline for calculating proportion
  mutate(`1` = `1` / total_length,  ## Calculate the proportions for each replicate for each site
         `2` = `2` / total_length)




## Abundance Array -----------------------------------
## An intermediate step for creating the array
df = step_a %>% 
  #left_join(rep_group, by = c("SITE" = "site")) %>%
  mutate(mean_catch = round(mean_catch, digits = 0)) %>%
  select(YEAR, SPECIES, SITE, mean_catch) %>%
  na.omit() %>%
  
  group_by(YEAR,SPECIES, SITE) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  ungroup() 


# Get unique values for each dimension

species = unique(df$SPECIES)
#rep_groups <- sort(unique(df$rep_group_simple))
rep_nums <- sort(unique(df$rep_num))
sites = sort(unique(df$SITE))

# Create an empty 4D array
mean_catch_array <- array(
  NA, 
  dim = c(length(sites), length(rep_nums),  length(species)),
  dimnames = list(sites = sites, rep_num = rep_nums,  SPECIES = species)
)

# Fill the array
for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  mean_catch_array[
    as.character(row$SITE),
    as.character(row$rep_num),

    as.character(row$SPECIES)
  ] <- row$mean_catch
}

save(mean_catch_array, file = "Data/VaryingIsotopesData/LateData/mean_catch_array_pre.RData")


## Prior setup


# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)


# Loading in all Jags Data from other scripts -------

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


# Number of sites, years, and species


## CMax for innits is the maximum number for any site/species combo

Cmax = apply(mean_catch_array, c(1,3), max,na.rm=T) ## Set up innits


## Set up variables for JAGS data
sites =  dim(mean_catch_array[,,1])[1]

replicates = (mean_catch_array%>% dim())[2]
species = (mean_catch_array %>% dim())[3]

#### Load in covariates -------------------


## New model covariates
load("Data/VaryingIsotopesData/LateData/ice_off.RData") # Ice off new model
load("Data/VaryingIsotopesData/LateData/hab_cov.RData") # Habitat new model
load("Data/VaryingIsotopesData/LateData/shoreline_length.RData") # Length new model
load("Data/VaryingIsotopesData/PreData/habs.RData")



## Priors
priors=list(); priors$R=1*diag(2); priors$k=2; priors$tau.mu=1.0E-3 


# Get unique values for each dimension


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

model.2000 = read_file("Data/Models/Model_2000.txt")
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



