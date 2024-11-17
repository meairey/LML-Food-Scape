## Set parameters for the model run
year_max = 2022
year_min = 2019

source("Landscapes/VaryingIsotopes/LateLandscapes/step1_LML_source_L.R")
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
  group_by(YEAR, SPECIES, SITE) %>%
  summarize(mean_catch = mean(n), mean_effort = mean(EFFORT)) %>%
  ungroup() %>%
  filter(SPECIES %in% species)




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
  left_join(rep_group, by = c("SITE" = "site")) %>%
  mutate(mean_catch = round(mean_catch, digits = 0)) %>%
  select(YEAR, SPECIES, rep_group_simple,  mean_catch) %>%
  na.omit()%>%
  
  group_by(YEAR,SPECIES, rep_group_simple) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  ungroup() 

# Get unique values for each dimension
years <- sort(unique(df$YEAR))
species <- unique(df$SPECIES)
rep_groups <- sort(unique(df$rep_group_simple))
rep_nums <- sort(unique(df$rep_num))

# Create an empty 4D array
mean_catch_array <- array(
  NA, 
  dim = c(length(rep_groups), length(rep_nums), length(years), length(species)),
  dimnames = list(rep_group = rep_groups, rep_num = rep_nums, YEAR = years, SPECIES = species)
)

# Fill the array
for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  mean_catch_array[
    as.character(row$rep_group_simple),
    as.character(row$rep_num),
    as.character(row$YEAR),
    as.character(row$SPECIES)
  ] <- row$mean_catch
}

# Print the array to check
print(mean_catch_array)

Cmax = apply(mean_catch_array, c(1,4), max,na.rm=T)
Cmax


## Set up variables for JAGS data
sites =  dim(mean_catch_array[,,1,1])[1]
year = (mean_catch_array[,,,1] %>% dim())[3]
replicates = (mean_catch_array[,,,1] %>% dim())[2]
species = (mean_catch_array %>% dim())[4]

## Habitat covariate 


hab = rep_group %>%
  select(rep_group_simple, hab) %>%
  na.omit() %>%
  mutate(sub = substr(hab, 1, 1),
         wood = substr(hab, 2, 2)) %>%
  mutate(wood = as.numeric(as.factor(wood)),
         sub  = as.numeric(as.factor(sub))) %>%
  select(rep_group_simple, sub, wood) %>%
  group_by(rep_group_simple, sub) %>%
  summarize(mean_w = mean(wood)) %>%
  unique()


ice_off = sample %>% left_join(read.csv("Data/ICE_OUT.csv"),by = c("YEAR" = "Year")) %>%
  filter(YEAR > 1997,
         GEAR == "BEF",
         GEAR_CODE == "NAF", 
         WATER == "LML", 
         MONTH < 7) %>%
  select(WATER, DATE_COL, MONTH, DAY_N, YEAR, Ice.out.date) %>%
  unique() %>%
  mutate(dif = DAY_N - Ice.out.date) %>%
  group_by(YEAR) %>%
  summarize(avg.day.ice = mean(dif) ) %>% 
  print(n = 100)
ice_off[24,2] = 30
ice_off[24,1] = 2022


ice_off =as.numeric((ice_off %>% filter(YEAR %in% c(year_min:year_max)))$avg.day.ice)




### RJAGS

# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

inits = function() list(N = Cmax)

jags_data = list(sites = sites,
                 year = year,
                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 site_proportion = shoreline_length[,1:2] ,
                 hab = as.numeric(hab$sub),
                 num_hab = 2, 
                 wood = as.numeric(hab$mean_w), 
                 ice_off = ice_off)

## Which parameters to monitor
jags_parameters= c("lambda", "p", "N","Ntotal")

## Which model to use
model = read_file("Data/Models/Nmix_species.txt")

jags_model <- jags.model(textConnection(model), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)
jags_output.late <- coda.samples(jags_model, 
                                 variable.names = jags_parameters, 
                                 n.iter = parms$n.iter, 
                                 thin = parms$n.thin,
                                 inits = inits)


chains <- as.mcmc.list(jags_output.late) ## set chains
chains_combined <- gtable_combine(chains) ## combine chains

dat = chains_combined %>%as.matrix() %>%  as.data.frame() #%>%

dat %>%
  ggplot(aes(x = Ntotal)) + 
  geom_histogram() + 
  scale_x_log10()

n.posts = 100

abund.dat =  chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("N", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = sum(value)) %>%

  left_join(legend, by = "group")

abund.dat %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot() +
  scale_y_log10()

legend


detection_probability = chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("p", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = mean(value)) %>%
  
  left_join(legend, by = "group") 


detection_probability %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot


### Adding in mass
## Mass components 



# Mass data
inits = function() list(N = Cmax, 
                        avg_mass  = rnorm(n_species))

tau = array(data = array(0, 
                         dim = c(n_species, 
                                 n.iso, n.iso)),
            dim = c(n_species, n.iso, n.iso)), ## Isotope model innit


jags_data = list(sites = sites,
                 year = year,
                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 site_proportion = shoreline_length[,1:2] ,
                 hab = as.numeric(hab$sub),
                 num_hab = 2, 
                 wood = as.numeric(hab$mean_w), 
                 ice_off = ice_off,
                 ## mass components
                 mass_observed_lake = observed_lengths.late,
                 n_obs_mass = n_obs_mass.late)

## Which parameters to monitor
jags_parameters= c("lambda", "p", "N","avg_mass")

## Which model to use
model = read_file("Data/Models/Nmix_species.txt")
new_model = read_file("Data/Models/nmix_weight.txt")
jags_model <- jags.model(textConnection(new_model), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)
jags_output.late <- coda.samples(jags_model, 
                                 variable.names = jags_parameters, 
                                 n.iter = parms$n.iter, 
                                 thin = parms$n.thin,
                                 inits = inits)
chains <- as.mcmc.list(jags_output.late) ## set chains
chains_combined <- gtable_combine(chains) ## combine chains

dat = chains_combined %>%as.matrix() %>%  as.data.frame() #%>%

dat %>%
  ggplot(aes(x = Ntotal)) + 
  geom_histogram() + 
  scale_x_log10()

n.posts = 100

abund.dat =  chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("N", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = sum(value)) %>%
  
  left_join(legend, by = "group")

abund.dat %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot() +
  scale_y_log10()

legend


detection_probability = chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("p", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","year","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = mean(value)) %>%
  
  left_join(legend, by = "group") 


detection_probability %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot()



## Clean length estimates
mass.dat = dat %>%
  select(contains("avg_mass", ignore.case = TRUE))  %>%
  as.data.frame() %>%
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]),
               names_to = "metric", values_to = "mass.avg") %>%
  separate(metric, into = c("metric1","metric"), sep = "_") %>%
  mutate(group = parse_number(metric)) %>%
  select(-metric1, -metric) %>%
    left_join(legend)

# Figure for looking at the mass data
mass.dat %>%


  ggplot(aes(x = species, y = exp(mass.avg)))  +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_log10()+

  ylab("Avg. Mass (g)") +
  xlab("Species") +
  labs(col = "Species")



### Adding in the isotopes

inits = function() list(N = Cmax, 
                        avg_mass  = rnorm(n_species))

## Which parameters to monitor
jags_parameters= c("lambda", "p", "N","avg_mass", "mu", "Sigma2")

## Which model to use

full_model = read_file("Data/Models/updated_full.txt")



jags_data = list(sites = sites,
                 year = year,
                 replicates = replicates, 
                 species = species ,
                 C = mean_catch_array, 
                 site_proportion = shoreline_length[,1:2] ,
                 hab = as.numeric(hab$sub),
                 num_hab = 2, 
                 wood = as.numeric(hab$mean_w), 
                 ice_off = ice_off,
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



jags_model <- jags.model(textConnection(full_model), 
                         data = jags_data, 
                         n.chains = parms$n.chains, 
                         inits = inits)
jags_output.late <- coda.samples(jags_model, 
                                 variable.names = jags_parameters, 
                                 n.iter = parms$n.iter, 
                                 thin = parms$n.thin,
                                 inits = inits)
dat



chains <- as.mcmc.list(jags_output.late) ## set chains
chains_combined <- gtable_combine(chains) ## combine chains

dat = chains_combined %>%as.matrix() %>%  as.data.frame() #%>%

dat %>%
  ggplot(aes(x = Ntotal)) + 
  geom_histogram() + 
  scale_x_log10()

n.posts = 100

abund.dat =  chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("N", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = sum(value)) %>%
  
  left_join(legend, by = "group")

abund.dat %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot() +
  scale_y_log10()

legend


detection_probability = chains_combined %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  select(contains("p", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site","year","group")) %>%
  mutate(group = as.numeric(group),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, group) %>%
  summarize(tot_abund = mean(value)) %>%
  
  left_join(legend, by = "group") 


detection_probability %>%
  ggplot(aes(x = species, y = tot_abund)) + 
  geom_boxplot()



## Clean length estimates
mass.dat = dat %>%
  select(contains("avg_mass", ignore.case = TRUE))  %>%
  as.data.frame() %>%
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]),
               names_to = "metric", values_to = "mass.avg") %>%
  separate(metric, into = c("metric1","metric"), sep = "_") %>%
  mutate(group = parse_number(metric)) %>%
  select(-metric1, -metric) %>%
  left_join(legend)

# Figure for looking at the mass data
mass.dat %>%
  
  
  ggplot(aes(x = species, y = exp(mass.avg)))  +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_log10()+
  
  ylab("Avg. Mass (g)") +
  xlab("Species") +
  labs(col = "Species")

# Clean mu estimates
mu.dat = dat %>%
  as.data.frame() %>% 
  select(contains("mu", ignore.case = TRUE)) %>%
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]), names_to = "metric", values_to = "mu.value") %>% 
  separate(metric, into = c("metric", "group", "isotope")) %>% 
  unite("names", c(metric, isotope)) %>%
  pivot_wider(names_from = names, values_from = mu.value) %>%
  mutate(group = as.numeric(group)) %>%
  left_join(legend)
mu.dat %>% 
  ggplot(aes(x = species, y = mu_2)) + 
  geom_boxplot()

