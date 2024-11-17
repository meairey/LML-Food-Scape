

## Simple N-mixture model

sites = 150 ## Spatial replicate
lamba = 2.5 # The expected number of individuals at any site
set.seed(123) ## Random seed


N = rpois(sites, lamba) ## Actual counts observed during survey
plot(table(N)) ## Plotting so you can see that the max number of individuals seen an any site is 7, but most often we see 2 animals per site


## Observation model 

## Just use a binomial

temp = 3 ## Temporal replicates

p = .6 ## Probability of detecting an individual


C = matrix(NA, sites, temp)
for(i in 1:temp) {
  C[, i] <- rbinom(sites, N, p)
}


Cmax = apply(C, 1, max)


## Loading in the data to run the analysis

jags_data = list(sites = sites,
                 temp = temp,
                 C = C)

str(jags_data)

jags_parameters= c("lambda", "p", "Ntotal")

library(rjags)
# options for running jags
parms <- list(); parms$n.iter=2*10^4; parms$n.burnin=1*10^3; 
parms$n.thin=10; parms$n.chains=2  

set.seed(123)
model = read_file("Data/Models/n_mixture_simple_example.txt")

inits = function() list(N = Cmax)
str(inits())

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

dat = chains_combined %>%as.matrix() %>%  as.data.frame()

(dat$Ntotal / 150) %>% mean(.)

### Trying the new structure with multiple different replicates per site but different averages per site... three dimensional C/observation matrix



## Fake data

lambda = c(1, 5, 7, 10)
sites = length(lambda)
replicates = 3
N = matrix(NA, nrow = sites, ncol = replicates)
for(i in 1:sites){

  N[i,] = rpois(replicates, lambda[i])
  }
}

apply(N,1, mean)
N


C = array(NA, dim = c(sites, replicates, temp))


for(i in 1:temp) {
  for(h in 1:sites){
    
  C[h, ,i] <- rbinom(replicates, N[h,], p)
  
  }
}


Cmax = apply(C, 1, max)


model = read_file("Data/Models/n_mixture_siterepyear.txt")

## Loading in the data to run the analysis

jags_data = list(sites = sites,
                 temp = temp,
                 replicates = replicates, 
                 C = C)

str(jags_data)

jags_parameters= c("lambda", "p", "Ntotal","N")


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

dat = chains_combined %>%as.matrix() %>%  as.data.frame()

dat %>%
  ggplot(aes(x = Ntotal)) + 
  geom_bar()







### Using code from step2_CatchData_L.R





rep_group = read.csv("Data/rep_groups.csv") %>%
 
  mutate(mean_length = mean(shoreline_length)) %>%
  mutate(shoreline_length = round(shoreline_length, digits = 0)) %>%
  mutate(std_shoreline = round(shoreline_length/mean_length, digits = 1))

step_b = step_a %>% 
  left_join(rep_group, by = c("SITE" = "site")) %>%
  filter(SPECIES == "CS") %>%
  mutate(mean_catch = round(mean_catch, digits = 0))



shoreline_length = rep_group %>% select(shoreline_length, rep_group_simple) %>%
  na.omit() %>% group_by(rep_group_simple) %>%
  mutate(count = c(1:length(shoreline_length))) %>%
  pivot_wider(names_from = count, values_from = shoreline_length) %>%
  ungroup() %>%
  select(-rep_group_simple) %>%
  mutate(total_length = rowSums(.)) %>%
  mutate(`1` = `1` / total_length, 
         `2` = `2` / total_length)


mean_length = mean(rep_group$shoreline_length) %>% round(., digits = 0)

y1 = step_b %>% 
  filter(YEAR == 2001) %>% 
  select(mean_catch, rep_group_simple) %>%
  group_by(rep_group_simple) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  na.omit() %>%

  pivot_wider(names_from = rep_num, values_from = mean_catch) %>% 
  
  ungroup() %>%
  select(-rep_group_simple) %>%
  as.matrix() 

y2 =  step_b %>% 
  filter(YEAR == 2003) %>% 
  select(mean_catch, rep_group_simple) %>%
  group_by(rep_group_simple) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  na.omit() %>%
  pivot_wider(names_from = rep_num, values_from = mean_catch) %>% 
  
  ungroup() %>%
  select(-rep_group_simple) %>%
  as.matrix() 


y3 =  step_b %>% 
  filter(YEAR == 2004) %>% 
  select(mean_catch, rep_group_simple) %>%
  group_by(rep_group_simple) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  na.omit() %>%
  pivot_wider(names_from = rep_num, values_from = mean_catch) %>% 
  
  ungroup() %>%
  select(-rep_group_simple) %>%
  as.matrix() 



C = array(c(y1, y2, y3), dim = c(13,2,3))
 
C[which(is.na(C))] <- NA 


sites =  dim(y1)[1]
year = 3
replicates = 2

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


### Trying out new model
jags_data = list(sites = sites,
                 year = year,
                 replicates = replicates, 
                 C = C, 
                 site_proportion = shoreline_length[,1:2] ,
                 hab = as.numeric(hab$sub),
                 num_hab = 2, 
                 wood = as.numeric(hab$mean_w), 
                 ice_off = ice_off)


Cmax = apply(C, 1, max,na.rm=T) ## Set up innit to help model converge. This sets initial values at the max abundance counted at any site



jags_parameters= c("lambda", "p", "N","Ntotal")

model = read_file("Data/Models/expected_counts_replicate_length.txt")

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

mean(dat$Ntotal)
## CC = ] 2792.199 (late) -> 1752.752 (early)
## CS [1] 5374.574 -->  6355.459
## SMB [1] 4383.767 --> [1] 3755.756

dat %>%
  ggplot(aes(x = p)) + 
  geom_histogram()
dat$Ntotal %>% mean()


step_b %>%
  group_by(YEAR) %>%
  summarize(mean_catch = sum(mean_catch))

step_b %>%
  group_by(YEAR, rep_group_simple) %>%
  summarize(mean_catch = sum(mean_catch))

## Temperatures? 
sample = read.csv("../AFRP/Data/FISH_SAMPLE_2024.csv") %>%
  filter(YEAR %in% c(2019, 2021, 2022), 
         MONTH %in% c(5, 6),
         WATER == "LML", GEAR_CODE == "NAF")
temp = sample %>% 
  left_join(rep_group, by = c("SITE_N"="site")) %>% 
  select( YEAR, TEMP_SAMP, rep_group_simple) %>%
  na.omit() %>%
  group_by(rep_group_simple, YEAR) %>%
  mutate(count= 1:length(TEMP_SAMP)) %>%
  
  pivot_wider(names_from = count, values_from = TEMP_SAMP) 


k =temp %>% 
  ungroup() %>% 
  select(YEAR, `1`, `2`) 


# Split the data into a list of matrices by YEAR
matrices_list <- k %>%
  group_by(YEAR) %>%
  group_split() %>%
  lapply(function(sub_df) as.matrix(sub_df[, c("1", "2")]))

# Combine the list into an array (3 matrices, 13 rows, 2 columns)
temp_array <- array(unlist(matrices_list), dim = c(13, 2, 3))
temp_array[9,2,3] = 16
temp_array[10,2,3] = 15



air_temp = read.csv("Data/1_AirTempMod_V3_MeanMonthlyTemp.csv")
air_temp %>% filter(YEAR %in% 2019)

sample = read.csv("Data/FISH_SAMPLE_2024.csv")

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

ice_off = as.numeric(ice_off$avg.day.ice[4:6])




step_a  %>%
  filter(YEAR == 2001) %>% select(SITE) %>%
  unique() %>% print(n = 1000)
                          
