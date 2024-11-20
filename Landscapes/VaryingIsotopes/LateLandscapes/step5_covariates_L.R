source("Landscapes/VaryingIsotopes/LateLandscapes/step1_LML_source_L.R")
## Setting up covariates for JAGS Early
n_years = year_max - year_min
n_species = length(species.late)
## Filtering out the year range
temp = temp.full %>%
  filter(YEAR %in% (year_min:year_max), YEAR != 2002)
covariate_T = temp$mean_y %>% scale(, center = 0) %>% as.numeric()

save(file = "Data/LateData/covariate_temp.RData", covariate_T)

## Habitat data covariates 
covariate_hab = (shoreline_length %>% filter(Water == "LML"))$Habitat %>%
  as.factor() %>%
  as.numeric()
save(file = "Data/LateData/habs_cov.RData", covariate_hab)

n_sites = length(covariate_hab)

z_covariate <- sample(c(0, 1), n_species * n_years * n_sites, replace = TRUE)
save(file = "Data/LateData/z_covariate.RData", z_covariate)


## Ice Off covariate -- new model 
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
save(ice_off, file = "Data/VaryingIsotopesData/LateData/ice_off_late.RData")



## Habitat variables
rep_group = read.csv("Data/rep_groups.csv") ## Round values for modifying p
hab = rep_group %>%
  select(rep_group, hab) %>%
  na.omit() %>%
  mutate(sub = substr(hab, 1, 1),
         wood = substr(hab, 2, 2)) %>%
  mutate(wood = as.numeric(as.factor(wood)),
         sub  = as.numeric(as.factor(sub))) %>%
  select(rep_group, sub, wood) %>%
  group_by(rep_group, sub) %>%
  summarize(mean_w = mean(wood)) %>%
  unique()
save(hab, file = "Data/VaryingIsotopesData/LateData/hab_cov.RData")



## Habitat length 
## Calculate what proportion of a shoreline is sampled during each replicate
shoreline_length = rep_group %>%
  select(shoreline_length, rep_group) %>%
  na.omit() %>% group_by(rep_group) %>%
  mutate(count = c(1:length(shoreline_length))) %>% ## index for pivoting wider
  pivot_wider(names_from = count, values_from = shoreline_length) %>% ## Pivot it so reps are columns
  ungroup() %>%
  select(-rep_group) %>% ## remove columns not needed in RJAGS
  mutate(total_length = rowSums(., na.rm = T)) %>% ## Total length of shoreline for calculating proportion
  mutate(`1` = `1` / total_length,  ## Calculate the proportions for each replicate for each site
         `2` = `2` / total_length)
save(shoreline_length, file = "Data/VaryingIsotopesData/LateData/shoreline_length.RData")
