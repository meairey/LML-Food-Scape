library(tidyverse)

source("Landscapes/VaryingIsotopes/PreLandscapes/step1_LML_source_pre.R")
## Load in weight frame from length_weight.R script


weight.frame = read.csv("Data/weight_frame.csv", header = T)

## using estimated weights
observed_lengths.pre= weight.frame %>% 
  arrange(SPECIES) %>%
  select(SPECIES, YEAR,weight_e) %>% 
  filter(SPECIES %in% species.pre) %>%
  filter(YEAR <= year_max & YEAR >= year_min, 
         is.na(weight_e) == F) %>%
  
  arrange(SPECIES, YEAR) %>% 
  group_by(SPECIES, YEAR) %>%
  slice_head(n = 10) %>%
  mutate(YEAR =1) %>%
  group_by(SPECIES) %>%
  mutate(individuals = c(1:length(SPECIES)))%>%
  pivot_wider(names_from = individuals, values_from = weight_e) %>% 
  as.data.frame() %>%
  select(-YEAR) %>%
  column_to_rownames(var = "SPECIES")

t(observed_lengths.pre) %>% 
  as.data.frame() %>%
  mutate(ind = 1:100) %>%
  select(ind, everything()) %>%
  pivot_longer(2:length(.[1,])) %>%
  group_by(name) %>%
  summarize(mean = mean(value))

save(file = "Data/VaryingIsotopesData/PreData/observed_lengths_pre.RData", observed_lengths.pre)

#how many length observations for each species
# Count the number of non-NA values in each row
n_obs_mass.pre = (observed_lengths.pre %>%
  rowwise() %>%
  mutate(n = sum(!is.na(c_across(everything())))) %>%
  ungroup())$n

save(file = "Data/VaryingIsotopesData/PreData/n_obs_mass_pre.RData", n_obs_mass.pre)





