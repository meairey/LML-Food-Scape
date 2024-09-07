library(tidyverse)

source("LML_source.R")
## Load in weight frame from length_weight.R script


weight.frame = read.csv("Data/weight_frame.csv", header = T)

## using estimated weights
observed_lengths.early = weight.frame %>% 
  arrange(SPECIES) %>%
  select(SPECIES, YEAR,weight_e) %>% 
  filter(SPECIES %in% species.early) %>%
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

save(file = "Data/EarlyData/observed_lengths_early.RData", observed_lengths.early)

#how many length observations for each species
# Count the number of non-NA values in each row
n_obs_mass.early = (observed_lengths.early %>%
  rowwise() %>%
  mutate(n = sum(!is.na(c_across(everything())))) %>%
  ungroup())$n

save(file = "Data/EarlyData/n_obs_mass_early.RData", n_obs_mass.early)
