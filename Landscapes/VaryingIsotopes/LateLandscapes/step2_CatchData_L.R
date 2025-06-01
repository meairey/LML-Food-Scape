library(tidyverse)
library(rjags)
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape/")
source("Landscapes/VaryingIsotopes/LateLandscapes/step1_LML_source_L.R")


## I'm a little worried that the multiple observation issue is going to catch up again here. But, for now, I'm simplifying this by just taking average catch and average effort for each site in each year. Not putting in each separately.

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
  filter(SPECIES %in% species) %>%
  group_by(SPECIES)



species.late = (step_a %>% 
  group_by(SPECIES) %>%
  summarize(total_catch = sum(mean_catch)) %>%
  filter(total_catch > 10))$SPECIES


### New model data setup --------------------------------------------


## Load in the data that says which sites to pair
rep_group = read.csv("Data/rep_groups.csv") %>%
  mutate(index = c(1:length(.[,1])))

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

shoreline_length = rep_group %>% arrange(site) %>%
  select(shoreline_length)


## Abundance Array -----------------------------------
## An intermediate step for creating the array
df = step_a %>% 
  left_join(rep_group, by = c("SITE" = "site")) %>%
  filter(SPECIES %in% species.late) %>%
  mutate(mean_catch = round(mean_catch, digits = 0)) %>%
  select(YEAR, SPECIES, index,  mean_catch) %>%
  na.omit()%>%
  
  group_by(YEAR,SPECIES, index) %>%
  mutate(rep_num = c(1:length(mean_catch))) %>%
  ungroup() 

# Get unique values for each dimension
years <- sort(unique(df$YEAR))
species <- unique(df$SPECIES)
rep_groups <- sort(unique(df$index))
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
    as.character(row$index),
    as.character(row$rep_num),
    as.character(row$YEAR),
    as.character(row$SPECIES)
  ] <- row$mean_catch
}

save(mean_catch_array, file = "Data/VaryingIsotopesData/LateData/mean_catch_array_late.RData")


## Effort frame from old model --------------------------------------
# Setup the catch data

species_cpue.data = step_a %>%
  select(-mean_effort) %>%
  mutate(mean_catch = round(mean_catch, digits = 0)) %>%
  pivot_wider(values_from = mean_catch, names_from = SITE) 


# Setup the effort data frame
effort = step_a %>% 
  ungroup() %>%
  select(-SPECIES, -mean_catch) %>%
  unique() %>%
  #mutate(mean_effort = normalize(mean_effort)) %>%
  mutate(mean_effort = mean_effort) %>%
  #mutate(mean_effort = scale(mean_effort, center = min(mean_effort), scale = max(mean_effort) - min(mean_effort))) %>%
  #mutate(mean_effort = round(mean_effort, digits = 0))%>%
  pivot_wider(values_from = mean_effort, names_from = SITE) %>%
  column_to_rownames(var = "YEAR") %>%
  mutate_all(~replace(., is.na(.), epsilon)) %>% # Make all NA's a very small value
  t()
# Save object
save(effort, file = "Data/VaryingIsotopesData/LateData/late_effortdata.RData")




