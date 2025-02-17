library(rjags)
library(tidyverse)

setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape/")
source("Landscapes/VaryingIsotopes/PreLandscapes/step1_LML_source_pre.R")

## I'm a little worried that the multiple observation issue is going to catch up again here. But, for now, I'm simplifying this by just taking average catch and average effort for each site in each year. Not putting in each separately.






### New Model ------------------------------------------


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


sample %>%
  filter(YEAR == 2000, MONTH < 7, WATER == "LML", GEAR_CODE == "NAF") %>%
  ggplot(aes(y = SITE_N, x = DAY_N)) +
  geom_point() +
  facet_wrap(~SITE_N, scales = "free_y")

step_a %>% 
  group_by(SPECIES) %>%
  summarize(total_catch = sum(mean_catch))



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
  dim = c(length(sites), 1,length(rep_nums),  length(species)),
  dimnames = list(sites = sites, "replicates" = 1 , "year" = rep_nums,  SPECIES = species)
)

# Fill the array
for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  mean_catch_array[
    as.character(row$SITE),
    1,
    as.character(row$rep_num),
    
    as.character(row$SPECIES)
  ] <- row$mean_catch
}

dim(mean_catch_array)

mean_catch_array %>% dim()

mean_catch_array = mean_catch_array[, , 1:5,,drop = FALSE ]
save(mean_catch_array, file = "Data/VaryingIsotopesData/PreData/mean_catch_array_pre.RData")


df %>%
  group_by(SITE) %>%
  summarize(max_rep = max(rep_num)) %>%
  ungroup() %>%
  summarize(summary = median(max_rep))









## Old Model --------------------------------------------------

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
save(effort, file = "Data/PreData/pre_effortdata.RData")

### Species data matrix 
# Get the list of unique species
species_list <- unique(step_a$SPECIES)
# Initialize an empty list to store matrices
species_matrices <- list()
# Loop over each species
for (species in species_list) {
  # Filter data for the current species
  species_data <- species_cpue.data %>%
    filter(SPECIES == species) %>%
    ungroup() %>%
    select(-SPECIES)  # Remove the species column
  
  # Reshape the data to have years as rows and sites as columns
  
  ## Hmm the bef.combined.txt seems to want years as columns and sites as rows.
  species_matrix <- as.matrix(species_data %>%
                                pivot_longer(cols = starts_with("BEF"), names_to = "Site", values_to = "Count") %>%
                                mutate(Count = replace_na(Count, 0)) %>% ## Added in to remove NAs for dnbinom
                                pivot_wider(names_from = Site, values_from = Count) %>%
                                column_to_rownames(var = "YEAR"))  %>% t() # For some reason this is here but it makes the matrix opposite of what I think it should be so I've commented it out?
  
  # Add the matrix to the list
  species_matrices[[species]] <- species_matrix
}


# Get unique dimensions
# I changed this around because of comment above that i think ive flipped the dimensions for the array Y
years <- colnames(species_matrices[[1]])
sites <- rownames(species_matrices[[1]])
species <- names(species_matrices)

# Initialize the array

# redoing the array to have years as columns and rows as sites
array_data <- array(0, dim = c(length(sites), length(years), length(species)),
                    dimnames = list(sites, years, species))
# Fill the array
for (i in seq_along(species)) {
  array_data[,,i] <- species_matrices[[i]]
}

save(file = "Data/PreData/Pre_arraydata.RData", array_data)






