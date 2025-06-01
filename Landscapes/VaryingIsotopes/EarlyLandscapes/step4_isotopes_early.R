## Checking dataframe
## Access Database

data = read.csv("../1.clean_isotope/SI_MEASUREMENT.csv")
sample.clean = read.csv("../1.clean_isotope/SI_SAMPLE.csv")
da.acc = left_join(data, sample.clean, by = "ISO_YSAMP_N") %>%
  select(-corrected,-Corrected)
## Missing samples from JML File

mi.sa = read.csv("Data/IsotopeComparisons/Clean/missing_sample.csv") %>% unique()
mi.me = read.csv("Data/IsotopeComparisons/Clean/missing_measurement.csv")
miss.dat = left_join(mi.me, mi.sa, by = "ISO_YSAMP_N") %>% filter(YEAR < 2005)

data = rbind(da.acc, miss.dat)

data %>%
  filter(GROUP %in% "FISH", WATER == "LML", YEAR < 2005) %>%
  select(TAXON, D13C, D15N) %>% 
  group_by(TAXON) %>%
  filter(n() >= 3) %>%
  select(TAXON) %>% unique()

## Load in access LML measurement file for other lengths

`%nin%` = Negate(`%in%`)

### Early isotope values ----------------------------
data = data %>%
  filter(GROUP %in% "FISH", WATER == "LML", YEAR < 2005) %>%
  select(TAXON, D13C, D15N) %>% 
  group_by(TAXON) %>%
  na.omit()  %>%
  filter(TAXON %in% species.early) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  arrange(TAXON) %>% 
  mutate(group = as.numeric(as.factor(TAXON)),
         community = 1) %>%
  select(community, group, D13C, D15N) %>%
  rename("iso_1" = "D13C",
         "iso_2" = "D15N") %>%
  
  ungroup() %>%
  
  mutate(iso_1 = scale(as.numeric(iso_1)),
         iso_2 = scale(as.numeric(iso_2))) %>%
  na.omit()



df= data %>% 
  group_by(group) %>% 
  mutate(num = c(1:length(group))) %>%
  ungroup() %>%
  #group_by() %>%
  complete(group, num) %>%
  rename("species" = "group")


n_species = df$species %>% unique() %>% length()




# Create an empty 3D array
arr <- array(NA, dim = c(n_species, max(df$num), n.iso))

# Iterate over species and fill the array
for (i in 1:n_species) {
  # Subset data for the current species
  subset_df <- subset(df, species == unique(df$species)[i])
  
  # Fill the array with delN and delC values
  arr[i, , ] <- matrix(c(subset_df$iso_1, subset_df$iso_2), ## I'm wondering if the order here influences the graph down the line because it looks like things are flipped around...? 
                       nrow = nrow(subset_df), byrow = FALSE)
}

# Save the array as an R object (RData file)
save(arr, file = "Data/VaryingIsotopesData/EarlyData/IsotopeArray_early.RData")


backtrace.early = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Group %in% "Fish") %>%
  select(Species, C, N) %>% 
  arrange(Species) %>% 
  mutate(spp = as.numeric(as.factor(Species)),
         community = 2, 
         C = as.numeric(C),
         N = as.numeric(N)) %>%
  select(community, spp, C, N) %>%
  rename("iso_1" = "C",
         "iso_2" = "N") %>%
  group_by(spp) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(community) %>%
  summarize(mean_C = mean((iso_1)), sd_C = sd(iso_1),
            mean_N = mean((iso_2)), sd_N = sd(iso_2)) 
save(backtrace.early, file = "Data/VaryingIsotopesData/EarlyData/backtrace_early.RData")


## Pre data because they use the same isotope data

backtrace.pre = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Group %in% "Fish") %>%
  select(Species, C, N) %>% 
  arrange(Species) %>% 
  mutate(spp = as.numeric(as.factor(Species)),
         community = 1, 
         C = as.numeric(C),
         N = as.numeric(N)) %>%
  select(community, spp, C, N) %>%
  rename("iso_1" = "C",
         "iso_2" = "N") %>%
  group_by(spp) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(community) %>%
  summarize(mean_C = mean((iso_1)), sd_C = sd(iso_1),
            mean_N = mean((iso_2)), sd_N = sd(iso_2)) 
save(backtrace.early, file = "Data/VaryingIsotopesData/PreData/backtrace_pre.RData")

