### Late Isotope values ---------------------



data.late = read.csv("Data/SI_MEASUREMENT.csv") %>%
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  filter(WATER == "LML",
         YEAR > 2018, 
         GROUP == "FISH",
         SAMPLE_TYPE == "TISSUE") %>%
  select(TAXON, D13C, D15N) %>% 
  unique() %>% 
  arrange(TAXON) %>%
  rename("iso_1" = "D13C",
         "iso_2" = "D15N",
         "group" = TAXON) %>%
  group_by(group) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit()  %>%
  as.data.frame() %>%
  mutate(community = 1, 
         iso_1 = scale(as.numeric(iso_1)),
         iso_2 = scale(as.numeric(iso_2))) %>%
  na.omit() %>%
  select(community, group, iso_1, iso_2) %>%
  filter(group %in% c("CC", "CS","LT","MM","PS","RS","SMB","SS","WS"))





### Filtering in the species (modifying from early data)

### Early isotope values ----------------------------
data.early = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Group %in% "Fish") %>%
  select(Species, C, N) %>% 
  group_by(Species) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  arrange(Species) %>% 
  mutate(  community = 1,
           group = Species) %>%
  select(community, group, C, N) %>%
  rename("iso_1" = "C",
         "iso_2" = "N") %>%
  
  ungroup() %>%
  na.omit()  %>%
  mutate(iso_1 = scale(as.numeric(iso_1)),
         iso_2 = scale(as.numeric(iso_2))) %>%
  na.omit() %>% 
  filter(group %in% c("BB"))




##  Combining data ------------------

data = rbind(data.early, data.late)


data %>% ggplot(aes(x = iso_1, iso_2, col = group)) + 
  geom_point() + stat_ellipse()

df = data %>% 
  arrange(group) %>%
  mutate(group = as.numeric(as.factor(group))) %>%
  group_by(group) %>% 
  mutate(num = c(1:length(group))) %>%
  ungroup() %>%
  complete(group, num) %>%
  rename("species" = "group")


n_species = length(unique(df$species))
n.iso = 2
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


save(file = "Data/VaryingIsotopesData/LateData/IsotopeArray_late.RData", arr)


## backtrace late

backtrace.late = read.csv("Data/SI_MEASUREMENT.csv") %>%
  
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  filter(WATER == "LML",
         YEAR > 2018, 
         GROUP == "FISH",
         SAMPLE_TYPE == "TISSUE") %>%
  select(TAXON, D13C, D15N) %>% 
  unique() %>% 
  arrange(TAXON) %>%
  rename("iso_1" = "D13C",
         "iso_2" = "D15N",
         "group" = TAXON) %>%
  group_by(group) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit()  %>%
  as.data.frame() %>%
  mutate(community = 2) %>%
  group_by(community) %>%
  summarize(mean_C = mean((iso_1)), sd_C = sd(iso_1),
            mean_N = mean((iso_2)), sd_N = sd(iso_2)) 
