### Early isotope values ----------------------------
data = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Group %in% "Fish") %>%
  select(Species, C, N) %>% 
  group_by(Species) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  arrange(Species) %>% 
  mutate(group = as.numeric(as.factor(Species)),
         community = 1) %>%
  select(community, group, C, N) %>%
  rename("iso_1" = "C",
         "iso_2" = "N") %>%
  
  ungroup() %>%
  na.omit()  %>%
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
save(arr, file = "IsotopeArray_early.RData")


backtrace.early = read.csv("Data/JML.Data.Master.csv") %>%
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
