library(raster) ## This package messes with the select function in library dplyr

library(snow)
library(tidyverse)
library(viridis)
## Source and Functions
source("Data/Models/functions.R")

## script setup
n_species = 10
n.posts = 20
cord_min = -7; cord_max = 7
# Resolution
xy_length = 50

## Legend
species= c("BB", "CC","CS","LT","MM","PS","RS","SMB","SS","WS")
species.pre = c("BB", "CC","CS","PS","SMB","SS","WS")
legend = data.frame(group = as.character(c(1:10)), 
                    species  = species, 
                    species.pre = c(species.pre, rep("NA", 3)))
#---------------------------------------------------------------------------



load("Data/UnvariedIsotopeData/posterior_early.unvaried.csv")
load("Data/UnvariedIsotopeData/posterior_late.RData")
load("Data/UnvariedIsotopeData/posterior_pre.RData")

## Comparison of posterior parameter distributions
posterior.pre = posterior.pre.unvaried %>%
  rename("group" = "species") %>% 
  left_join(legend) %>%
  mutate(community = 1) %>%
  mutate(species = species.pre) %>%
  ungroup() %>%
  select(-species.pre, -group)


posterior.early = posterior.early.unvaried %>%
  rename("group" = "species") %>% 
  left_join(legend) %>%
  mutate(community = 2) %>%
  ungroup() %>%
  select(-species.pre,-group) 


posterior.late = posterior.late.un %>%
  rename("group" = "species") %>% 
  left_join(legend) %>%
  mutate(community = 3) %>%
  ungroup() %>%
  select(-species.pre, -group)

posterior.late %>% group_by(species) %>%
  summarize(mean_1 = mean(mu_1), mean_2 = mean(mu_2))
  

posterior = rbind(posterior.pre, posterior.early, posterior.late)

posterior %>% 
 

  ggplot(aes(x = as.factor(community), y = tot_abund)) + 
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~species, scales = "free_y") +
  xlab("Period") +
  ylab("Total Abundance") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 

##

posterior %>% 
  
  ggplot(aes(x = as.factor(community), y = mass.avg)) + 
  geom_boxplot() +
  #scale_y_log10() +
  facet_wrap(~species, scales = "free_y") +
  xlab("Period") +
  ylab("Average Mass (g)") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 


# Ellipse generation ----------

## Coord limits


workers = 5


### Mean posterior heatmap for visualization ------

## New type of heat map with averages from posterior




posterior = rbind(posterior.pre, posterior.early, posterior.late) %>% 
  group_by(community, species) %>%
  summarize(Sigma_1_1 = mean(Sigma_1_1), 
            Sigma_2_1 = mean(Sigma_2_1), 
            Sigma_1_2 = mean(Sigma_1_2),
            Sigma_2_2 = mean(Sigma_2_2),
            mass.avg = mean(mass.avg), 
            tot_abund = mean(tot_abund), 
            mu_C = mean(mu_1),
            mu_N = mean(mu_2))  %>%
  select(community, species, Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2, mass.avg, tot_abund, mu_C, mu_N) %>%
  ungroup() %>%
  arrange(species) %>%
  mutate(group = as.numeric(as.factor(species))) %>%
  complete(group, community) %>%
  select(-species) %>%
  rename("species" = "group") %>%
  mutate(post = 1)

posterior_isotopes.summary = posterior %>% filter(community == 3) %>%
  select(species, post, mu_C, mu_N, Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2)
  
legend


posterior %>% print(n = 100)



## ellip dataframe for the mean landscape (visualization only)
ellip.mean = ellip.data(xy_length, 1, length(unique(posterior$species)), length(unique(posterior$community)))
# Cluster to create heights
clo = makeCluster(workers) # Create workers
clusterEvalQ(clo, c(library(mvtnorm), library(SIBER), library(dplyr))) # Set up packages on the cluster
clusterExport(clo, list = c("posterior", "ellipfunc")) # #Export dependants
ellip.mean$string = parRapply(clo,ellip.mean, # run the ellipfun across the coords
                              function(x) ellipfunc(xax = x[3],
                                                    yax = x[4],
                                                    post_n = x[2],
                                                    spp = x[1],
                                                    community_n = x[5]))
stopCluster(clo) ## stop workers/cluster

## 
ellip.mean %>% filter(spp == 10, community == 3)

# Filter the data frame string
ellip.mean.filtered = ellip.mean %>% 
  na.omit() %>%
  filter_ellip.data(., 10, 3) 



ellip.mean %>% group_by(community, spp) %>%
  summarize(total = sum(string)) %>%
  arrange(community, -total) %>%
  print(n = 100)

ellip.mean.filtered %>% filter(spp == 10, community == 3)

ellip.mean.filtered.un = ellip.mean.filtered
save(file = "Data/UnvariedIsotopeData/MeanEllipses.RData", ellip.mean.filtered.un)
load(file = "Data/UnvariedIsotopeData/MeanEllipses.RData")
# Facet wrap heat maps to see individual ellipses



ellip.mean.filtered %>% 

  mutate(spp = as.character(spp)) %>%
  left_join(legend,by = c("spp" = "group")) %>% 
  ggplot(aes(x= xax, y = yax, col =as.factor(community))) +
  geom_jitter(alpha = .15) + 
  facet_wrap(~species) 




## Filtered heatmap with landscape for all species
ellip.mean.filtered %>%
  #left_join(back.trace) %>%
  #mutate(xax = xax*sd_C + mean_C, yax = yax*sd_N + mean_N) %>%
  group_by(community, xax, yax) %>%
  summarize(vol = sum(string)) %>%
  ungroup() %>%
  ggplot(aes(x = xax, y = yax, col = (vol))) +
  geom_point(size = 5) +

  scale_color_viridis() +
  theme_minimal() + 
  labs(fill = "Z height") +
  xlab("d13C") +
  ylab("d15N") + 
  facet_wrap(~community)



posterior %>%
  arrange(species)

# Full landscape --------------------------------------

## Ellip coordinates ------------
## Redefine posterior because we had to overwrite it for the past chunk. The only other solution I can think of here is creating a second function that does the average visualization
posterior = rbind(posterior.pre, posterior.early, posterior.late) %>%
  
  rename("mu_C" = "mu_1", ## need to make sure have conforming names with function
         "mu_N" = "mu_2")%>%
  select(community, post,  species, Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2, mass.avg, tot_abund, mu_C, mu_N) %>%
  ungroup() %>%
  arrange(species) %>%
  mutate(group = as.numeric(as.factor(species))) %>%
  complete(group, community, post) %>%
  select(-species) %>%
  rename("species" = "group")  
         





ellip = ellip.data(xy_length, n.posts, 10,3)



# Cluster to create heights
clo = makeCluster(workers) # set up cluster
clusterEvalQ(clo, c(library(mvtnorm), library(SIBER), library(dplyr))) # export libraries to the cluster
clusterExport(clo, list = c("posterior", "ellipfunc")) # export the functions and data frames to the cluster
ellip$string = parRapply(clo,ellip, ## Iterate through ellip dataframe with the ellipfunc
                         function(x) ellipfunc(xax = x[3],
                                               yax = x[4],
                                               post_n = x[2],
                                               spp = x[1],
                                               community_n = x[5]))
stopCluster(clo) # Stop workers


## Filtering the ellip data

filtered_data= filter_ellip.data(ellip, n_species,3)

## Save .RData files
ellip.un = ellip
filtered_data.un = filtered_data
save(file = "Data/UnvariedIsotopeData/unfiltered_ellipses.RData", ellip.un)
save(file = "Data/UnvariedIsotopeData/filtered_ellipses.RData", filtered_data.un)

## Load full landscape ----------------
load("Data/UnvariedIsotopeData/filtered_ellipses.RData")


ellip.filtered = filtered_data %>%
  mutate(spp = as.character(spp)) %>%
  left_join(legend, by = c("spp" = "group")) %>%
  rename("code" = "species")




# Volume metrics and boxplot ---------------------------------
# Total volume of community 
total_vols = ellip.filtered %>%
  group_by(code, post, community) %>%
  summarise(total_vol = sum(string)) 


## Total volume of each species individually
total_vols %>%
  as.data.frame() %>%
  complete(code, post, community) %>%
  replace_na(list(total_vol = 0)) %>%
  ggplot(aes(x = code, y = total_vol +1, col = community %>% as.factor())) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Total Volume") +
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_log10() +
  theme_minimal(base_size = 14) +
  labs(col = "Period") +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late"))


# Total volume of the entire community


total_vols %>% 
  ungroup() %>%
  group_by(post, community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = community %>% as.factor(), y = total_vol, col = as.factor(community))) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  ylab("Total Volume") +
  xlab("Community") +
  theme(legend.position = "none") 
  





# Peak N and C Densities ------------------------------
## This has been back transformed



axis_slices = ellip.filtered %>% 
  group_by(post, xax, yax, community) %>%
  summarize(total_volume = sum(string)) %>%
  #left_join(back.trace) %>% 
 # mutate(xax = xax * sd_C + mean_C, 
  #       yax = yax * sd_N + mean_N) %>%
  group_by(post,community) %>%
  slice_max(total_volume) %>%
  select(post, xax, yax, community) %>% 
  pivot_longer(c("xax", "yax"), names_to = "Type", values_to = "Values")




axis_slices %>% 
  mutate(Type = str_replace(Type, "xax", "D13C"),
         Type = str_replace(Type, "yax", "D15N")) %>%
  ggplot(aes(x = as.factor(community), y = Values, col = as.factor(community))) + 
  geom_boxplot()+ 
  ylab("Axis value") + 
  facet_wrap(~ Type, scales = "free")  +
  labs(col = "Period") +
  scale_x_discrete("Community", labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") 




# Distance from center of ellipse to peak density ----------------------------
#### Needs to be modified in the distances function




# set up matrix for filtering within the for loop
z.mod = ellip.filtered %>% 
  group_by(post, xax, yax, community) %>%
  summarize(total_volume = sum(string)) %>%
 # left_join(back.trace) %>% ## Back transforming the axes
 # mutate(xax = xax * sd_C + mean_C, 
    #     yax = yax * sd_N + mean_N) %>%
  ungroup() %>%
  select(total_volume, xax, yax, post, community) 

## Trouble shooting


## for loop
dist_list = list()
dist_list1 = list()
for(i in 1:length(unique(z.mod$community))){
  for(h in 1:length(unique(z.mod$post))){
    
    
    
    z = z.mod %>%
      filter(community == i, post == h) %>%
      pivot_wider(values_from = total_volume, names_from = xax) %>%
      column_to_rownames(var = "yax") %>%
      mutate_all(~ ifelse(is.na(.), 0, .))
    
    # Apply the function to your Z matrix
    local_maxima <- find_local_maxima(z) 
    # Get coordinates of local maxima
    maxima_coords <- which(local_maxima, arr.ind = TRUE) ## 
    max = data.frame(xax = (z.mod$xax %>% unique())[maxima_coords[,1]], 
                     yax = (z.mod$yax %>% unique())[maxima_coords[,2]], 
                     max = "max") 
    
    
    max.frame.joined = left_join(z.mod, max)
    
    
    
    
    dist_list1[[h]] =  dist(max %>% select(-max), method = "euclidean") %>% as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column(., var = "row_index") %>%
      pivot_longer(1:length(.[,1])+1, names_to = "col_index", values_to = "values") %>%
      summarize(mean = mean(values), sd = sd(values))
  }
  
  dist_list[[i]] = dist_list1
}

dist_matrix = data.frame(value = c(unlist(dist_list[[1]]), unlist(dist_list[[2]]), unlist(dist_list[[3]]))) %>% 
  mutate(metric = rep(c("mean", "sd"), n.posts*3),
         post = rep(rep(1:n.posts, each = 2), 3), 
         community= rep(c(1,2, 3), each = (n.posts*2 )))

## Distance between all points in the landscape

dist_matrix %>% 
  ggplot(aes(x = as.factor(community), y = value, col = as.factor(community))) + 
  geom_boxplot() + 
  facet_wrap(~metric) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Value") + xlab("Period")  +
  labs(col = "Period") +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late"))

## Distance between next highest peak? ----------------------------
# I dont think you want to transform these values back into isotopically meaningful ones... it would be misleading if the variation isn't standardized...

## Demonstration of what it is doing



## Next highest point peak ----------------

z.mod %>% 
  group_by(community, post) %>%
  arrange(post, -total_volume) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values") %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Distance to next highest peak") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Pre", "2" = "Early", "3" = "Late")) +
  facet_wrap(~metric, scales = "free")


dist.graph.dat = z.mod %>% 
  
  group_by(community, post) %>%
  arrange(post, -total_volume) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  filter(community == 1, post == 1) %>%
  ungroup() %>%
  rename("string" = "total_volume")


ellip.filtered %>% filter(post == 1, community == 1)  %>%
  #left_join(back.trace)  %>%
  # mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
  #       yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  group_by(community, post, xax, yax) %>%
  summarize(sum = sum(string)) %>%
  ungroup() %>%
  ggplot(aes(x = xax, y = yax)) + 
  geom_tile(aes(fill = sum), alpha = .8) + 
  geom_point(data = dist.graph.dat[1:100,], mapping = aes(x = xax, y = yax)) + 
  geom_path(data = dist.graph.dat[1:100,], mapping = aes(x = xax, y = yax, col =string )) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  ylab("d15N") + xlab("d13C")


## Next highest species peak ---------------

ellip.filtered %>% 
  group_by(post, spp, community) %>%
  slice_max(string) %>%
  #left_join(backtrace)  %>%
  # mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
  # yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values") %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Distance to next highest peak") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Pre", "2" = "Early", "3" = "Late")) +
  facet_wrap(~metric, scales = "free")


dist.graph.dat = ellip.filtered %>% 
  group_by(post, spp, community) %>%
  slice_max(string) %>%
  #left_join(back.trace)  %>%
  #mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
  #       yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  filter(community == 2, post == 3) %>%
  ungroup() %>%
  mutate(index = c(1:10))

ellip.filtered %>% filter(post == 3, community == 2)  %>%
  #left_join(back.trace)  %>%
  # mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
  #       yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  group_by(community, post, xax, yax) %>%
  summarize(sum = sum(string)) %>%
  ungroup() %>%
  ggplot(aes(x = xax, y = yax)) + 
  geom_tile(aes(fill = sum), alpha = .8) + 
  geom_point(data = dist.graph.dat[1:100,], mapping = aes(x = xax, y = yax)) + 
  geom_path(data = dist.graph.dat[1:100,], mapping = aes(x = xax, y = yax, col =-index ), size = 1) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  ylab("d15N") + xlab("d13C") +
  scale_fill_gradient(high = "white", low = "black")


## Standard deviation of community -------------------------

z.mod %>% 
  group_by(post,community) %>%
  summarize(sd = sd(total_volume)) %>%
  ggplot(aes(x = as.factor(community), y = sd, col = as.factor(community))) +
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  xlab("Period") + 
  ylab("Standard Deviation") + 
  scale_x_discrete(labels = c("1" = "Early", "2" = "Late"))


## Standard deviation of individual species 

ellip.filtered %>%
  group_by(post, community, code) %>%
  summarise(mean_sd = sd(string)) %>%
  ggplot(aes(x = code, y = log(mean_sd))) + 
  geom_boxplot() + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust=1)) + 
  ylab("Ellipse Standard Deviation") + 
  xlab("Species") + 
  facet_wrap(~community) 



## Rugosity -----------------------------------------





cell_size = (abs(cord_min) + abs(cord_max)) / xy_length

#backtrace = back.trace %>% 
 # mutate(cell_size = (abs(cord_min) + abs(cord_max)) / 100) %>%
  #mutate(cell_size.C = cell_size*sd_C, 
   #      cell.size.N = cell_size * sd_N)






# Calculate the planar area for each cell


rug = list()
rug1 = list()

for(i in 1:length(unique(z.mod$community))){
  # Define cell dimensions
  #cell_height <- backtrace$cell.size.N[i]
  #cell_width = backtrace$cell_size.C[i] # Example height
  cell_height = cell_size
  cell_width = cell_size
  planar_area <- cell_height * cell_width
  for(h in 1:length(unique(z.mod$post))){
    
    elevation_matrix = z.mod %>%
      filter(community == i, post == h) %>%
      pivot_wider(values_from = total_volume, names_from = xax) %>%
      column_to_rownames(var = "yax") %>%
      mutate_all(~ ifelse(is.na(.), 0, .))
    
    # Calculate surface area
    surface_area <- calculate_surface_area(elevation_matrix)
    
    # Compute planar area
    # Note: Multiply by the number of cells to get total planar area
    total_cells <- ncell(elevation_matrix)
    total_planar_area <- total_cells * planar_area
    
    # Calculate rugosity
    rugosity <- surface_area / total_planar_area
    
    rug1[[h]] = rugosity
    
  }
  rug[[i]] = rug1
}




rugosity = data.frame(values = unlist(rug),
                      post = rep(c(1:n.posts), rep = length(unique(z.mod$community))),
                      community = rep(unique(z.mod$community), each = n.posts))

rugosity %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Rugosity") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Pre", "2" = "Early", "3" = "Late")) 

