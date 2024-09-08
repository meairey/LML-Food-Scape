## Source and Functions
source("Data/Models/functions.R")
load("Data/posterior_early.csv")

n_species = 10

# Ellipse generation ----------

## Coord limits

cord_min = -7; cord_max = 7
workers = 5
# Resolution
xy_length = 100


### Mean posterior heatmap for visualization ------

## New type of heat map with averages from posterior
#posterior = posterior.early
#posterior.save = posterior
posterior = posterior.early %>% 
  group_by(species) %>%
  summarize(Sigma_1_1 = mean(Sigma_1_1), 
            Sigma_2_1 = mean(Sigma_2_1), 
            Sigma_1_2 = mean(Sigma_1_2),
            Sigma_2_2 = mean(Sigma_2_2),
            length.avg = mean(mass.avg), 
            tot_abund = mean(tot_abund), 
            mu_C = mean(mu_1),
            mu_N = mean(mu_2)) %>%
  mutate(post = 1) %>%
  select(post, species, Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2, length.avg, tot_abund, mu_C, mu_N) 



## ellip dataframe for the mean landscape (visualization only)
ellip.mean = ellip.data(100, 1, n_species)
# Cluster to create heights
clo = makeCluster(workers) # Create workers
clusterEvalQ(clo, c(library(mvtnorm), library(SIBER), library(dplyr))) # Set up packages on the cluster
clusterExport(clo, list = c("posterior", "ellipfunc")) # #Export dependants
ellip.mean$string = parRapply(clo,ellip.mean, # run the ellipfun across the coords
                              function(x) ellipfunc(xax = x[3],
                                                    yax = x[4],
                                                    post_n = x[2],
                                                    spp = x[1]))
stopCluster(clo) ## stop workers/cluster


# Filter the data frame string
ellip.mean.filtered = ellip.mean %>% 
  filter_ellip.data(., n_species) 

# Facet wrap heat maps to see individual ellipses
ellip.mean.filtered %>% 
  #filter(spp != 1) %>%
  left_join(legend,by = c("spp" = "group")) %>% 
  #filter(spp == 1) %>%
  ggplot(aes(x= xax, y = yax, col = log10(string))) +
  geom_point() + 
  facet_wrap(~species.early) 

## Filtered heatmap with landscape for all species
ellip.mean.filtered %>%
  group_by(xax, yax) %>%
  summarize(vol = sum(string)) %>%
  #filter(xax < 3) %>%
  ggplot(aes(x = xax, y = yax, fill = log(vol))) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() + 
  labs(fill = "Z height") + 
  ggtitle("Nick's Circles")



## Trying out backtracing

ellip.mean.filtered %>%
  mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
         yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  group_by(xax, yax) %>%
  summarize(vol = sum(string)) %>%
  filter(xax < 3) %>%
  ggplot(aes(x = xax, y = yax, fill = log(vol))) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() + 
  labs(fill = "Z height") +
  xlab("d13C") +
  ylab("d15N")



## Full landscape --------------------------------------

# Ellip coordinates ------------

#posterior = posterior.save
posterior = posterior.early %>% 
  
  group_by(species, post) %>%
  summarize(Sigma_1_1 = mean(Sigma_1_1), 
            Sigma_2_1 = mean(Sigma_2_1), 
            Sigma_1_2 = mean(Sigma_1_2),
            Sigma_2_2 = mean(Sigma_2_2),
            length.avg = mean(length.avg), 
            tot_abund = mean(tot_abund), 
            mu_C = mean(mu_1),
            mu_N = mean(mu_2)) 


posterior = posterior.early

ellip = ellip.data(xy_length, n.posts, posterior$species%>%unique()%>%length())

# Cluster to create heights
clo = makeCluster(workers) # set up cluster
clusterEvalQ(clo, c(library(mvtnorm), library(SIBER), library(dplyr))) # export libraries to the cluster
clusterExport(clo, list = c("posterior", "ellipfunc")) # export the functions and data frames to the cluster
ellip$string = parRapply(clo,ellip, ## Iterate through ellip dataframe with the ellipfunc
                         function(x) ellipfunc(xax = x[3],
                                               yax = x[4],
                                               post_n = x[2],
                                               spp = x[1]))
stopCluster(clo) # Stop workers


## Filtering the ellip data

filtered_data= filter_ellip.data(ellip, n_species)

ellip.filtered = filtered_data %>%
  mutate(csv = 1)

## Saving ellip early and ellip late to compare
#ellip.early = ellip.filtered



#ellip.late = ellip.filtered
## Single lake -----------





ellip = rbind(ellip.late %>% mutate(csv = 2), ellip.early) %>%
  left_join(legend, by = c( "spp" ="group")) %>%
  mutate(code = case_when(csv == 1 ~ species.early, 
                          csv == 2 ~ species.late)) %>%
  select(-species.early, -species.late)




# Volume metrics and boxplot ---------------------------------
# Total volume of community 
total_vols = ellip %>%
  group_by(code, post, csv) %>%
  summarise(total_vol = sum(string)) 


## Total volume of each species individually
total_vols %>%
  ggplot(aes(x = code, y = total_vol, col = csv %>% as.factor())) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Total Volume") +
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_log10() +
  theme_minimal(base_size = 14) +
  labs(col = "Period") + 
  scale_color_discrete(labels = c("1" = "Early", "2" = "Late"))


# Total volume of the entire community


total_vols %>% 
  ungroup() %>%
  group_by(post, csv) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = csv %>% as.factor(), y = total_vol, col = as.factor(csv))) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  scale_x_discrete("Period",labels = c("1" = "Early", "2" = "Late")) +
  ylab("Total Volume") +
  theme(legend.position = "none")



# Peak N and C Densities ------------------------------
## This has been back transformed



axis_slices = ellip %>% 
  group_by(post, xax, yax, csv) %>%
  summarize(total_volume = sum(string)) %>%
  left_join(backtrace) %>% 
  mutate(xax = xax * sd_C + mean_C, 
         yax = yax * sd_N + mean_N) %>%
  group_by(post,csv) %>%
  slice_max(total_volume) %>%
  select(post, xax, yax, csv) %>% 
  pivot_longer(c("xax", "yax"), names_to = "Type", values_to = "Values")




axis_slices %>% 
  mutate(Type = str_replace(Type, "xax", "D13C"),
         Type = str_replace(Type, "yax", "D15N")) %>%
  ggplot(aes(x = as.factor(csv), y = Values, col = as.factor(csv))) + 
  geom_boxplot()+ 
  ylab("Axis value") + 
  facet_wrap(~ Type, scales = "free") + 
  scale_x_discrete("Period", labels = c("1" = "Early", "2" = "Late")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")




# Distance from center of ellipse to peak density ----------------------------
#### Needs to be modified in the distances function




# set up matrix for filtering within the for loop
z.mod = ellip %>% 
  group_by(post, xax, yax, csv) %>%
  summarize(total_volume = sum(string)) %>%
  left_join(backtrace) %>% 
  mutate(xax = xax * sd_C + mean_C, 
         yax = yax * sd_N + mean_N) %>%
  ungroup() %>%
  select(total_volume, xax, yax, post, csv) 

## for loop
dist_list = list()
dist_list1 = list()
for(i in 1:length(unique(z.mod$csv))){
  for(h in 1:length(unique(z.mod$post))){
    
    
    
    z = z.mod %>%
      filter(csv == i, post == h) %>%
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

dist_matrix = data.frame(value = c(unlist(dist_list[[1]]), unlist(dist_list[[2]]))) %>% 
  mutate(metric = rep(c("mean", "sd"), n.posts*2),
         post = rep(rep(1:n.posts, each = 2), 2), 
         csv= rep(c(1,2), each = n.posts * 2)) 

dist_matrix %>% 
  ggplot(aes(x = as.factor(csv), y = value, col = as.factor(csv))) + 
  geom_boxplot() + 
  facet_wrap(~metric) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Value") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Early", "2" = "Late"))

## Distance between next highest peak? ----------------------------

posterior.early %>% select(post, species, mu_1, mu_2)


ellip %>% 
  group_by(post, spp, csv) %>%
  slice_max(string) %>%
  left_join(backtrace)  %>%
  mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
         yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  ungroup() %>% 
  group_by(csv, post) %>%
  arrange(post, -string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values") %>%
  ggplot(aes(x = as.factor(csv), y = values, col = as.factor(csv))) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Value") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Early", "2" = "Late")) +
  facet_wrap(~metric, scales = "free")


## Demonstration of what it is doing


dist.graph.dat = ellip %>% 
  group_by(post, spp, csv) %>%
  slice_max(string) %>%
  left_join(backtrace)  %>%
  mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
         yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  ungroup() %>% 
  group_by(csv, post) %>%
  arrange(post, -string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  filter(csv == 2, post == 1) %>%
  ungroup()



ellip %>% filter(post == 1, csv == 2)  %>%
  left_join(backtrace)  %>%
  mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
         yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  group_by(csv, post, xax, yax) %>%
  summarize(sum = sum(string)) %>%
  ungroup() %>%
  ggplot(aes(x = xax, y = yax)) + 
  geom_tile(aes(fill = sum), alpha = .8) + 
  geom_point(data = dist.graph.dat, mapping = aes(x = xax, y = yax)) + 
  geom_path(data = dist.graph.dat, mapping = aes(x = xax, y = yax, col =string )) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  ylab("d15N") + xlab("d13C")






## Standard deviation of community -------------------------

z_comm %>% 
  group_by(post,csv) %>%
  summarize(sd = sd(tot)) %>%
  ggplot(aes(x = as.factor(csv), y = sd, col = as.factor(csv))) +
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  xlab("Period") + 
  ylab("Standard Deviation") + 
  scale_x_discrete(labels = c("1" = "Early", "2" = "Late"))


## Standard deviation of individual species 

ellip %>%
  group_by(post, csv, code) %>%
  summarise(mean_sd = sd(string)) %>%
  ggplot(aes(x = code, y = log(mean_sd))) + 
  geom_boxplot() + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust=1)) + 
  ylab("Ellipse Standard Deviation") + 
  xlab("Species") + 
  facet_wrap(~csv) 



## Rugosity -----------------------------------------





cell_size = (abs(cord_min) + abs(cord_max)) / 100

backtrace = backtrace %>% 
  mutate(cell_size = (abs(cord_min) + abs(cord_max)) / 100) %>%
  mutate(cell_size.C = cell_size*sd_C, 
         cell.size.N = cell_size * sd_N)






# Calculate the planar area for each cell


rug = list()
rug1 = list()

for(i in 1:length(unique(z.mod$csv))){
  # Define cell dimensions
  cell_height <- backtrace$cell.size.N[i]
  cell_width = backtrace$cell_size.C[i] # Example height
  planar_area <- cell_height * cell_width
  for(h in 1:length(unique(z.mod$post))){
    
    elevation_matrix = z.mod %>%
      filter(csv == i, post == h) %>%
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
                      post = rep(c(1:n.posts), rep = length(unique(z.mod$csv))),
                      csv = rep(unique(z.mod$csv), each = n.posts))

rugosity %>%
  ggplot(aes(x = as.factor(csv), y = values, col = as.factor(csv))) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Rugosity") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "Early", "2" = "Late")) 

