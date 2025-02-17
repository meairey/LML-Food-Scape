## Trying to speed up the calculations

library(mvtnorm)
library(MASS)

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
species.early = c("BB", "CC","CS","MM","PS","SMB","SS","WS") # 8 long
species.late =  c( "CC","CS","LT","MM","PS","RS","SMB","SS","WS") # 9 long
species.pre = c("BB", "CC","CS","PS","SMB","SS","WS") #7 long
legend = data.frame(group = as.character(c(1:9)), 
                    C  = species.late, 
                    A = c(species.pre, rep("NA", 2)), 
                    B = c(species.early, rep("NA",1))) %>%
  pivot_longer(2:4, names_to = "community", values_to = "species") %>%
  arrange(community) %>%
  mutate(community = as.numeric(as.factor(community)))


legend
#---------------------------------------------------------------------------

load("Data/VaryingIsotopesData/posterior_early.csv")
load("Data/VaryingIsotopesData/posterior_late.csv")
load("Data/VaryingIsotopesData/posterior_pre.csv")


## Comparison of posterior parameter distributions
posterior.pre = posterior.pre %>%
  rename("group" = "species") %>% 
  mutate(community = 1) %>%
  left_join(legend, by = c("community", "group")) %>%
  ungroup() 


posterior.early = posterior.early %>%
  rename("group" = "species") %>% 
  mutate(community = 2) %>%
  left_join(legend, by = c("community", "group")) 

posterior.late = posterior.late %>%
  rename("group" = "species") %>% 
  mutate(community = 3) %>%
  left_join(legend, by = c("community", "group")) 

posterior = rbind(posterior.pre, posterior.early, posterior.late) %>%
  mutate(group = as.numeric(as.factor(species)))


# Function to calculate boxplot percentiles
bp.pctiles = function (x, probs = c(0.025, 0.25, 0.5, 0.75, .975)) {
  r <- quantile(x, probs = probs, na.rm = TRUE)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


posterior %>% 
  mutate(group = (species)) %>%
  ggplot(aes(x = as.factor(community), y = tot_abund)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", aes(width=0.4)) +
  facet_wrap(~species, scales = "free_y") +
  xlab("Period") +
  ylab("Total Abundance") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  scale_y_log10() 


##

posterior %>% 
  mutate(group =species) %>%
  #left_join(legend) %>%
  ggplot(aes(x = as.factor(community), y = mass.avg)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", aes(width=0.4)) +
  #scale_y_log10() +
  facet_wrap(~species, scales = "free_y") +
  xlab("Period") +
  ylab("Average Mass (g)") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 


posterior %>% 
  mutate(biomass = tot_abund * mass.avg) %>%
  mutate(group = species) %>% 
  ggplot(aes(x = as.factor(community), y = biomass)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", aes(width=0.4)) +
  #scale_y_log10() +
  facet_wrap(~species, scales = "free_y") +
  xlab("Period") +
  ylab("Biomass") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 


## Community biomass
posterior %>% 
  mutate(biomass = tot_abund * mass.avg) %>%
  group_by(community, post) %>%
  mutate(comm_biomass = sum(biomass)) %>% 
  ggplot(aes(x = as.factor(community), y = comm_biomass)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", aes(width=0.4)) + 
  xlab("Period") +
  ylab("Biomass") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 



# Ellipse generation ----------

## Coord limits


workers = 5



### Mean posterior heatmap for visualization ------

## script setup
n_species = 10
n.posts = 3000        
cord_min = -7; cord_max = 7
# Resolution
xy_length =100


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

## ellip dataframe for the mean landscape (visualization only)
ellip.mean = ellip.data(xy_length, 1, 
                        length(unique(posterior$species)),
                        length(unique(posterior$community))) 

p = posterior %>%
  as.data.frame() %>%
  select(species, community, post , everything()) %>%
  mutate(ma.tot = (mass.avg^.75) * tot_abund) %>%
  rename("spp" = "species")




combined <- left_join(p, ellip.mean, by = c("spp", "community", "post")) %>%
  rowwise() %>% # Ensure row-wise computation
  mutate(
    z = dmvnorm(
      c(xax, yax), # Ensure xax and yax are numeric and conform to the sigma dimensions
      mean = c(mu_C, mu_N), # Ensure mu_C and mu_N are numeric
      sigma = matrix(
        c(Sigma_1_1, Sigma_1_2, Sigma_2_1, Sigma_2_2), # Build the covariance matrix
        nrow = 2, byrow = TRUE
      )
    )
  ) %>%
  ungroup() %>%
  select(spp, community, post, xax, yax, z) %>%
  rename("string" = "z") %>%
  group_by(community, spp, post)%>%
  mutate( string = string / max(string) * 0.4 ) %>%
  filter(string >= quantile(string, 0.4, na.rm = T)) 



  
  
## Trying 


dog = posterior %>%
  slice(rep(1:n(), each = xy_length^2)) %>%
  mutate(xax = ellip$xax, yax = ellip$yax)

ellip$cat = 1
  
#  library(MASS)  # For 'ginv' if covariance matrix inversion fails
library(car)
  # Define the parameters
mu_C <- 0      # Mean x-coordinate
mu_N <- 0      # Mean y-coordinate
cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix
  
# Example points (x, y)
points <- data.frame(x = c(0, 1, 2, 20), y = c(0, 1, 2, 100))
  
# Inverse covariance matrix
inv_cov_matrix <- solve(cov_matrix)  # Or use ginv() for pseudo-inverse
  
# Compute Mahalanobis distances
points$mahal_dist <- apply(points, 1, function(point) {
  diff <- c(point["x"] - mu_C, point["y"] - mu_N)
  sqrt(t(diff) %*% inv_cov_matrix %*% diff)
})
  
# Critical value for 95% confidence
critical_value <- sqrt(qchisq(0.95, df = 2))  # Square root for distance
  
# Check if points fall within the ellipse
points$inside_ellipse <- points$mahal_dist <= critical_value
  
print(points)


library(ellipse)
  
# Parameters
center <- c(0, 0)  # Mean vector (mu_C, mu_N)
cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix
  
# Generate ellipse points at 95% confidence level
ellipse_data <- data.frame(ellipse(cov_matrix, centre = center, level = 0.95))
  
# Plot the ellipse
library(ggplot2)
ggplot() +
    geom_path(data = ellipse_data, aes(x = x, y = y)) +
    coord_equal() +
    theme_minimal() + 
    geom_point(aes(x = points$x, points$y))




###


  
## Testing if this breaks the full version --------------------
  ## Trying out backtracing
  ## script setup
n_species = 10
n.posts = 5
cord_min = -7; cord_max = 7
# Resolution
xy_length = 100
  

  
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



p = posterior %>%
  as.data.frame() %>%
  select(species, community, post , everything()) %>%
  mutate(ma.tot = (mass.avg^.75) * tot_abund) %>%
  rename("spp" = "species")




execution_time <- system.time({

combined <- left_join(ellip,p, by = c("spp", "community", "post")) %>%
  rowwise() %>% # Ensure row-wise computation
  mutate(
    z = dmvnorm(
      c(xax, yax), # Ensure xax and yax are numeric and conform to the sigma dimensions
      mean = c(mu_C, mu_N), # Ensure mu_C and mu_N are numeric
      sigma = matrix(
        c(Sigma_1_1, Sigma_1_2, Sigma_2_1, Sigma_2_2), # Build the covariance matrix
        nrow = 2, byrow = TRUE
      )
    )
  ) %>%
  ungroup() %>%
  select(spp, community, post, xax, yax, z) %>%
  rename("string" = "z") %>%
  group_by(community, spp, post)%>%
  mutate( string = string / max(string) * 0.4 ) %>%
  filter(string >= quantile(string, 0.4, na.rm = T)) 



# Start timing a code block

  
})
save(combined, file = "Data/save")
print(execution_time)
combined$post %>% unique()

## Chopping this up into chunks
chunks = 500 ## 45 seconds per posterior draw... so 500 draws == 6.25 hrs. So 2000 draws == 

chunks = data.frame(start = c(1, 3), end = c(2, 4))

list.combined = list()

for(i in 1:length(chunks[,1])){
  

  p = posterior %>%
   # filter(post > chunks[i,1] & post <= chunks[i,2]) %>%
  as.data.frame() %>%
  select(species, community, post , everything()) %>%
  mutate(ma.tot = (mass.avg^.75) * tot_abund) %>%
  rename("spp" = "species")
  
  e = ellip %>%
    filter(post > chunks[i,1] & post <= chunks[i,2])
  
  
  
  combined <- left_join(e,p, by = c("spp", "community", "post")) %>%
  rowwise() %>% # Ensure row-wise computation
  mutate(
    z = dmvnorm(
      c(xax, yax), # Ensure xax and yax are numeric and conform to the sigma dimensions
      mean = c(mu_C, mu_N), # Ensure mu_C and mu_N are numeric
      sigma = matrix(
        c(Sigma_1_1, Sigma_1_2, Sigma_2_1, Sigma_2_2), # Build the covariance matrix
        nrow = 2, byrow = TRUE
      )
    )
  ) %>%
  ungroup() %>%
  select(spp, community, post, xax, yax, z) %>%
  rename("string" = "z") %>%
  group_by(community, spp, post)%>%
  mutate( string = string / max(string) * 0.4 ) %>%
  filter(string >= quantile(string, 0.4, na.rm = T)) 
  
  list.combined[[i]] = combined
  
}


list.combined




## --------------------------

> print(execution_time)
user  system elapsed 
392.50   11.95  409.90 --> 7min for 100 posterior draws
.07/min per posterior draw --> 2000draws = 2000 mins




filter_ellip.data = function(x, spp, community) {
  
  # Create a list to store results for each community
  ellips = vector("list", community)
  
  # Calculate max values and thresholds for filtering in one go
  thresholds = x %>%
    group_by(community, post, spp) %>%
    summarise(max_string = max(string), .groups = 'drop') %>%
    mutate(remove = max_string * 0.1)
  
  # Filter the main data based on the calculated thresholds
  ellip_full = x %>%
    inner_join(thresholds, by = c("community", "post", "spp")) %>%
    filter(string >= remove) %>%
    select(-max_string, -remove) # Remove unnecessary columns
  
  return(ellip_full)
}

filtered = filter_ellip.data(combined,10,3)
dim(filtered)
dim(combined)

filtered %>% group_by(post, community,xax, yax) %>% 
  summarize(tot = sum(string)) %>%
  ggplot(aes(x = xax, y = yax, col = tot)) + 
  geom_point() + 
  facet_wrap(~community) + 
  scale_color_viridis_c()

combined
  ungroup() %>%
  group_by(community, post, xax, yax) %>%
  summarize(tot = sum(string, na.rm = T)) %>%
  ggplot(aes(x = xax, y = yax, col = tot)) + 
  geom_point() + 
  facet_wrap(~community) + 
  scale_color_viridis_c()


