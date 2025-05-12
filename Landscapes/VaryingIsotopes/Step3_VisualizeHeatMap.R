# Set random seed for reproducibility
set.seed(123)
## Load in libraries
library(snow)
library(raster)
library(viridis)
library(MASS) # for mvrnorm
library(mvtnorm)
library(tidyverse)
library(gridExtra)
library(wesanderson)
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

species.pre = c("BB", "CC","CS","PS","SMB","SS","WS") #7 long

species.late =  c( "CC","CS","LT","MM","PS","RS","SMB","SS","WS") # 9 long

## Backtrasforming values for plotting

load("Data/VaryingIsotopesData/EarlyData/backtrace_early.RData")
load("Data/VaryingIsotopesData/LateData/BacktraceTrace.RData")
load("Data/VaryingIsotopesData/PreData/backtrace_pre.RData")
backtrace.pre = backtrace.early
backtrace.pre$community =1

back.trace = rbind(backtrace.pre,backtrace.early, backtrace.late)

## Legend for joining together the posteriors - includes group/community combinations
legend = data.frame(group = as.character(c(1:9)), 
                    A = c(species.pre, rep("NA", 2)),
                    B = c(species.early, rep("NA",1)),
                    C  = species.late
                     
                    )%>%
  pivot_longer(2:4, names_to = "community", values_to = "species") %>%
  arrange(community) %>%
  mutate(community = as.numeric(as.factor(community)))
  
legend

bayes_cri <- function(x) {
  data.frame(
    y = mean(x),  # Mean of posterior samples
    ymin = quantile(x, 0.025),  # 2.5% quantile (Lower bound)
    ymax = quantile(x, 0.975)   # 97.5% quantile (Upper bound)
  )
}



mean.legend = legend$species %>% unique() %>% as.data.frame() %>% 
  na.omit() %>% 
  filter(. != "NA") %>%
  mutate(group = as.character(as.numeric(as.factor(.)))) %>%
  arrange(group) %>%
  rename("species" = ".")
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

posterior %>% 
  select(community, species, group) %>%
  unique() %>%
  print(n = 100)
  

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

## Species specific biomass

posterior %>% 
  mutate(biomass = tot_abund * mass.avg) %>%
  group_by(community,species, post) %>%
  mutate(comm_biomass = sum(biomass)) %>% 
  ggplot(aes(x = as.factor(community), y = comm_biomass)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", aes(width=0.4)) + 
  xlab("Period") +
  ylab("Biomass") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  facet_wrap(~species, scales= "free_y") +
  scale_y_log10()


### Mean posterior heatmap for visualization ------

## script setup
n_species = 10
n.posts = 3000


## New type of heat map with averages from posterior

posterior.mean = rbind(posterior.pre, posterior.early, posterior.late) %>% 
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

## Trying to use posteriors from community 2 model in pre-initiation model

pm1 = posterior.mean %>% 
  filter(community == 1) %>% 
  select(species, community, mass.avg, tot_abund, post) %>%
  mutate(Sigma_1_1 = (posterior.mean %>% filter(community == 2))$Sigma_1_1) %>%
  mutate(Sigma_2_1 = (posterior.mean %>% filter(community == 2))$Sigma_2_1) %>%
  mutate(Sigma_1_2 = (posterior.mean %>% filter(community == 2))$Sigma_1_2) %>%
  mutate(Sigma_2_2 = (posterior.mean %>% filter(community == 2))$Sigma_2_2) %>%
  mutate(mu_C = (posterior.mean %>% filter(community == 2))$mu_C) %>%
  mutate(mu_N = (posterior.mean %>% filter(community == 2))$mu_N) 


posterior.mean = rbind(posterior.mean %>% filter(community!= 1), pm1) %>%
  as.data.frame() %>%
  mutate(across(c(Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2, mu_C, mu_N), ~ case_when(is.na(tot_abund) ~ NA_real_, TRUE ~ .)))


  
## Ellipse - should contain 40% of ellipse data 

confidence_level <- 0.40 # Confidence level for the ellipse
threshold <- qchisq(confidence_level, df = 2) # df = 2 for 2D
# Height function
z_func <- function(X1, X2, mu, sigma) {
  z <- dmvnorm(
    c(X1, X2),  # Combine X1 and X2 into a numeric vector
    mean = mu,  # Mean vector
    sigma = sigma  # Covariance matrix
  )
  
  return(z)
}


## Generating points and heights -----------------------------
# Find the threshold for the given confidence level

n_points.mean <- 1000000 # Number of points to generate


## List to fill

points.list.mean = list()

points.gen.mean =  posterior.mean %>% 
  unite("lookup", c(community, species, post), remove = F) %>%
  arrange(community, species, post) %>%
  na.omit()


execution_time <- system.time({
for(i in 1:dim(points.gen.mean)[1]){
  ## Create the subset of the posterior your're going to use
  subset = points.gen.mean %>% 
  filter(lookup == points.gen.mean$lookup[i])

  mu = c(subset$mu_C, subset$mu_N) # Ensure mu_C and mu_N are numeric
  sigma = matrix(
        c(subset$Sigma_1_1, subset$Sigma_1_2, subset$Sigma_2_1, subset$Sigma_2_2), # Build the covariance matrix
        nrow = 2, byrow = TRUE)
  
  inv_sigma <- solve(sigma)
  
  # Generate random points from the multivariate normal distribution

  points <- mvrnorm(n = n_points.mean, mu = mu, Sigma = sigma)

  
  
  # Calculate the Mahalanobis distance for each point
  mahalanobis_dist <- apply(points, 1, function(pt) t(pt - mu) %*% inv_sigma %*% (pt - mu))
  
  # Keep only points within the ellipse
  points_in_ellipse <- points[mahalanobis_dist <= threshold, ]  %>%
    as.data.frame() %>%
  mutate(x1_grid = round(V1 / grid_size) * grid_size,
         x2_grid = round(V2 / grid_size) * grid_size,
         mean_area = ellipse_area) %>%
    unique()
  
  
  
    # Compute ellipse area
  eigen_values <- eigen(sigma)$values
  a <- sqrt(eigen_values[1]) * threshold
  b <- sqrt(eigen_values[2]) * threshold
  ellipse_area <- pi * a * b
  
  ## Calculate z string
  
  z <- apply(points_in_ellipse, 1, function(point) z_func(point[1], point[2], mu, sigma))
  
  
  
  ellipse_points <- data.frame(points_in_ellipse) %>%
    mutate(lookup  = subset$lookup)   %>%
   mutate(z_string = z, mean_area = ellipse_area)
  
  
  points.list.mean[[i]] = ellipse_points
  

  
}
  
## Snap points to a grid and them summarize across that grid
grid_size <- .01

points.frame.mean <- bind_rows(points.list.mean)%>% 
  separate(lookup, into = c("community", "group", "post")) %>%
  left_join(posterior.mean %>%
              select(community, species, post, mass.avg, tot_abund) %>%
              mutate(community = as.character(community), 
                     group = as.character(species), 
                     post = as.character(post)), 
            by = c("community", "group", "post")) %>%
  select(-species) %>%
  rename(xax = X1, yax = X2) %>%
  mutate(x1_grid = round(xax / grid_size) * grid_size,
         x2_grid = round(yax / grid_size) * grid_size) %>% ## Snaps to grid  
  mutate(z_string.mod = z_string * tot_abund * (mass.avg ^ .75)) %>%
  ungroup() 
})


print(execution_time)
#save(points.frame.mean, file = "Data/VaryingIsotopesData/points_frame.mean.RData")

load(file =  "Data/VaryingIsotopesData/points_frame.mean.RData")
## Looking at what the total volume of each species is
v = points.frame.mean  %>% 
  group_by(community, post, group, mean_area, x1_grid, x2_grid) %>%
  summarize(z.mean = mean(z_string), 
            z.mean.mod = mean(z_string.mod))

v %>% ungroup() %>%
  group_by(community, group, mean_area) %>%
  mutate(cel.vol = z.mean * grid_size) %>%
  summarize(total_volume = sum(cel.vol)) %>%
  ggplot(aes(x = mean_area, y = total_volume)) + 
  geom_point() 


v %>% group_by(community, post, group, mean_area) %>%
  summarize(count = n()) %>%
  mutate(ration = count / mean_area) %>%
  ggplot(aes(x = mean_area, y = ration)) + 
  geom_point()



points.frame.mean  %>% 
  group_by(community, post, group, mean_area) %>%
  summarize(count = n()) %>% 
  ungroup() %>%
  ggplot(aes(x = count)) +
  geom_histogram()

 ## Heat map


v %>% 
  group_by(community, post, group, mean_area) %>%
  summarise(total_volume = sum(z.mean)) %>%
  mutate(total_volume = total_volume * grid_size) %>%
  ggplot(aes(x = group, y = total_volume, col = group)) + 
  geom_point() +
  facet_wrap(~community)




## Filtered heatmap with landscape for all species
v %>%
  na.omit() %>%
  #left_join(back.trace) %>%
  #mutate(xax = xax*sd_C + mean_C, yax = yax*sd_N + mean_N) %>%
  group_by(community, post, x1_grid, x2_grid) %>%
  summarize(z_string = sum(z.mean)) %>%
  ungroup() %>%
  ggplot(aes(x = x1_grid, y = x2_grid, fill = (z_string))) +
  #geom_point() +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() + 
  labs(fill = "Z") +
  xlab("d13C") +
  ylab("d15N") + 
  facet_wrap(~community, labeller = labeller(community = c("1" = "Pre-Initiation", "2"= "Post-Initiation", "3" = "Modern Observation"))) 


# Facet wrap heat maps to see individual ellipses
points.frame.mean %>% 
  mutate(group = as.character(group)) %>%
  left_join(mean.legend,by = c("group" = "group")) %>% 
  na.omit() %>%
  ggplot(aes(x= x1_grid, y = x2_grid, col =as.factor(community))) +
  geom_point(alpha = .15) + 
  facet_wrap(~species) 