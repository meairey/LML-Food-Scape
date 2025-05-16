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
         





#### --------------- trying to change so that the isotopes are the same between community 1 and community 2

## Trying to use posteriors from community 2 model in pre-initiation model

p1 = posterior %>% 
  filter(community == 1) %>% 
  select(species, community, mass.avg, tot_abund, post) %>%
  mutate(Sigma_1_1 = (posterior %>% filter(community == 2))$Sigma_1_1) %>%
  mutate(Sigma_2_1 = (posterior %>% filter(community == 2))$Sigma_2_1) %>%
  mutate(Sigma_1_2 = (posterior %>% filter(community == 2))$Sigma_1_2) %>%
  mutate(Sigma_2_2 = (posterior %>% filter(community == 2))$Sigma_2_2) %>%
  mutate(mu_C = (posterior %>% filter(community == 2))$mu_C) %>%
  mutate(mu_N = (posterior %>% filter(community == 2))$mu_N) 


posterior = rbind(posterior %>% filter(community!= 1), p1) %>%
  as.data.frame() %>%
  mutate(across(c(Sigma_1_1, Sigma_2_1, Sigma_1_2, Sigma_2_2, mu_C, mu_N), ~ case_when(is.na(tot_abund) ~ NA_real_, TRUE ~ .)))

n.posts = 1000
points.gen =  posterior %>% 
  unite("lookup", c(community, species, post), remove = F) %>%
  arrange(community, species, post) %>%
  filter(post <= n.posts) %>%
  na.omit()## Run below loop over 2000 posteriors
save(points.gen, file = "Data/VaryingIsotopesData/pointsgen.RData")

###

## List to fill

points.list = list()
n_points = 500000
grid_size = .01
## Run the loop and time it
execution_time <- system.time({
for(i in 1:dim(points.gen)[1]){
  ## Create the subset of the posterior your're going to use
  subset = points.gen %>% 
  filter(lookup == points.gen$lookup[i])

  mu = c(subset$mu_C, subset$mu_N) # Ensure mu_C and mu_N are numeric
  sigma = matrix(
        c(subset$Sigma_1_1, subset$Sigma_1_2, subset$Sigma_2_1, subset$Sigma_2_2), # Build the covariance matrix
        nrow = 2, byrow = TRUE)
  
  inv_sigma <- solve(sigma)
  
  # Generate random points from the multivariate normal distribution

  points <- mvrnorm(n = n_points, mu = mu, Sigma = sigma)

  
  # Compute ellipse area
  eigen_values <- eigen(sigma)$values
  a <- sqrt(eigen_values[1]) * threshold
  b <- sqrt(eigen_values[2]) * threshold
  ellipse_area <- pi * a * b
  
  
  # Calculate the Mahalanobis distance for each point
  mahalanobis_dist <- apply(points, 1, function(pt) t(pt - mu) %*% inv_sigma %*% (pt - mu))
  
  # Keep only points within the ellipse
  points_in_ellipse <- points[mahalanobis_dist <= threshold, ] %>% 
    as.data.frame() %>%
  mutate(x1_grid = round(V1 / grid_size) * grid_size,
         x2_grid = round(V2 / grid_size) * grid_size) %>% 
    select(x1_grid, x2_grid) %>% unique() 
    
    
  
  z <- apply(points_in_ellipse, 1, function(point) z_func(point[1], point[2], mu, sigma))
  
  
  
  ellipse_points <- data.frame(points_in_ellipse) %>%
    mutate(lookup  = subset$lookup,
          mean_area = ellipse_area) %>%
   mutate(z_string = z)
  
  
  points.list[[i]] = ellipse_points 
  
  print(i/dim(points.gen)[1]*100)
  
}
  

grid_size = .01
points.frame <- bind_rows(points.list)%>% 
  separate(lookup, into = c("community", "group", "post")) %>%
  left_join(posterior %>%
              select(community, species, post, mass.avg, tot_abund) %>%
              mutate(community = as.character(community), 
                     group = as.character(species), 
                     post = as.character(post)), 
            by = c("community", "group", "post")) %>%
  select(-species) %>%
  mutate(z_string.mod = z_string * tot_abund * (mass.avg ^ .75)) %>%
  ungroup() 

})

print(execution_time)
## Just for reference
#i = 173 after 3.5hrs... so i = 1000

save(points.frame, file = "Data/VaryingIsotopesData/points_frame.RData")

## START HERE load points -----------------------------------------------
load(file = "Data/VaryingIsotopesData/points_frame.RData")

colnames(points.frame)
points.frame$post %>% unique()




## Checking if species have different volumes inherently
points.frame %>% 
  filter(community == 3, post == 1) %>%
  group_by(community, post, group) %>%
  summarise(total_volume = sum(z_string)) %>%
  ggplot(aes(x = group, y = total_volume)) + 
  geom_boxplot() + 
  facet_wrap(~community)

points.frame %>% 
  filter(community == 3, post == 1) %>%
  group_by(community, post, group) %>%
  summarise(total_volume = sum(z_string)) 

points.frame$post %>% unique() %>% length()


points.frame %>% 
  filter(community == 3, post == 1) %>%
  group_by(community, x1_grid, x2_grid) %>%
  summarize(total_z = sum(z_string)) %>%
  ggplot(aes(x = x1_grid, y = x2_grid, col = (total_z))) + 
  geom_point()



points.frame = points.frame %>%
  mutate(species = group) %>%
  mutate(z_string = z_string.mod)
  
# Volume metrics and boxplot ---------------------------------
# Total volume of community 
total_vols = points.frame %>%
  group_by(group, post, community) %>%
  summarise(total_vol = sum(z_string)) %>%
  left_join(mean.legend)


## Total volume of each species individually
species_vols= total_vols %>%
  as.data.frame() %>%
  complete(species, post, community) %>%
  replace_na(list(total_vol = 0)) 


#save(species_vols, file = "Data/VaryingIsotopesData/species_vols.RData")

species_Vols_graph = species_vols %>%
  ggplot(aes(x = community %>% as.factor(),
             y = total_vol +1,
             col = community %>% as.factor())) +
  #geom_boxplot(outlier.shape = NA, key_glyph = "rect") +
  ylab("Total Volume") +
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_log10() +
  theme_minimal(base_size = 14) +
  labs(col = "Period") +
  scale_x_discrete("Period",
                   labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.") ) +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none", 
        axis.title.x = element_blank()) +
  scale_color_manual("Community",
                    values = wes_palette("Darjeeling1", 
                                         type = "discrete", n = 3),
  labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) + 
  facet_wrap(~species, scales = "free_y") +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 

species_Vols_graph
save(species_Vols_graph, file = "Data/VaryingIsotopesData/species_vols_graph.RData")
load(file = "Data/VaryingIsotopesData/species_vols_graph.RData")


## Contemplating the incorporation of sensitive vs. tolerant natives in this conversation
summary.vol = species_vols %>% 
  filter(species == "SMB") %>%
  left_join(suc, by = c("species" = "CODE")) %>%
  group_by(suc, post, community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ungroup() %>%
  group_by(suc, community) %>%
  summarize(mean = quantiles_95(total_vol)[3], min = quantiles_95(total_vol)[2],
            max = quantiles_95(total_vol)[4])

species_spec.vol = species_vols %>% 
  filter(species %nin% c("SMB")) %>%
  left_join(suc, by = c("species" = "CODE")) %>%
  group_by(suc, species, community) %>%
  summarize(mean =  quantiles_95(total_vol)[3], min = quantiles_95(total_vol)[2],
            max = quantiles_95(total_vol)[4])


  
  #### Load in the legend
load("Data/IsotopeComparisons/simmr_legend.RData")

load("Data/Legend.col.RData") ## Has colors for plotting individuals species


legend_sci = read.csv("Data/IsotopeComparisons/legend.csv")




legend = legend %>% left_join(legend_sci, by = c("CODE" = "code"))


ggplot() +
  geom_bar(data = species_spec.vol, aes(x = suc, y = mean, fill = species), stat = "identity") +
  facet_wrap(~community, labeller =labeller(community = c("1" = "Pre-Removal", "2" = "Post-Initiation", "3" = "Modern Observation"))) +
  theme_minimal(base_size = 14) +
  scale_x_discrete("Response", labels = c("ne" = "Declined", 
                              "p" = "Recovered",
                              "na" = "None")) +
  scale_fill_manual("Species", values = (legend %>% filter(CODE %in% c("BB","CC","CS","LT","MM","PS","RS", "SS", "WS")))$color,
                    labels =(legend %>% filter(CODE %in% c("BB","CC","CS","LT","MM","PS","RS","SS", "WS")))$scientific ) +
 # geom_hline(data = summary.vol, aes( yintercept = mean), lty = 2) +
  ylab("Niche Volume") +
  theme(legend.position = "bottom")
  
  

  

summary.vol$mean[1] - summary.vol$mean[3]

cat = species_spec.vol %>%
  group_by(community) %>%
  summarize(total = sum(mean))

cat$total[1] - cat$total[3]

smb_loss = 19631487 
community_gain = 2097469

(community_gain / smb_loss) * 100

# Total volume of the entire community





total_vols_graph = total_vols %>% 
  ungroup() %>%
  group_by(post, community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = community %>% as.factor(),
             y = total_vol, col = as.factor(community))) +
 # geom_boxplot() +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  ylab("Total Volume") +
  xlab("Community") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  scale_color_manual("Community", values = wes_palette("Darjeeling1", type = "discrete", n = 3),
  labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 

  
save(total_vols_graph, file = "Data/VaryingIsotopesData/total_vols_graph.RData")
total_vols_graph

total_vols_native = total_vols %>% 
  filter(species %nin% c("SMB", "RS", "MM")) %>%
  ungroup() %>%
  group_by(post, community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = community %>% as.factor(),
             y = total_vol, col = as.factor(community))) +
 # geom_boxplot() +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  ylab("Total Volume") +
  xlab("Community") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  scale_color_manual("Community", values = wes_palette("Darjeeling1", type = "discrete", n = 3),
  labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 

total_vols_native




# Distance from center of ellipse to peak density ----------------------------
#### Needs to be modified in the distances function






## Distance between next highest peak? ----------------------------
# I dont think you want to transform these values back into isotopically meaningful ones... it would be misleading if the variation isn't standardized...

## Demonstration of what it is doing





## Next highest species peak ---------------

nhsp.graph = points.frame %>% 
  mutate(species = group) %>%
  group_by(post, species, community) %>%
  slice_max(z_string) %>%

  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -z_string) %>%
  mutate(dist = sqrt((x1_grid - lag(x1_grid))^2 + (x2_grid - lag(x2_grid))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values") %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  #geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("DNHP") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  facet_wrap(~metric, scales = "free",
             labeller = labeller("metric" = c("mean" = "Mean DNHP", "sd" = "SD DNHP"))) +

 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )   + 
  scale_color_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 3),
                     labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) 

nhsp.graph



save(nhsp.graph, file = "Data/VaryingIsotopesData/nhsp.graph.RData")

## Distance graph

dist.graph.dat = points.frame %>% 
  
  group_by(post, species, community) %>%
  slice_max(z_string) %>%
  rename(xax = x1_grid, 
         yax = x2_grid) %>%
  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -z_string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  filter(community == 2, post == 3) %>%
  ungroup() %>%
  mutate(index = c(1:8))



points.frame %>% filter(post == 3, community == 2)  %>%
  filter(group != 8) %>%
  rename(xax = x1_grid, 
         yax = x2_grid) %>%
  group_by(community, post, xax, yax) %>%
  #summarize(sum = sum(z_string)) %>%
  ungroup() %>%
  ggplot(aes(x = as.numeric(xax), y = as.numeric(yax))) + 
  geom_point(aes(col = (z_string)), alpha = .8) + 
  geom_point(data = dist.graph.dat[,], mapping = aes(x = xax, y = yax)) + 
  geom_path(data = dist.graph.dat[,], mapping = aes(x = xax, y = yax, col = -index ), size = 1) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  ylab("d15N") + xlab("d13C") +
  scale_fill_gradient(high = "white", low = "black")

## Checking credible intervals on the DNHP differences

points.frame %>% 
    rename(xax = x1_grid, 
         yax = x2_grid) %>%
  group_by(post, species, community) %>%
  slice_max(z_string) %>%

  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -z_string) %>%
  mutate(dist = sqrt((xax - lag(xax))^2 + (yax - lag(yax))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values") %>%
  ungroup() %>%
  pivot_wider(names_from = community, values_from = values) %>%
  mutate(diff_12 = `2` - `1`,
         diff_31 = `3` - `1`,
         diff_23 = `3` - `2`) %>%
  group_by(metric) %>%
  summarize(mean12 = mean(diff_12),
            low12 = quantile(diff_12, .025),
            up12= quantile(diff_12, .95),
            mean13 = mean(diff_31),
            low13 = quantile(diff_31, .025),
            up13 = quantile(diff_31, .95))



## Standard deviation of peak height


cat = points.frame %>% 
  ungroup() %>%
  filter(group != 8) %>%
  group_by(post,group, community) %>%
  slice_max(z_string) %>% 
  group_by(post, community) %>%
  summarize(max_height = max(z_string), 
            mean_height = mean(z_string), 
            sd_height = sd(z_string))

cat %>% 
  ggplot(aes(x = community, y = mean_height)) +

 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 



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

points.frame %>%
  group_by(post, community, species) %>%
  summarise(mean_sd = sd(z_string)) %>%
  ggplot(aes(x = species, y = log(mean_sd))) + 
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

cell_size = grid_size



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
      mutate_all(~ ifelse(is.na(.), 0, .)) %>% 
      select(-post, -community)
    
    # Calculate surface area
    surface_area <- calculate_surface_area(as.matrix(elevation_matrix))
    
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

#save(rugosity, file = "Data/VaryingIsotopesData/rugosity.RData")
load(file = "Data/VaryingIsotopesData/rugosity.RData" )
rugosity_graph = rugosity %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  #geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("Rugosity") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  scale_color_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 3),labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 

rugosity_graph



## What about a metric that looks at the standard deviation in peak height? is this the same as rugosity?








load(file = "Data/VaryingIsotopesData/simmr_vol.RData")




simmr_vol_graph = simmr_vol %>%
#  filter(species != "SMB") %>%
  mutate(weighted.avg = total_vol * value) %>% ## total volume per species/post/community * weighted average of that species' benthic contribution
  group_by(post, community) %>%
  mutate(total_web_vol = sum(total_vol)) %>% ## Total volume of the web
  ungroup() %>%
  group_by(post, community, total_web_vol) %>%
  summarize(wa = sum(weighted.avg, na.rm = T)) %>% ## Add together the weighted volume for each species
  ungroup() %>%
  mutate(benthic.contr = wa / total_web_vol) %>% ## divide it by the total web volume to get proportion/percent contribution
  ggplot(aes(x = as.factor(community), y = 100*benthic.contr, col = as.factor(community))) +
  stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )  +
  theme_minimal(base_size = 14) +
  ylab("Benthic Contribution (%)") +
  scale_x_discrete("Period",
                   labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.") ) +

  scale_color_manual("Community",
                    values = wes_palette("Darjeeling1", 
                                         type = "discrete", n = 3),
  labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  theme(legend.position = "non")

simmr_vol %>%
#  filter(species != "SMB") %>%
  mutate(weighted.avg = total_vol * value) %>%
  group_by(post, community) %>%
  mutate(total_web_vol = sum(total_vol)) %>%
  ungroup() %>%
  group_by(post, community, total_web_vol) %>%
  summarize(wa = sum(weighted.avg, na.rm = T)) %>%
  ungroup() %>%
  mutate(benthic.contr = wa / total_web_vol) %>% 
  select(-total_web_vol, -wa) %>%
  pivot_wider(names_from = community, values_from = benthic.contr) %>% 
  mutate(diff13 = `3` - `2`,
         diff12 = `2` - `1`,
         diff23 = `3` - `2`) %>%
  select(diff13, diff12, diff23) %>%
  reframe(mean13 = mean(diff13), low13 = quantile(diff13, .025), up13 = quantile(diff13, .975),
            mean12 = mean(diff12), low12 = quantile(diff12, .025), up12 = quantile(diff12, .975),
            mean23 = mean(diff23), low23 = quantile(diff23, .025), up23 = quantile(diff23), .975)


save(simmr_vol_graph, file = "Data/VaryingIsotopesData/simmr_vol_graph.RData")  


load(file = "Data/VaryingIsotopesData/simmr_vol_graph.RData")
simmr_vol_graph



### Grid arrange all figures -------------------

library(gridExtra)
grid.arrange(total_vols_graph, simmr_vol_graph,rugosity_graph, nhsp.graph, ncol = 2)

species_Vols_graph
rugosity_graph

rugosity_graph
nhsp.graph




### old Metrics ----------------------------

## Next highest point peak ----------------

z.mod %>% 
  group_by(community, post) %>%
  arrange(post, -total_volume) %>%
  mutate(dist = sqrt((xax- lag(xax))^2 + (yax - lag(yax))^2)) %>%
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

points.frame %>% filter(post == 1, community == 1)  %>%
  #left_join(back.trace)  %>%
  # mutate(xax = xax * backtrace.early$sd_C + backtrace.early$mean_C,
  #       yax = yax * backtrace.early$sd_N + backtrace.early$mean_N) %>%
  group_by(community, post, x1_grid, x2_grid) %>%
  summarize(sum = sum(z_string)) %>%
  ungroup() %>%
  ggplot(aes(x = x1_grid, y = x2_grid)) + 
  #geom_tile(aes(fill = sum), alpha = .8) + 
  geom_point(aes(col = sum)) +
  geom_point(data = dist.graph.dat[1:10000,], mapping = aes(x = xax, y = yax)) + 
  geom_path(data = dist.graph.dat[1:10000,], mapping = aes(x = xax, y = yax, col =string )) +
  theme_minimal(base_size = 14) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  ylab("d15N") + xlab("d13C")



# Peak N and C Densities ------------------------------
## This has been back transformed



axis_slices = points.frame %>% 
  group_by(post, x1_grid, x2_grid, community) %>%
  summarize(total_volume = sum(z_string)) %>%
  left_join(back.trace %>% 
              mutate(community = as.character(community))) %>% 
  mutate(xax = x1_grid * sd_C + mean_C, 
        yax = x2_grid * sd_N + mean_N) %>%
  group_by(post,community) %>%
  slice_max(total_volume) %>%
  select(post, xax, yax, community) %>% 
  pivot_longer(c("xax", "yax"), names_to = "Type", values_to = "Values")

(axis_slices %>% filter(Type == "xax", community == 1))$Values


axis_slice_graph = axis_slices %>% 
  mutate(Type = str_replace(Type, "xax", "D13C"),
         Type = str_replace(Type, "yax", "D15N")) %>%
  ggplot(aes(x = as.factor(community), y = Values, col = as.factor(community))) + 
 # geom_boxplot()+ 
  ylab("PPT") + 
  facet_wrap(~ Type, scales = "free")  +
  labs(col = "Period") +
  scale_x_discrete("Community",labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  scale_color_manual( values = wes_palette("Darjeeling1", 
                                         type = "discrete", n = 3),
                    labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 


axis_slice_graph 


# set up matrix for filtering within the for loop
## This wont be backtransformed because I don't want to influence the dist with the natural variation in baseline
z.mod = points.frame %>% 
  group_by(post, x1_grid, x2_grid, community) %>%
  summarize(total_volume = sum(z_string)) %>%
  ungroup() %>%
  select(total_volume, x1_grid, x2_grid, post, community) %>%
      rename(xax = x1_grid, 
             yax = x2_grid)

## for loop
dist_list = list()
dist_list1 = list()
options(expressions = 5000) 
for(i in 1:length(unique(z.mod$community))){
  for(h in 1:length(unique(z.mod$post))){
    
    
    
    z = z.mod  %>%
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

## Distance between all local maxima in the landscape

dist_matrix %>% 
  ggplot(aes(x = as.factor(community), y = value, fill = as.factor(community))) + 
  geom_boxplot() + 
  facet_wrap(~metric) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") +
  ylab("Distance between maxima") + xlab("Period")  +
  labs(col = "Period") +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  scale_fill_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 3),labels = c("1" = "Pre","2" = "Early", "3" = "Late")) 




