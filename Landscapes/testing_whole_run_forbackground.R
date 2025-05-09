## List to fill
# Set random seed for reproducibility
set.seed(123)
## Load in libraries

library(raster)

library(MASS) # for mvrnorm

library(tidyverse)
## Source and Functions

## function

## Calculate surface area ---------------------


calculate_surface_area <- function(matrix) {
  nrow <- nrow(matrix)
  ncol <- ncol(matrix)
  
  # Initialize surface area
  surface_area <- 0
  
  # Calculate horizontal cell faces (area = cell_width * cell_height)
  horizontal_face_area <- cell_width * cell_height * (nrow - 1) * (ncol - 1)
  
  # Compute surface area from elevation differences
  for (i in 1:(nrow - 1)) {
    for (j in 1:(ncol - 1)) {
      dz_dx <- matrix[i + 1, j] - matrix[i, j]
      dz_dy <- matrix[i, j + 1] - matrix[i, j]
      
      # Diagonal distances (hypotenuse for each cell's side)
      hypotenuse_x <- sqrt(cell_width^2 + dz_dx^2)
      hypotenuse_y <- sqrt(cell_height^2 + dz_dy^2)
      
      # Add the areas of vertical faces
      surface_area <- surface_area + hypotenuse_x * cell_height + hypotenuse_y * cell_width
    }
  }
  
  # Add horizontal faces area
  surface_area <- surface_area + horizontal_face_area
  
  return(surface_area)
}


## Legends 



## Legend
species.early = c("BB", "CC","CS","MM","PS","SMB","SS","WS") # 8 long

species.pre = c("BB", "CC","CS","PS","SMB","SS","WS") #7 long

species.late =  c( "CC","CS","LT","MM","PS","RS","SMB","SS","WS") # 9 long



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

mean.legend = legend$species %>% unique() %>% as.data.frame() %>% 
  na.omit() %>% 
  filter(. != "NA") %>%
  mutate(group = as.character(as.numeric(as.factor(.)))) %>%
  arrange(group) %>%
  rename("species" = ".")

## Posterior

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
  












points.list = list()
total.vol.list = list() ## Added this in to make a list of total volumes
nhsp.list = list()
rugosity.list = list()
n_points = 500000
#n_points = 100
grid_size = .01
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape")
load("Data/VaryingIsotopesData/pointsgen.RData")

## script setup
n_species = 10
n.posts = 20
cord_min = -7; cord_max = 7
# Resolution
xy_length = 50



#points.gen = points.gen %>% filter(post == 1)
  
  
  
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





## Run the loop and time it ------------------------------------------------


points.gen = points.gen %>% filter(post %in% c(1:2))

execution_time <- system.time({
  
for(h in 1:length(unique(points.gen$post))){
  
points.gen.post = points.gen %>% filter(post == h)
  
points.frame = NULL
  
  for(i in 1:dim(points.gen.post)[1]){
    ## Create the subset of the posterior your're going to use
    subset = points.gen.post %>% 
    filter(lookup == points.gen.post$lookup[i])
  
    
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
    
  
  ## Calculate metrics per posterior draw

  grid_size = .01 ## Grid size 
  
  ## Points frame to work with per posterior draw
  points.frame <- bind_rows(points.list)%>% 
    separate(lookup, into = c("community", "group", "post")) %>%
    left_join(posterior %>%
                select(community, group, post, mass.avg, tot_abund, species) %>%
                mutate(community = as.character(community), 
                       group = as.character(group), 
                       post = as.character(post)), 
              by = c("community", "group", "post")) %>%
    #select(-species) %>%
    mutate(z_string = z_string * tot_abund * (mass.avg ^ .75)) %>%
    ungroup()
  
points.frame %>% select(community, group,species) %>%
  unique()
  
  
  ## Total Volume Metric - can be split out later
  total.vol.list[[h]] = points.frame %>%
    group_by(species,group, post, community) %>%
    summarise(total_vol = sum(z_string))# %>%
    #left_join(mean.legend)
  
  
    
  
  
  ## Next highest species peak
  nhsp.list[[h]] = points.frame %>% 
  mutate(species = group) %>%
  group_by(post, species, community) %>%
  slice_max(z_string) %>%
  ungroup() %>% 
  group_by(community, post) %>%
  arrange(post, -z_string) %>%
  mutate(dist = sqrt((x1_grid - lag(x1_grid))^2 + (x2_grid - lag(x2_grid))^2)) %>%
  select(dist,everything()) %>%
  summarize(mean = mean(dist, na.rm = T), sd = sd(dist, na.rm = T)) %>%
  pivot_longer(mean:sd, names_to = "metric", values_to = "values")
  
  ## Rugosity
  z.mod = points.frame %>% 
  group_by(post, x1_grid, x2_grid, community) %>%
  summarize(total_volume = sum(z_string)) %>%
  ungroup() %>%
  select(total_volume, x1_grid, x2_grid, post, community) %>%
      rename(xax = x1_grid, 
             yax = x2_grid)
  
  
    
  cell_size = (abs(cord_min) + abs(cord_max)) / xy_length
  
  cell_size = grid_size
  
  
  
  # Calculate the planar area for each cell



  rug = list()

  rug.community = list()
  for(j in 1:length(unique(z.mod$community))){

    cell_height = cell_size
    cell_width = cell_size
    planar_area <- cell_height * cell_width
    
      
      elevation_matrix = z.mod %>%
        filter(community ==j, post ==h) %>%
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
      
      rug.community[[j]] = data.frame(community = j, post = h, rug = rugosity)
      
   
    }
  
    rugosity.list[[h]] = bind_rows(rug.community)

  
  
  }

  
total.vol = bind_rows(total.vol.list) ## dataframe with volume estimates for each species
nhsp = bind_rows(nhsp.list)
rugosity.summary = bind_rows(rugosity.list)



})


print(execution_time)

## 482 seconds for 2 posterior draws 

((execution_time/h * 1000)/60/60)[3] /24


save(total.vol,file = "Data/VaryingIsotopesData/total.vol.frame.RData")
save(nhsp, file = "Data/VaryingIsotopesData/nhsp.RData")
save(rugosity.summary, file = "Data/VaryingIsotopesData/rugosity.summary.RData")

total.vol %>% group_by(community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = community, y = total_vol)) + 
  geom_boxplot()

nhsp %>% 
  ggplot(aes(x = community, y = values)) + 
  geom_boxplot() + 
  facet_wrap(~metric)

rugosity.summary %>%
  ggplot(aes(x = as.factor(community), y = rug)) + 
  geom_boxplot()

total.vol %>% 
  filter(species == "SMB")


