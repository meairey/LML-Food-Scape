## List to fill
# Set random seed for reproducibility
set.seed(123)
## Load in libraries
library(snow)
library(raster)
library(viridis)
library(MASS) # for mvrnorm
library(mvtnorm)
library(tidyverse)


points.list = list()
n_points = 500000
grid_size = .01
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape")
load("Data/VaryingIsotopesData/pointsgen.RData")


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
## Just for reference
#i = 173 after 3.5hrs... so i = 1000

save(points.frame, file = "Data/VaryingIsotopesData/points_frame.RData")