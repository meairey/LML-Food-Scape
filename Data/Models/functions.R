### Functions -----------

## Ellip function

xax =1;yax=1
spp = 9; post_n =1; community_n =2
ellipfunc = function( xax,yax,spp,post_n, community_n){
  
  dat.ellipfunc = posterior %>%
    ungroup() %>% filter(species == spp,## Filter out the species 
                                       post == post_n, 
                                       community == community_n) ## Filter out the posterior draw
  
  ## For each posterior draw for each species, calculate the height of the ellipse
  z =as.vector(dat.ellipfunc$tot_abund) * # abundance
    (as.vector(dat.ellipfunc$mass.avg)^.75) * # average size
    dmvnorm(c(xax,yax),
            mean = as.vector(c(dat.ellipfunc$mu_C, dat.ellipfunc$mu_N)), ## I changed this from mean to mu but am not sure if its an issue with which one of the packages I load in first
            sigma = matrix(as.vector(c(dat.ellipfunc$Sigma_1_1,
                                       dat.ellipfunc$Sigma_2_1, 
                                       dat.ellipfunc$Sigma_1_2,
                                       dat.ellipfunc$Sigma_2_2)),2,2))
  return(z)
}

## Multiple CSV filter
filter_ellip.data = function(x,spp, community){
  
  ## Define pop density and body size of species 
  for(z in 1:community){
    for(i in 1:n.posts){
      for(h in 1:spp){
        max = as.numeric(x %>% 
                           filter(post == i, spp == h, community == z) %>%
                           summarise(max = max(string)))
        remove = max*.1
        if((h == 1) & (i == 1)){
          ellip = x %>% filter(spp == h,
                               post == i,
                               community == z,
                               string >= remove)
        }else{
          new = x %>% filter(spp == h, 
                             post == i,
                             community == z,
                             string >= remove) 
          ellip = rbind(ellip, new)
        }
      }
    }
    
    if(z == 1){
      ellip.1 = ellip
    } else if (z == 2){
      ellip.2 = ellip
    } else {
      ellip.3 = ellip
    } 
    
    
  }
  
  
  
  ellip.full = rbind(ellip.1, ellip.2, ellip.3)
    return(ellip.full)
}

xy_length = 10

post_draws = 1
spp = 10
community_N = 3
## Fitting ellipses
## I think this will be done best by looping through species?
ellip.data = function(xy_length, post_draws, spp, community_N){
  
  ## Define pop density and body size of species 
  
  
  ## Gather information about species in dataset
  # Num. of species
  
  
  
  ## Create a series of x, y coordinate points to plot across. 
  ### Note - this may need to be adjusted based on ellipses locations
  x.p = seq(cord_min,cord_max,length.out=xy_length); y.p = x.p 
  
  # posterior draws and ellipse size
  n.posts <- post_draws; p.ell <- 0.90
  
  
  
  
  df.ellip = data.frame(
    spp = rep(rep(1:spp, each = length(x.p)^2*n.posts), community_N),
    post = rep(rep(rep(1:n.posts, each =length(x.p)^2),spp),  community_N), 
    xax = rep(rep(rep(rep(x.p, length(x.p)),n.posts),spp),  community_N),
    yax = rep(rep(rep(rep(y.p,each = length(x.p)),n.posts),spp),  community_N),
    community = rep(1:community_N, each = (length(x.p)^2)*spp*n.posts)
    )
  
  
  
  return(df.ellip)
}




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



find_local_maxima <- function(z) {  ## function that finds the local maximum
  maxima <- matrix(FALSE, nrow(z), ncol(z))
  for (i in 2:(nrow(z) - 1)) {
    for (j in 2:(ncol(z) - 1)) {
      if (z[i, j] > z[i - 1, j] && z[i, j] > z[i + 1, j] &&
          z[i, j] > z[i, j - 1] && z[i, j] > z[i, j + 1]) {
        maxima[i, j] <- TRUE
      }
    }
  }
  maxima
}
