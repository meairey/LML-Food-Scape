source("Landscapes/EarlyLandscapes/step1_LML_source.R")
## Setting up covariates for JAGS Early
n_years = year_max - year_min

## Filtering out the year range
temp = temp.full %>%
  filter(YEAR %in% (year_min:year_max), YEAR != 2002)
covariate_T = temp$mean_y %>% scale(, center = 0) %>% as.numeric()

save(file = "Data/EarlyData/covariate_temp.RData", covariate_T)

## Habitat data covariates 
covariate_hab = (shoreline_length %>% filter(Water == "LML"))$Habitat %>%
  as.factor() %>%
  as.numeric()
save(file = "Data/EarlyData/habs_cov.RData", covariate_hab)

n_sites = length(covariate_hab)

z_covariate <- sample(c(0, 1), n_species * n_years * n_sites, replace = TRUE)
save(file = "Data/EarlyData/z_covariate.RData", z_covariate)
