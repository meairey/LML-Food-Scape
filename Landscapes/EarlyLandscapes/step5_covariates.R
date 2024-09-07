



## Filtering out the year range
temp = temp.full %>%
  filter(YEAR %in% (year_min:year_max), YEAR != 2002)


## Habitat data 

hab = (shoreline_length %>% filter(Water == "LML"))$Habitat

step_a$z_covariate <- sample(c(0, 1), nrow(step_a), replace = TRUE)