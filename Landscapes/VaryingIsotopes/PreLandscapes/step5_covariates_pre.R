source("Landscapes/VaryingIsotopes/PreLandscapes/step1_LML_source_pre.R")
## Setting up covariates for JAGS Early
n_years = 1
n_species = length(species)




## New Model -----------------------------------------------------




## Ice Off covariate -- new model 
ice_off = sample %>% left_join(read.csv("Data/ICE_OUT.csv"),by = c("YEAR" = "Year")) %>%
  filter(YEAR > 1997,
         GEAR == "BEF",
         GEAR_CODE == "NAF", 
         WATER == "LML", 
         MONTH < 7) %>%
  select(WATER, DATE_COL, MONTH, DAY_N, YEAR, Ice.out.date) %>%
  unique() %>%
  mutate(dif = DAY_N - Ice.out.date) %>%
  group_by(YEAR) %>%
  summarize(avg.day.ice = mean(dif) ) %>% 
  print(n = 100)
ice_off[24,2] = 30
ice_off[24,1] = 2022


ice_off =as.numeric((ice_off %>% filter(YEAR %in% c(year_min:year_max)))$avg.day.ice)
save(ice_off, file = "Data/VaryingIsotopesData/PreData/ice_off_pre.RData")



## Habitat variables
## New model
habs = shoreline_length %>% filter(Water == "LML") %>% arrange(SITE_N) %>%
  select(SITE_N,Habitat)%>%
  mutate(sub = substr(Habitat, 1, 1),
         wood = substr(Habitat, 2, 2)) %>%
  mutate(wood = as.numeric(as.factor(wood)),
         sub  = as.numeric(as.factor(sub))) 

save(habs, file = "Data/VaryingIsotopesData/PreData/habs.RData")






### Old Model --------------------------------------------

## Filtering out the year range
temp = temp.full %>%
  filter(YEAR %in% (year_min:year_max), YEAR != 2002)
covariate_T = temp$mean_y %>% scale(, center = 0) %>% as.numeric()

save(file = "Data/PreData/covariate_temp.RData", covariate_T)

## Habitat data covariates 
covariate_hab = (shoreline_length %>% filter(Water == "LML"))$Habitat %>%
  as.factor() %>%
  as.numeric()
save(file = "Data/PreData/habs_cov.RData", covariate_hab)

n_sites = length(covariate_hab)

z_covariate <- sample(c(0, 1), n_species * n_years * n_sites, replace = TRUE)
save(file = "Data/PreData/z_covariate.RData", z_covariate)


