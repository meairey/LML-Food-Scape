library(tidyverse)
## Load in data frames
fish = read.csv("../2. clean_AFRPData/FISH_MEASUREMENT_LML.csv") 
sample = read.csv("../2. clean_AFRPData/FISH_SAMPLE_2024.csv")
sites = read.csv("../AFRP/MA2276_Code/Data/SITES.csv")
shoreline_length = read.csv("../AFRP/MA2276_Code/Data/BEFsites_LengthAndHabitat.csv")

## Values
epsilon <- 1e-6
n.iso = 2

## Covariates ---------

## Temperature data 

temp.full = read.csv("Data/1_AirTempMod_V3_MeanMonthlyTemp.csv")  %>%
  filter(MONTH %in% c(4:11)) %>%
  group_by(YEAR) %>%
  summarize(mean_y = mean(MN_MNTH_TMP_C, na.rm = T)) 


temp.full[20,1] = 2020
temp.full[20,2] = 12.8
temp.full[21,1] = 2021
temp.full[21,2] = 12.8


## Currently this is set up to be run twice. Once for the early period and once for the late period. I think this is to deal with some of the uncertainty of the smelt population in the mid 2000s, unless I come up with a better solution for modelling populations.

### LML Model Draft ---------------------------------------------------

## LML Data 
LML.data = read.csv("Data/LML_data_cpue.csv")


## What species to include in the early period
species.early = c("BB", "CC","CS","LT","MM","PS","RS","SMB","SS","WS")
## what species to include in the late period
#species.late = c("CC","CS","MM","PS","SMB","WS")  ## Going to try to simulate data for these species in the interim
species.late = c("BB", "CC","CS","LT","MM","PS","RS","SMB","SS","WS")

## Species pre
species.pre = c("BB", "CC","CS","LT","PS","SMB","SS","WS")

## Set parameters for the model run
year_max = 2000
year_min = 2000
species = species.pre
## Graphing 
legend = data.frame(group = c(1:8), 
                    species  = species.pre)

