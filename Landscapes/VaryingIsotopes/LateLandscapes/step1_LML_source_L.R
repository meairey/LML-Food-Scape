library(tidyverse)
## Load in data frames
fish = read.csv("../AFRP/Data/FISH_MEASUREMENT_LML.csv") 
sample = read.csv("../AFRP/MA2276_Code/Data/FISH_SAMPLE_2022.csv")
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
species.late = c( "CC","CS","LT","MM","PS","RS","SMB","SS","WS")

## Set parameters for the model run
year_max = 2023
year_min = 2019
species = species.early
## Graphing 
legend = data.frame(group = c(1:9), 
                    species = species.late )


## Backtracing -------------
## Early Period -- 
backtrace.early = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Group %in% "Fish") %>%
  select(Species, C, N) %>% 
  arrange(Species) %>% 
  mutate(spp = as.numeric(as.factor(Species)),
         community = 1, 
         C = as.numeric(C),
         N = as.numeric(N)) %>%
  select(community, spp, C, N) %>%
  rename("iso_1" = "C",
         "iso_2" = "N") %>%
  group_by(spp) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(community) %>%
  summarize(mean_C = mean((iso_1)), sd_C = sd(iso_1),
            mean_N = mean((iso_2)), sd_N = sd(iso_2)) 
## Late period -- 
backtrace.late = read.csv("Data/SI_MEASUREMENT.csv") %>%
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  filter(WATER == "LML",
         YEAR > 2018, 
         GROUP == "FISH",
         SAMPLE_TYPE == "TISSUE") %>%
  select(TAXON, D13C, D15N) %>% 
  unique() %>% 
  arrange(TAXON) %>%
  rename("iso_1" = "D13C",
         "iso_2" = "D15N",
         "group" = TAXON) %>%
  group_by(group) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  na.omit()  %>%
  as.data.frame() %>%
  mutate(community = 2) %>%
  group_by(community) %>%
  summarize(mean_C = mean((iso_1)), sd_C = sd(iso_1),
            mean_N = mean((iso_2)), sd_N = sd(iso_2)) 


back.trace = rbind(backtrace.early, backtrace.late)
