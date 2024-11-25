## Baselines and SIMMR Script
# Libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
`%nin%` = Negate(`%in%`)


## Script

JML_baselines = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Species %in% c("Snails", "Periphyton", "Zoops","Heptageniidae")) %>%
  separate(Date, into = c("MONTH","DAY","YEAR"))%>%  
  mutate(WATER = "LML")  %>%
  rename(TAXON = Species, D13C = C, D15N = N, GROUP = Group)  %>%
  mutate(GROUP = case_when(TAXON == "Periphyton" ~ "algae", TAXON == "Snails" ~ "SNAIL", TAXON == "Zoops" ~ "ZOOP",
                           TAXON == "Heptageniidae" ~ "INSECT")) %>% 
  select(WATER, YEAR, TAXON, GROUP, D13C, D15N)


baselines.all =  read.csv("Data/SI_MEASUREMENT.csv") %>%
  filter(GROUP %in% c("PERI", "ZOOP", "INSECT", "SNAIL"), 
         str_detect(ISO_YSAMP_N, "LML"), 
         
         SAMPLE_TYPE == "TISSUE") %>% 
  separate(ISO_YSAMP_N, into = c("GEAR", "WATER","YEAR","SAMP")) %>%
  select(WATER, YEAR, TAXON, GROUP, D13C, D15N) %>%
  mutate(YEAR = case_when(GROUP %in% c("ZOOP", "PERI") ~ "2023",
                          GROUP %nin% c("ZOOP","PERI") ~ YEAR)) %>%
  rbind(JML_baselines) %>%
  mutate(D13C = as.numeric(D13C),
         D15N = as.numeric(D15N)) %>%
  mutate(community = case_when(YEAR < 2005 ~ "early", 
                               YEAR > 2019 ~ "late")) 
## Snails only
snail = baselines.all %>% filter(GROUP %in%  c("SNAIL")) %>%
  group_by(GROUP, community) %>%
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N), sd_C = sd(D13C), sd_N = sd(D15N))
snail.detailed  = baselines.all %>% filter(GROUP %in%  c("SNAIL")) %>%
  group_by(GROUP, community)
## Mayflies only
insect = baselines.all %>%
  mutate(TAXON = toupper(TAXON)) %>% filter(GROUP %in% c("INSECT"), TAXON %in% c("HEPTAGENIIDAE","EPHEMEROPTERA")) %>%
  group_by(GROUP, community) %>% 
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N), sd_C = sd(D13C), sd_N = sd(D15N))
## Zooplankton only
zoop = baselines.all %>% filter(GROUP == "ZOOP")  %>%
  group_by(GROUP, community) %>%
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N), sd_C = sd(D13C), sd_N = sd(D15N))

## Periphyton only

peri = baselines.all %>% filter(GROUP %in% c("algae", "PERI")) %>%
  group_by(GROUP, community) %>%
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N, na.rm = T), sd_C = sd(D13C), sd_N = sd(D15N, na.rm =T))
peri

## bind together all three potential baselines
baselines = rbind(snail, insect, zoop) %>%
  na.omit()

save(baselines, file = "Data/IsotopeComparisons/baselines.RData")


# Creating the baselines that go into SIMMR
baseline_simmr = baselines %>% filter(GROUP %in% c("INSECT", "ZOOP"))

## Creating a baselin for SIMMR that averages the insect and snail baseline

baseline_simmr = baselines.all %>% 
  mutate(hab = case_when(GROUP %in% c("SNAIL", "INSECT")~"B", GROUP == "ZOOP" ~ "P") ,
         TAXON = toupper(TAXON)) %>%
  filter(GROUP %in% c("SNAIL", "ZOOP") | TAXON %in% c("HEPTAGENIIDAE","EPHEMEROPTERA")) %>%
  ungroup() %>%
  select(-GROUP) %>%
  rename(GROUP = hab) %>%
  select(GROUP, everything()) %>%
  group_by(GROUP, community) %>%
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N, na.rm = T), sd_C = sd(D13C), sd_N = sd(D15N, na.rm =T))

## Save the baseline_simmr frame

save(baseline_simmr, file = "Data/IsotopeComparisons/baseline_simmr.RData")

#baseline_simmr[4,4] = 1.1 ## Artificially lowering the d15N of zooplankton to see how it influences the diet contributions


load(file = "Data/IsotopeComparisons/IsotopeData.RData")

data_simmr = full %>%
  select(C, N, period, Species)  %>%
  rename(iso1 = C, 
         iso2 = N, 
         community = period, 
         group = Species) 

library(simmr)

CI_list = list()


for(i in 1:length(data_simmr$community %>% unique())){
  
  data_s = data_simmr %>% filter(community == (data_simmr$community %>% unique())[i]) %>% ungroup() %>% as.data.frame() %>%
    mutate(iso1 = as.numeric(iso1), 
           iso2 = as.numeric(iso2)) %>%
    na.omit()
  
  baseline_s = baseline_simmr %>% filter(community == (data_simmr$community %>% unique())[i]) %>% ungroup() %>% as.data.frame()
  
  S_N = c("Benthic", "Pelagic")
  
  
  simmr = simmr_load(
    mixtures = data_s[,1:2] %>% as.matrix(),
    source_names = S_N %>% as.matrix(),
    source_means = baseline_s[,3:4],
    source_sds = baseline_s[,5:6], 
    group = data_s[,4] %>% as.matrix()
  )
  
  title_graph = ggplot() + ggtitle(paste(i))
  
  print(title_graph)
  plot(simmr)
  
  
  simmr_out <- simmr_mcmc(simmr)
  #summary(simmr_out,
  # type = "diagnostics",
  # group = 1
  # )
  #posterior_predictive(simmr_out, group = 1)
  
  #prior_viz(simmr_out)
  #plot(simmr_out, type = "histogram")
  
  compare_groups(simmr_out,
                 groups = 1:length(unique((data_s$group))),
                 
                 source_name = "Benthic")
  
  
  
  CI_list[[i]] = summary(simmr_out, 
                         group = 1:length(unique(data_s$group)), 
                         type = "quantiles")
  
  
}


save(CI_list, file = "Data/IsotopeComparisons/SimmrOutput.RData")




