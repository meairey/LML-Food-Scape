library(tidyverse)
library(wesanderson)
`%nin%` = Negate(`%in%`)

## Loading in measurement information
fish_lml = read.csv("../AFRP/Data/FISH_MEASUREMENT_LML.csv") %>%
  filter(FISH_N != "") %>%
  select(YSAMP_N:WEIGHT) %>% unique()

lml_sample = read.csv("Data/FISH_SAMPLE_2024.csv")

lml_sample %>% filter(YSAMP_N == "GLN.LML.2023.003")

## Checking dataframe
## Access Database
data.acc = read.csv("../1.clean_isotope/SI_MEASUREMENT.csv") %>%
  filter(corrected %nin% c("delete duplicate","duplicate delete")) %>%
  select(-corrected)

data.acc %>% filter(TAXON == "NRD") %>%
  filter(grepl( "LML",ISO_FISH_N))


sample.clean = read.csv("../1.clean_isotope/SI_SAMPLE.csv")
da.acc = left_join(data.acc, sample.clean, by = "ISO_YSAMP_N") %>%
  filter(#GROUP == "FISH", 
         str_detect(ISO_YSAMP_N, "LML"), 
         
        # GROUP %nin% c("PHYTO", "PERI")
         ) %>%
  select(ISO_YSAMP_N, ISO_FISH_N, ITEM_N, NUM, D13C, D15N, CATEGORY, GROUP, TAXON, SAMPLE_TYPE, WATER, YEAR, SITE_N)
         



## Missing samples from JML File

mi.sa = read.csv("Data/IsotopeComparisons/Clean/missing_sample.csv") %>% unique()
mi.me = read.csv("Data/IsotopeComparisons/Clean/missing_measurement.csv")
miss.dat = left_join(mi.me, mi.sa, by = "ISO_YSAMP_N") %>% filter(YEAR < 2005) %>%
  select(ISO_YSAMP_N, ISO_FISH_N, ITEM_N, NUM, D13C, D15N, CATEGORY, GROUP, TAXON, SAMPLE_TYPE, WATER, YEAR, SITE_N)




## Load in more recent data
recent.data = read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  filter(grepl("LML",ISO_YSAMP_N)| ISO_YSAMP_N == "CHECK" | ISO_YSAMP_N == "UPDATE")
recent.sample = read.csv("../1.clean_isotope/isotope_sample.csv") %>%
  filter(WATER == "LML")

rec.dat = left_join(recent.data, recent.sample, by = "ISO_YSAMP_N") %>%
  left_join(fish_lml, by = c("ISO_FISH_N" = "FISH_N")) %>%
  select(ISO_YSAMP_N, ISO_FISH_N, ITEM_N, NUM, D13C, D15N, CATEGORY, GROUP, TAXON, SAMPLE_TYPE, WATER, YEAR, SITE_N) %>%
  filter(WATER == "LML")

rec.dat %>% filter(GROUP == "PERI") %>% 
  filter(WATER == "LML")

rec.dat %>% 
  filter(TAXON == "SS") %>%
  left_join(fish_lml, by = c("ISO_FISH_N" = "FISH_N"))

rec.dat %>% filter(TAXON == "ST")

## Updated samples that need ysamp in the isotope_sample sheet


recent.data %>% filter(ISO_YSAMP_N %nin% recent.sample$ISO_YSAMP_N) %>%
  select(ISO_YSAMP_N) %>%
  unique()

## This contains a subset of fish that are missing YSAMP info and thus not getting included?


updating.samples = recent.data %>% 
  filter(ISO_YSAMP_N %in% c("UPDATE", "CHECK", "")) %>%
  left_join(fish_lml, by = c("ISO_FISH_N" = "FISH_N")) %>%
  select(YSAMP_N, everything()) %>%
  left_join(lml_sample) 
  
samples.updated = updating.samples %>% 
  select(YSAMP_N, ISO_YSAMP_N, YEAR, SITE_N, DAY_N, DSAMP_N) %>%
  unique() %>%
  arrange(YEAR)


#write.csv(updating.samples, "updatingmeasurement.csv")

#write.csv(samples.updated, "updatingsamples.csv")

new.samps = read.csv("updatingsamples.csv")

new.samps %>% 
  group_by(YSAMP_N) %>%
  summarize(n = n())

cat = updating.samples %>% 
  select(-ISO_YSAMP_N) %>% 
    
    left_join(new.samps, by = "YSAMP_N")

write.csv(cat, "cat.csv")
## Adding in the additional fish from Kim 1/21/2025
ad.sa = read.csv("Data/IsotopeComparisons/Clean/sample_additional_fish.csv") %>%
  unique() %>%
  select( -YSAMP_N)

ad.ma = read.csv("Data/IsotopeComparisons/Clean/updated.measurement.csv") %>%
  filter(GROUP == "FISH")

## Pulling together the whole data
ad.dat = left_join(ad.ma, ad.sa, by = "ISO_YSAMP_N") %>%
  left_join(fish_lml, by = c("ISO_FISH_N" = "FISH_N"))  %>%
  select(ISO_YSAMP_N, ISO_FISH_N, ITEM_N, NUM, D13C, D15N, CATEGORY, GROUP, TAXON, SAMPLE_TYPE, WATER, YEAR, SITE_N) %>%
  filter(WATER == "LML")


## Full measurement file

measurement = rbind(da.acc, miss.dat, rec.dat, ad.dat)
## I think the clean updated data sheet in the 1. clean istope folder contains the additional fish I should check though
measurement = rbind(da.acc, miss.dat, rec.dat)

measurement$GROUP %>% unique()


measurement %>% filter(GROUP == "PERI", WATER == "LML")


measurement %>%
  filter(GROUP == "FISH") %>%
  ggplot(aes(x = SITE_N,y = TAXON)) + 
  geom_point(alpha = .1) + 
  facet_wrap(~YEAR)

## Get lengths out of JML file

jml.lengths = read.csv("Data/IsotopeComparisons/Clean/JML_weight_data.csv") %>%
  separate(Length, into = c("A", "B", "C","D","E","F","G","H","I")) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  mutate(Sum = rowSums(across(A:I, ~ as.numeric(.)), na.rm = TRUE)) %>%
  mutate(avg.length = Sum / Number) %>%
  select(ISO_FISH_N, avg.length, Number)







measurement %>% group_by(TAXON) %>%
  filter(YEAR < 2005) %>%
  summarize(count = n())






#Set up data ------------


data = left_join(measurement, fish_lml, by = c("ISO_FISH_N" = "FISH_N")) %>%
  left_join(jml.lengths, by = c("ITEM_N" = "ISO_FISH_N")) %>%
  mutate(LENGTH = case_when(is.na(avg.length) == F ~ avg.length, is.na(avg.length) == T ~ LENGTH)) %>%
  select(-avg.length) %>% 
  filter(is.na(D15N)==F & 
           is.na(D13C) == F & 
           WATER %in% c("LML"), 
         GROUP == "FISH")  %>%
  mutate(D13C = as.numeric(D13C), 
         D15N = as.numeric(D15N)) %>% 
  mutate(period = case_when(YEAR < 2006 ~ "early" , YEAR > 2018 ~ "late"))



data %>%
  group_by(TAXON, period) %>%
  summarize(n = n()) %>%
  print(n = 100)



save(data,  file= "Data/IsotopeComparisons/IsotopeData.RData")
full = data
#write.csv( full,"LML_iso_corrected.csv", row.names = FALSE) ## Old probably delete

## JML Baselines

baselines =  measurement %>%
  filter(
         GROUP %nin% c("FISH", "CRAYFISH"), 
         WATER == "LML") %>%
  select(WATER, YEAR, GROUP, TAXON, D13C, D15N ) %>%
  na.omit()


save(baselines, file = "Data/IsotopeComparisons/baselines.RData")






## Table summary of samples

load(file= "Data/IsotopeComparisons/IsotopeData.RData")



## Trying a new format

data %>% 
  mutate(period = case_when(YEAR < 2006 ~ "early" , YEAR > 2018 ~ "late")) %>%
  ggplot(aes(x = LENGTH, y = D13C, col = period)) + 
  geom_smooth(method = "lm", se = F) + 
  geom_point() + 
  facet_wrap(~TAXON, scales = "free_x")

data %>% 
  mutate(period = case_when(YEAR < 2006 ~ "early" , YEAR > 2018 ~ "late")) %>%
  ggplot(aes(x = LENGTH, y = D15N, col = period)) + 
  geom_smooth(method = "lm", se = F) + 
  geom_point() + 
  facet_wrap(~TAXON, scales = "free_x")



### Looking at differences in LT stable isotopes by capture location


LT = read.csv("Data/IsotopeComparisons/JML.Data.Master.csv") %>%
  filter(Species == "LT") 

LT %>% ggplot(aes(x = as.numeric(C), y = as.numeric(N), col = Area)) + 
  geom_point() + 
  stat_ellipse(level = .4)


measurement %>% filter(GROUP == "INSECT") %>%
  filter(YEAR %in% c(2021, 2023)) %>%
  ggplot(aes(x= as.numeric(D13C), y = D15N, col = SITE_N)) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~YEAR)


## Adding in the additional fish from Kim 1/21/2025



         

ad.sa = read.csv("Data/IsotopeComparisons/Clean/sample_additional_fish.csv") %>%
  unique() %>%
  select( -YSAMP_N)

ad.ma = read.csv("Data/IsotopeComparisons/Clean/updated.measurement.csv")


  
## habitat characteristics 
hab = read.csv("Data/IsotopeComparisons/Clean/BEFsites_LengthAndHabitat.csv")


## Pulling together the whole data
ad.dat = left_join(ad.ma  %>%
  filter(GROUP == "FISH"), ad.sa, by = "ISO_YSAMP_N") %>%
  left_join(fish_lml, by = c("ITEM_N" = "FISH_N")) %>% 
  left_join(hab)

ad.dat %>% select(ISO_FISH_N)

## Summary of how many sites for each habitat were sampled
ad.dat %>% 
  select(SITE_N, Habitat) %>%
  unique() %>%
  group_by(Habitat) %>%
  summarize(n())

## Summary of how many fish per each site were sampled
ad.dat %>% group_by(SITE_N) %>%
  filter(TAXON == "SMB") %>%
  filter(is.na(SITE_N) != T) %>%
  summarize(n()) %>%
  print(n = 100)
  
## Small SMB niches across habitats
ad.dat %>%
  filter(LENGTH < 100, 
         TAXON == "SMB", 
         SITE_N != "BEF.LML.005") %>%
  ggplot(aes(x = D13C, y = D15N, col = Habitat, size = LENGTH)) +
  geom_point() + 
  stat_ellipse(level = .4) 

## Large SMB niches across habitats
ad.dat %>%
  filter(LENGTH > 100, 
         TAXON == "SMB", 
         SITE_N != "BEF.LML.005") %>%
  ggplot(aes(x = D13C, y = D15N, col = Habitat)) +
  geom_point() + 
  stat_ellipse(level = .4) 

## All SMB across habitats
ad.dat %>%
  filter( 
         TAXON == "SMB", 
         SITE_N != "BEF.LML.005") %>%
  ggplot(aes(x = D13C, y = D15N, col = Habitat)) +
  geom_point() + 
  stat_ellipse(level = .4) 

## SMB across sites within habitats
ad.dat %>%
  filter(LENGTH < 100, 
         TAXON == "SMB", 
         #SITE_N != "BEF.LML.005"
         ) %>%
  ggplot(aes(x = D13C, y = D15N, col = SITE_N)) +
  geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~Habitat)


## Looking at isotopes across locations
ad.dat %>% 
  filter(TAXON == "SMB") %>%
  ggplot(aes(x = SITE_N, y = D15N, col = LENGTH)) + 
  geom_point()


ad.dat %>% 
    filter(TAXON == "SMB") %>%
  ggplot(aes(x = SITE_N, y =  D13C)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90)) 


## Looking at trophic shifts within site 5

ad.dat %>%
  #filter(SITE_N == "BEF.LML.005") %>%
  filter(TAXON == "SMB") %>%
  ggplot(aes(x = LENGTH, y = D13C, col = SITE_N)) + 
  geom_point() + 
  geom_smooth(method = lm, se = F)

ad.dat %>%
  #filter(SITE_N == "BEF.LML.005") %>%
  ggplot(aes(x = LENGTH, y = D15N, col = SITE_N)) + 
  geom_point() + 
  geom_smooth(method = lm, se = F)

ad.dat %>%
  ggplot(aes(x = D13C, y = D15N, col = TAXON)) + 
  geom_point() + 
  stat_ellipse(level = .4) + 
  facet_wrap(~SITE_N)

## Different sized fish

## It really looks like fish smaller than 200 mm look like one thing and then the larger fish are different
## So - ive made the breaks at 200
ad.dat %>% 
  filter(TAXON == "SMB") %>%
  mutate(TAXON = case_when(LENGTH < 100 ~ "SMB1",
                           #LENGTH >= 100 & LENGTH < 150 ~ "SMB2", 
                           LENGTH >=100 & LENGTH < 200 ~ "SMB3", 
                           LENGTH >= 200 ~ "SMB4")) %>%
  filter(is.na(TAXON) == F) %>%
  ggplot(aes(x = D13C, y = D15N, col = TAXON)) + 
  geom_point() + 
  stat_ellipse(level = .4) + 
  theme_minimal(base_size = 14) +
  scale_color_manual("Size Group", values = wes_palette("Darjeeling1", type = "discrete", n = 3), labels = c("< 100 mm" , "< 100 & < 200", "> 200")) +
  theme(legend.position = "top")
  
ad.dat %>% 
  filter(TAXON == "SMB") %>%
  count()




## Cleaning up and adding in Matt's data to the previous data sheets that I created (ad.ma and ad.sa)

## Looking at matt's data


## Pulling together the whole data


invert.taxa = read.csv("../reimagined-invertebrates/Data/CSVs/taxon_frame.csv")

ad.dat.inv  = left_join(ad.ma  %>%
  filter(CATEGORY == "INVERT"), ad.sa, by = "ISO_YSAMP_N") %>%
  left_join(hab)  %>%
    mutate(TAXON = tolower(TAXON)) %>%
  left_join(invert.taxa, by = c("TAXON")) 



ad.dat.inv %>% select(TAXON) %>%
  mutate(TAXON = tolower(TAXON)) %>%
  left_join(invert.taxa, by = c("TAXON")) %>% 
  filter(is.na(FAMILY)==T) %>% unique()


invert.taxa %>%
  filter(TAXON == "sphaeriidae") 
  


ad.dat.inv %>%
  group_by(FAMILY) %>%
  mutate(count = n()) %>%
  filter(count > 5) %>%
  ggplot(aes(x = D13C, y = D15N, col = SITE_N)) + 
  geom_point() + 
  stat_ellipse(level = .4) + 
  facet_wrap(~FAMILY)

## Which taxa found at most sites

ad.dat.inv %>% 
  select(SITE_N, FAMILY) %>%
  unique() %>%
  group_by(FAMILY) %>%
  summarize(total_sites = n()) %>%
  filter(total_sites > 10)

## Taxa at > 10 sites
ad.dat.inv %>%
  filter(FAMILY %in% c("cambaridae", "ephemereidae", "heptageniidae", "talitridae")) %>%
  ggplot(aes(x = D13C, y = D15N, col = SITE_N)) + 
  geom_point() + 
  stat_ellipse(level = .4) + 
  facet_wrap(~FAMILY)

