## Data frame generation
sample = read.csv("Data/SI_SAMPLE.csv")

measurement = read.csv("Data/SI_MEASUREMENT.csv") %>%
  filter(GROUP == "FISH", 
         str_detect(ISO_YSAMP_N, "LML"), 
         ISO_FISH_N != "", 
         SAMPLE_TYPE == "TISSUE")



fish_lml = read.csv("../AFRP/Data/FISH_MEASUREMENT_LML.csv") %>%
  filter(FISH_N != "") %>%
  select(YSAMP_N:WEIGHT) %>% unique()

#Set up data ------------
data = left_join(measurement, sample) %>% 
  filter(is.na(D15N)==F & 
           is.na(D13C) == F & 
           WATER %in% c("LML")) 

no_date = read.csv("Data/JML.Data.Master.csv") %>%
  as.data.frame() %>%
  filter( Date == "",Species %in% c("SMB", "LT"),) %>%
  separate(FISH_N, into = c("species","water","YEAR","gear","dsamp")) %>%
  mutate(YEAR = as.Date(YEAR, format = "%m%d%y")) %>%
  separate(YEAR, into = c("YEAR", "MONTH", "DAY")) %>% 
  select(water, YEAR, MONTH,  species, C, N) %>%
  rename(TAXON = species, D13C = C, D15N = N, WATER = water) %>%
  mutate(GROUP = "FISH") %>% 
  select(WATER, YEAR, TAXON, GROUP, D13C, D15N) 

date = read.csv("Data/JML.Data.Master.csv") %>%  
  mutate(WATER = "LML") %>%
  filter(Species %in% c("SMB", "LT"), Date != "") %>%
  separate(Date, into = c("MONTH","DAY","YEAR")) %>% 
  select(WATER, YEAR, MONTH, Species, C, N) %>%
  rename(TAXON = Species, D13C = C, D15N = N) %>% 
  mutate(GROUP = "FISH") %>%
  select(WATER, YEAR, TAXON, GROUP, D13C, D15N)

JML_baselines = read.csv("Data/JML.Data.Master.csv") %>%
  filter(Species %in% c("Snails", "Periphyton", "Zoops","Heptageniidae")) %>%
  separate(Date, into = c("MONTH","DAY","YEAR"))%>%  
  mutate(WATER = "LML")  %>%
  rename(TAXON = Species, D13C = C, D15N = N, GROUP = Group)  %>%
  mutate(GROUP = case_when(TAXON == "Periphyton" ~ "algae", TAXON == "Snails" ~ "SNAIL", TAXON == "Zoops" ~ "ZOOP",
                           TAXON == "Heptageniidae" ~ "INSECT")) %>% 
  select(WATER, YEAR, TAXON, GROUP, D13C, D15N)

full = data %>% select(WATER, YEAR, TAXON, GROUP, D13C, D15N) %>%
  rbind(date) %>%
  rbind(no_date) %>%
  rbind(JML_baselines) %>%
  filter(TAXON != "" | GROUP == "ZOOP") 

#write.csv( full,"LML_iso_corrected.csv", row.names = FALSE)


LML_recent = measurement %>%
  left_join(fish_lml, by = c("ISO_FISH_N" = "FISH_N")) %>%
  separate(ISO_YSAMP_N, into = c("G", "WATER","YEAR","SAMP"))

## Early trends in d13C and d15N
length_trends = read.csv("Data/JML.Data.Master.csv") %>% filter(Length != "") %>%
  separate(Length, into = c("a","b","c","d","e","f","g","h","i")) %>%
  pivot_longer(a:i, names_to = "individual", values_to = "Length") %>%
  select(Species, Group, Date, C, N, Length) %>%
  na.omit() %>% 
  group_by(Species, Group, Date, C, N) %>%
  summarize(mean_length = mean(as.numeric(Length))) 

## Combined 
a = LML_recent %>% select(TAXON, GROUP, D13C, D15N, LENGTH) %>% mutate(period = "late") %>%
  rename(C = D13C, N = D15N, 
         mean_length = LENGTH,
         Species= TAXON, 
         Group = GROUP)


### Updated full 
full = length_trends %>% ungroup() %>% select(Species, Group, C, N, mean_length) %>%
  filter(Species %nin% c("CR","mussel")) %>%
  mutate(period = "early") %>% 
  rbind(a)

save(full, file= "Data/IsotopeComparisons/IsotopeData.RData")




