fish_sample = read.csv("../AFRP/Data/FISH_SAMPLE_2024.csv") %>%
  select(YSAMP_N, WATER, SITE_N, MONTH, YEAR)

sic_sample = read.csv("../AFRP/Data/isotope_sample.csv") %>%
  select(ISO_YSAMP_N, WATER, SITE_N, MONTH, YEAR) %>%
  rename(YSAMP_N = ISO_YSAMP_N)


fish_measurement = read.csv("../AFRP/Data/FISH_MEASUREMENT_LML.csv")


sample = rbind(fish_sample, sic_sample)  

  
  

si_measurement =  read.csv("Data/SI_MEASUREMENT.csv") %>%
  filter(str_detect(ISO_YSAMP_N, "LML")) %>%
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  separate(ISO_FISH_N, into = c("code", "water", "date", "gear", "num")) %>%
  filter(SAMPLE_TYPE == "TISSUE",
         GROUP == "FISH")

si_measurement %>%
  ggplot(aes(x = D13C, y = D15N, col = YEAR)) + 
  geom_point() + 
  stat_ellipse()  + 
  facet_wrap(~code, scales = "free")


si_measurement %>% left_join(fish_measurement, by = c("ITEM_N" = "FISH_N"))







data.late = read.csv("Data/SI_MEASUREMENT.csv") %>%
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
  mutate(community = 1, 
         iso_1 = scale(as.numeric(iso_1)),
         iso_2 = scale(as.numeric(iso_2))) %>%
  na.omit() %>%
  select(community, group, iso_1, iso_2)

library(wesanderson)
data.late %>%
  ggplot(aes(x = iso_1, y = iso_2, col = group)) + 
  theme_minimal() +
  geom_point() + 
  stat_ellipse() +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 12, type = "continuous"))
