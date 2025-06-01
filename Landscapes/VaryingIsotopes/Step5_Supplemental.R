## How many individuals caught in each period
`%nin%` = Negate(`%in%`)

LML.data %>%
  separate(SPECIES, into = c("SPECIES", "AGE")) %>%
  group_by(SPECIES, YEAR) %>%
  summarize(sum_species = n()) %>%
  mutate(period = case_when(YEAR == 2000 ~ "pre", 
                            YEAR %in% c(2001:2005) ~ "post",
                            YEAR %in% c(2019:2023) ~ "modern")) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(SPECIES, period) %>%
  summarize(total_years = n(), 
            total_catch = sum(sum_species)) %>%
  print(n = 100)


fish = read.csv("../2. clean_AFRPData/FISH_MEASUREMENT_LML.csv")
sample = read.csv("../2. clean_AFRPData/FISH_SAMPLE_2024.csv")

data = left_join(fish, sample, by = "YSAMP_N") %>%
  filter(MONTH < 7, 
         GEAR == "BEF", 
         GEAR_CODE == "NAF",
         YEAR != 2002)


data  %>%
  
    filter(SPECIES %nin% c("FI", "RSN","SPL","UID","UIF","cs","NF", "TP", "LWF", "MAM")) %>%
  group_by(SPECIES, YEAR) %>%
  summarize(sum_species = n()) %>%
  mutate(period = case_when(YEAR == 2000 ~ "pre", 
                            YEAR %in% c(2001:2004) ~ "post",
                            YEAR %in% c(2019:2023) ~ "modern"),
         num_year = case_when(YEAR == 2000 ~ 1, 
                            YEAR %in% c(2001:2004) ~ 3,
                            YEAR %in% c(2019:2023) ~ 3)) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(SPECIES, period, num_year) %>%
  summarize(total_years = n(), 
            total_catch = sum(sum_species)) %>%
  mutate(proportion_years = total_years / num_year,
         total_catch = total_catch/num_year) %>%
  select(proportion_years, everything()) %>%
 # filter(SPECIES == "BB")
  #select(SPECIES, proportion_years) %>% unique() %>% print(n = 100)
  ungroup() %>%
  select(-total_years, -proportion_years, -num_year) %>%
  pivot_wider(names_from = period, values_from = total_catch, values_fill = 0) %>%
  select(SPECIES, pre, post, modern)
  print(n = 100)

write.csv(test, "summary_catch.csv")
