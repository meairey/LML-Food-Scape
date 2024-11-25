## Set up functions and load in critical data for visualizations

#### Quantile function for plotting credible intervals
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
#### Inverse of %in% to exclude things from strings
`%nin%`= Negate(`%in%`)

#### Load in all isotope data then convert strings to numerics
load("Data/IsotopeComparisons/IsotopeData.RData")
full = full %>% 
  mutate(C = as.numeric(C), 
         N = as.numeric(N))

#### Load in baselines
load("Data/IsotopeComparisons/baselines.RData")

load("Data/IsotopeComparisons/baseline_simmr.RData")

#### Load in the overlap data
load("Data/IsotopeComparisons/overlap.df.RData")


#### Load in the area data
load("Data/IsotopeComparisons/med_area.RData")

#### Load in the legend
load("Data/IsotopeComparisons/simmr_legend.RData")


## Differences in body size between the two sampling periods

## Note that RWF lengths are messed up and I have to see why theyre not getting found in the fish_measurement csv
full %>%
  ggplot(aes(x = Species, y = log(mean_length), fill = period)) + 
  geom_boxplot() + 
  theme_minimal()


## Do the mean lengths during each period match the mean length of the isotope fish




## Trends in length X
full %>% 
  ggplot(aes(y = (mean_length), x = C, col = Species)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(~period)




## View baselines 

baselines %>% 
  ggplot(aes(x = (as.factor(community)), y = mean_d15N, col = GROUP , group = GROUP)) + 
  geom_point() + geom_line() +
  theme_minimal() +
  ylab("d15N") + 
  xlab("Community")


#### View webs



full_corrected = full %>%
  left_join(baseline_simmr %>% filter(GROUP == "B"), by = c("period" = "community")) %>%
  mutate(D13C_c = C - mean_d13C,
         D15N_c = (N - mean_d15N)/3.4 + 1)

full_corrected %>% 
  filter(Species %nin%c("NRD", "RWF","BND","BB","LLS", "ST")) %>%
  ggplot(aes(x = D13C_c, y = D15N_c, col = period)) + geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~Species, scales = "free")



full_corrected %>% 
  mutate(Species = case_when(Species == "SMB" & mean_length < 200 ~ "SMB1",
                             Species == "SMB" & mean_length >= 200 ~ "SMB2",
                             Species != "SMB" ~ Species)) %>%
  #filter(group %nin%c("NRD", "RWF","BND","BB","LLS", "ST")) %>%
  ggplot(aes(x = D13C_c, y = D15N_c, col = Species)) + geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~period) + 
  ylab("Trophic Position") +
  xlab("D13C (corrected)") + 
  theme_minimal() 


### Niche Area 
med_area %>%
  rename(Species = CODE) %>%
  ggplot(aes(x = (med_area), y = as.factor(Species), fill = as.factor(community))) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) +
  theme_minimal(base_size = 12) + 
  xlab("Niche Area 95% CI") +
  labs(y = NULL) +
  facet_wrap(~community, labeller = labeller(community = c("1" = "Early", "2"  = "Late"))) +
  scale_fill_manual("Community", values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  theme(legend.position = "none")




library(wesanderson)


med_area %>%
  #filter(!grepl("SMB", CODE)) %>%
  ## Filter for comparisons in both early and late
  filter(CODE %in% c("PS", "MM","CS","CC", "SS", "RS","LT")) %>%
  group_by(community,period, post_n) %>%
  summarize(mean_area = mean(med_area)) %>%
  
  ggplot(aes(y = (mean_area), x = period, fill = period)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) +
  #geom_boxplot()+
  theme_minimal(base_size = 14) + 
  xlab("Niche Area") +
  labs(y = NULL)


## Overlap Boxplots




exclude = legend$CODE %>% unique() 
exclude = exclude[which(exclude %in% c("PS", "MM","CS","CC","SS","RS","LT", "SMB") ==F)]

overlap.df %>%
  #filter(!grepl("SMB", sorted_comparison)) %>%
  mutate(period = case_when(Community == 1 ~ "early", Community == 2 ~ "late")) %>%
  filter(!grepl(exclude, sorted_comparison)) %>%
  group_by(Community, post, period) %>%
  summarize(mean_value = mean(Values)) %>%
  ggplot(aes(x = period, y = mean_value, fill = period)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 14) + 
  ylab("Average Pairwise Overlap") +
  scale_fill_manual("Perid", values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  scale_x_discrete(labels = c("Early", "Late")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())
