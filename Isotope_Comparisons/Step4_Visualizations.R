## Set up functions and load in critical data for visualizations
library(tidyverse)
library(wesanderson)
library(viridis)

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
full = data 

#### Load in baselines
load("Data/IsotopeComparisons/baselines.RData")

load("Data/IsotopeComparisons/baseline_simmr.RData")

#### Load in the overlap data
load("Data/IsotopeComparisons/overlap.df.RData")


#### Load in the area data
load("Data/IsotopeComparisons/med_area.RData")

#### Load in the legend
load("Data/IsotopeComparisons/simmr_legend.RData")

load("Data/Legend.col.RData") ## Has colors for plotting individuals species
## Species success

suc = data.frame(CODE = unique(legend$CODE), 
                 suc = c("ne", "na", "ne", 
                         "ne", "p", "p",
                         "p","p","ne","na",
                         "SMB1","SMB2","p",
                         "na", "na","na",
                         "ne" ))

## Differences in body size between the two sampling periods

## Note that RWF lengths are messed up and I have to see why theyre not getting found in the fish_measurement csv
full %>%
  ggplot(aes(x = TAXON, y = log(LENGTH), fill = period, col = period)) + 
  geom_boxplot() + 
  theme_minimal() + geom_jitter()


## Do the mean lengths during each period match the mean length of the isotope fish




## Trends in length X
full %>% 
  ggplot(aes(x = (LENGTH), y = D13C, col = TAXON)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme_minimal() +
  scale_x_log10() +
  facet_wrap(~period)


## Trends in length Y
full %>% 
  ggplot(aes(x = (LENGTH), y = D15N, col = TAXON)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme_minimal() +
  scale_x_log10() +
  facet_wrap(~period)




## View baselines Z

baselines %>% 
  ggplot(aes(x = (as.factor(community)), y = mean_d15N, col = GROUP , group = GROUP)) + 
  geom_point() + geom_line() +
  theme_minimal() +
  ylab("d15N") + 
  xlab("Community")



## Summary of fish in webs

summary = full %>% 
  select(-ISO_YSAMP_N) %>%
  unique() %>%
  group_by(period, TAXON) %>%
  summarize(count = n()) %>%
  pivot_wider(names_from = period, values_from = count) %>%
  mutate_all(~ ifelse(is.na(.), 0, .))

#write.csv(summary, file = "Data/IsotopeComparisons/Clean/summary.csv")


#### View webs



full_corrected = full %>%
  left_join(baseline_simmr %>% filter(GROUP == "B"), by = c("period" = "community")) %>%
  mutate(D13C_c = D13C - mean_d13C,
         D15N_c = (D15N - mean_d15N)/3.4 + 1)

full_corrected %>% 
  filter(TAXON %nin%c("NRD", "RWF","BND","BB","LLS", "ST")) %>%
  ggplot(aes(x = D13C_c, y = D15N_c, col = period)) + geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~TAXON, scales = "free")



full_corrected %>% 
  mutate(TAXON = case_when(TAXON == "SMB" & LENGTH < 100 ~ "SMB1",
                             TAXON == "SMB" & LENGTH >= 100 ~ "SMB2",
                             TAXON != "SMB" ~ TAXON)) %>%
  filter(is.na(TAXON) != T) %>%
  #filter(group %nin%c("NRD", "RWF","BND","BB","LLS", "ST")) %>%
  ggplot(aes(x = D13C_c, y = D15N_c, col = TAXON)) + geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~period) + 
  ylab("Trophic Position") +
  xlab("D13C (corrected)") + 
  theme_minimal() 

legend_sci = read.csv("Data/IsotopeComparisons/legend.csv")
colnames(legend_sci)

### Niche Area - Overview Graph
med_area %>%
  left_join(legend_sci, by = c("CODE" = "code")) %>%

  rename(TAXON = CODE) %>%
  filter(TAXON != "SMB") %>%
  ggplot(aes(x = (med_area), y = as.factor(scientific), fill = as.factor(community))) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) +
  theme_minimal(base_size = 12) + 
  xlab("SEAc 95% CI") +
  labs(y = NULL) +
  facet_wrap(~community, labeller = labeller(community = c("1" = "Post Initiation", "2"  = "Modern Observation"))) +
  scale_fill_manual("Community", values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  theme(legend.position = "none") +
  scale_x_log10()

### Niche area -- summary of credible differences between the two periods
med_area %>%
  ungroup %>%
  select(group, post_n, period, med_area) %>%
  unique() %>% 
  #filter(period == "early", post_n == 1)
  pivot_wider(names_from = period, values_from = med_area)%>%
  mutate(difference = late - early) %>%
  group_by(group) %>%
  summarize(mean = mean(difference, na.rm = T), 
            low =  quantile(difference,.025, na.rm = T),
            high = quantile(difference,.975, na.rm = T)) %>%
  left_join(legend) %>%
  left_join(legend_sci, by = c("CODE" = "code")) %>%
  filter(CODE != "SMB") %>% 
  select(scientific, mean, low, high) %>%
  unique() %>%
  mutate(sig = case_when(sign(low) == sign(high)~"*", sign(low) != sign(high) ~ NA)) %>%
  ggplot(aes(x = "Credible Difference", y = scientific, col = mean, shape = sig)) + 
  geom_point(size = 5) + 
  theme_minimal(base_size = 14) +
  scale_color_viridis_c() +
  theme(axis.title = element_blank()) +
  guides(shape = "none") + 
  labs(col = expression(Delta*"SEAc"))


## Calculating difference in SEAc between periods
med_area %>%
  #filter(!grepl("SMB", CODE)) %>%
  ## Filter for comparisons in both early and late  
  filter(CODE %in% c("PS", "MM","CS","CC", "SS", "RS","LT")) %>%
  group_by(community,period, post_n) %>%
  summarize(mean_area = mean(med_area)) %>%
  ungroup() %>% 
  select(-community) %>%
  pivot_wider(names_from = period, values_from = mean_area) %>%
  
  mutate(difference = late - early) %>%
  summarize(mean = mean(difference, na.rm = T), 
            low =  quantile(difference,.025, na.rm = T),
            high = quantile(difference,.975, na.rm = T)) 
  ggplot(aes(y = (mean_area), x = period, fill = period)) + 
 # stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) +
  #geom_boxplot()+
  theme_minimal(base_size = 14) + 
  xlab("Niche Area") +
  labs(y = NULL)
  
  
## Comparing niche charactreistics to species success
  


med_area %>%
  filter(CODE %nin% c("SMB", "SMB1", "SMB2")) %>%
  left_join(suc) %>%
  group_by(CODE, post_n, suc) %>%
  summarize(m_s = mean(med_area)) %>%
  ungroup() %>% 
  group_by(post_n, suc) %>%
  summarize(m = mean(m_s)) %>%
  ungroup() %>% pivot_wider(names_from = suc, values_from = m) %>%
  mutate(dif = p - ne) %>%
  summarize(mean = mean(dif) , lower = quantile(dif, 0.025), upper = quantile(dif, .975))


m_species = med_area %>%
  #filter(CODE %in% c("PS", "MM","CS","CC", "SS", "RS","LT")) %>%
  filter(CODE %nin% c("SMB", "SMB1", "SMB2")) %>%
  left_join(suc) %>%
  filter(suc != "na") %>%
  group_by( suc, CODE, community) %>%
  mutate(m = mean(med_area)) %>%
  unite("ID", c(suc, community)) 
m_quant = med_area %>%
  #filter(CODE %in% c("PS", "MM","CS","CC", "SS", "RS","LT")) %>%
  filter(CODE %nin% c("SMB", "SMB1", "SMB2")) %>%
  left_join(suc) %>%
  filter(suc != "na") %>%
  ungroup() %>%
  group_by( suc, post_n, community) %>%
  summarize(m = mean(med_area)) %>%
  unite("ID", c(suc, community)) 

ggplot(data = m_quant,aes( x = ID, y = m)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) +
  geom_point(data = m_species, aes(x= ID, y = m, col = CODE))+
  geom_line(data = m_species, aes(x = ID, y = m, col = CODE, group = CODE)) +
  theme_minimal() + 
  scale_x_discrete("Response|Period", labels = c("ne_1" = "(-) Post Initiation",
                                                 "ne_2" = "(-) Modern Observation",
                                                 "p_1" = "(+) Post Initiation",
                                                 "p_2" = "(+) Modern Observation")) + 
  ylab("SEAc") + 
  scale_color_manual("Species", values = (legend %>% 
                       select(CODE, color) %>%
                       unique() %>%
                       filter(CODE %in% c(m_species$CODE)))$color)



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


## Difference between periods only showing credible differences

overlap.df %>% 
  arrange(sorted_comparison, post) %>%
  separate(sorted_comparison, into = c("S1", "S2")) %>%
  left_join(legend_sci, by = c(S1 = "code")) %>%
  mutate(s1_sci = factor(scientific, levels = ov.fa)) %>%
  select(-scientific) %>%
  left_join(legend_sci, by = c("S2" = "code")) %>%
  mutate(s2_sci =  factor(scientific, levels = ov.fa)) %>%
  select(-scientific) %>%
  mutate(Values = round(as.numeric(Values), digits = 3))  %>%
  ungroup() %>%  
  unique() %>%
  pivot_wider(names_from = Community, values_from = Values) %>%
  group_by(S1, s1_sci, s2_sci, S2, post) %>%
  summarize(diff = `2`-`1`) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(S1, s1_sci, s2_sci, S2) %>%
  summarize(mean = mean(diff), 
            lower = quantile(diff, .025), 
            upper = quantile(diff, .95)) %>%
  ungroup() %>% 
  mutate(sig = sign(lower) == sign(upper)) %>%
  filter(sig == T) %>%
  filter(S1 != "SMB" & S2 != "SMB") %>%
  ggplot(aes(x = s1_sci, y = s2_sci, fill = mean*100 )) +
  geom_tile() +
  scale_fill_viridis_c()+
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5)) +
  labs(fill = expression(Delta * "Overlap (%)")) + 
  guides(size = "none")

## Plot of overlap in the two periods
early = overlap.df %>%
  filter(Community == 1) %>%
  separate(sorted_comparison, into = c("S1", "S2")) %>%
  mutate(s1 = S2, s2 = S1) %>%
  filter(s1 > s2) %>%
  select(s1, s2, Community, post, Values)
late = overlap.df %>%
  filter(Community == 2) %>%
  separate(sorted_comparison, into = c("S1", "S2")) %>%
  filter(S1 < S2) %>%
  rename(s1 = S1, s2 = S2)


ov.fa  = (legend_sci%>% 
  filter(code %nin% c("BND", "SMB","RT", "NRD")) %>%
  arrange(code))$scientific


unique.overlap = rbind(early, late) %>%
  filter(s1 != "SMB", s2 != "SMB") %>%
  group_by(Community, s1, s2) %>%
  summarize(value = mean(Values)) %>%
  left_join(legend_sci, by = c(s1 = "code")) %>%
  mutate(s1_sci = factor(scientific, levels = ov.fa)) %>%
  select(-scientific) %>%
  left_join(legend_sci, by = c("s2" = "code")) %>%
  mutate(s2_sci =  factor(scientific, levels = ov.fa)) %>%
  select(-scientific)

unique.overlap %>%
  ggplot(aes(x = reorder(s2_sci,as.factor(ov.or)), y = reorder(s1_sci, as.factor(ov.or)), fill = value)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Overlap (%)")

suc %>% colnames()
rbind(early, late) %>%
  filter(s1 != "SMB", s2 != "SMB") %>%
  #group_by(Community, s1, s2) %>%
 #summarize(value = mean(Values)) %>%
  left_join(legend_sci, by = c(s1 = "code")) %>%
  mutate(s1_sci = factor(scientific, levels = ov.fa)) %>%
  select(-scientific) %>%
  left_join(legend_sci, by = c("s2" = "code")) %>%
  mutate(s2_sci =  factor(scientific, levels = ov.fa)) %>%
  select(-scientific) %>%
  left_join(suc, by = c("s1" = "CODE")) %>%
  rename(suc1 = suc) %>%
  left_join(suc, by = c("s2" = "CODE")) %>%
  rename(suc2 = suc) %>%
  filter(s1 %in% c("SMB1", "SMB2") | s2 %in% c("SMB1", "SMB2")) %>%
  mutate(duplicate = suc1) %>%
  mutate(suc1 = case_when(duplicate %in% c("SMB1", "SMB2") ~ suc1, 
                          duplicate %nin% c("SMB1", "SMB2") ~ suc2)) %>%
  mutate(suc2 = case_when(duplicate %in% c("SMB1", "SMB2") ~ suc2, 
                          duplicate %nin% c("SMB1", "SMB2" )~ duplicate )) %>%
  group_by(post, suc1, suc2, Community) %>%
  summarise(mean = mean(Values)) %>%
  filter(suc2 != "na") %>%
  filter(!(suc1 == "SMB1" & suc2 == "SMB2"),
         !(suc1 == "SMB2" & suc2 == "SMB1")) %>%
  ggplot(aes(x = interaction(suc1, suc2), y = mean)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4)) + 
  facet_wrap(~Community, labeller = labeller("Community"= c("1" = "Post Initiation", 
                                                          "2" = "Modern Observtion"))) + 
  theme_minimal(base_size = 12) +
  #theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete("SMB Age|Species Response", labels = c("SMB1.ne" = "(J) | (-)", 
                                              "SMB2.ne" = "(A) | (-)",
                                              "SMB1.p" = "(J) | (+)",
                                              "SMB2.p" = "(A) | (+)")) + 
  ylab("Pairwise Overlap")
  





## incorporating success between periods




## Visualizing posterior mean ellipses

legend = legend.a %>% left_join(col_join)

legend = legend %>% arrange(CODE) %>% 
  filter(CODE != "NRD", 
         CODE != "SMB",
         CODE != "BND",
         !(CODE == "RWF" & community == 1)) 
for(h in 1:2){
  
  dat = data_setup(data_siber %>% filter(group != 12), h)
  
  
  spp=length(names(dat[[2]]))
  
  dat.corr = dat[[3]] %>%

    left_join(baseline_simmr %>%
                mutate(community = as.numeric(as.factor(community))) %>%
                filter(GROUP == "B")) %>%
    mutate(D13C_c = iso1 - mean_d13C,
           D15N_c = (iso2 - mean_d15N)/3.4 + 1)
  
  p = ggplot() + geom_point(dat.corr, mapping = aes(x = D13C_c, y = D15N_c, 
                                                    color = as.factor(group)), alpha = .25) + 
    theme_minimal(base_size = 14)
  
  ellip = array(0, dim=c((n.points),length(c(0,1)),spp))
  
  for(i in 1:spp){
    ellipse_data =  ellipse::ellipse(x = matrix(c(mean(dat[[2]][[i]][,1]),
                                                  mean(dat[[2]][[i]][,2]),
                                                  mean(dat[[2]][[i]][,3]),
                                                  mean(dat[[2]][[i]][,4])),
                                                2,2), 
                                     centre = c(mean(dat[[2]][[i]][,5]),
                                                mean(dat[[2]][[i]][,6])), 
                                     level = .4, 
                                     npoints = n.points)
    ellip[,,i] = ellipse_data
  }
  

    mean_d13C = (baseline_simmr %>%
             mutate(community = as.numeric(as.factor(community))) %>%
             filter(GROUP == "B",
                    community == h))$mean_d13C 
   
    mean_d15N = (baseline_simmr %>%
                         mutate(community = as.numeric(as.factor(community))) %>%
                         filter(GROUP == "B",
                                community == h))$mean_d15N

  
  x.vec = as.vector((ellip[,1,])) - mean_d13C
  y.vec = (as.vector((ellip[,2,]))  - mean_d15N)/3.4 + 1
  
  p = p +
    geom_path(aes(x = x.vec,
                  y = y.vec,
                  color = (rep(as.factor(unique((dat[[3]]$group))),each = 1000))),
              lwd = 1) +
    
    ylab("Trophic Position") + xlab("d13C (corrected)")   + 
    scale_color_manual(values = (legend %>% filter(community == h))$color, 
                       labels =(legend %>% filter(community == h))$CODE,
                       name = "Species")+ 
    theme(text = element_text(size = 13)) +
    ylim(1.75,4.25) + xlim(-11,5)
  
  
  print(p)
}
  
