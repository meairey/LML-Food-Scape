## Set up functions and load in critical data for visualizations
library(tidyverse)
library(viridis)
library(wesanderson)
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


legend_sci = read.csv("Data/IsotopeComparisons/legend.csv")




legend = legend %>% left_join(legend_sci, by = c("CODE" = "code"))

## Species success

suc = data.frame(CODE = sort(unique(legend$CODE)), 
                 suc = c("ne", ##bB
                         "na", ## BND - NA
                         "ne", ##CC
                         "ne", ## CS
                         "na", ## LLS
                         "p", ## LT
                         "p", ## MM
                         "na", ## NRD
                         "p", ## PS
                         "p", ## RS
                         "ne", ## RWF,
                         "na", ## SMB
                         "SMB1", ## SMB1
                         "SMB2", ## SMB2
                         "SMB3", ## SMB3
                         "p", ## SS
                          "ne", ## ST
                         "na" )) %>% ## WS
  filter(CODE != "SMB")

su_order <- c("A. nebulosus", "C. commersonii","L. cornutus", "P. cylindraceum", "S. atromaculatus", "S. fontinalis", 
              "C. cognatus", "L. gibbosus", "O. mordax", "S. namaycush", "S. salar","U. limi")



italic_labels = c(expression(italic(A.~nebulosus)),
                  #expression(italic(C.~commersonii)),
                  expression(italic(L.~cornutus)), 
                  expression(italic(P.~cylindraceum)),
                  expression(italic(S.~atromaculatus)),
                  expression(italic(S.~fontinalis)),
                  expression(italic(C.~cognatus)),
                  expression(italic(L.~gibbosus)),
                  expression(italic(O.~mordax)),
                  expression(italic(S.~namaycush)),
                  #expression(italic(S.~salar)),
                  expression(italic(U.~limi))
                  )



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


## Overlap between different size classes of SMB in LML between the historical and modern periods
full_corrected %>% 
  mutate(TAXON = case_when(TAXON == "SMB" & LENGTH < 100 ~ "SMB1",
                             
                           TAXON == "SMB" & LENGTH >= 100 & LENGTH < 200 ~ "SMB2",
                           TAXON == "SMB" & LENGTH >= 200 ~ "SMB3",
                             TAXON != "SMB" ~ TAXON)) %>%
  filter(TAXON %in%  c("SMB1", "SMB2", "SMB3")) %>%
  filter(is.na(TAXON) != T) %>%
  #filter(group %nin%c("NRD", "RWF","BND","BB","LLS", "ST")) %>%
  ggplot(aes(x = D13C_c, y = D15N_c, col = TAXON)) + geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~period, labeller = labeller(period = c("early" = "Historical", "late" = "Modern"))) + 
  ylab("Trophic Position") +
  xlab("D13C (corrected)") + 
  theme_minimal(base_size = 18) +
  scale_color_manual("Age Class",values = wes_palette("Darjeeling1", n = 3), labels = c("Juvenile", "Maturing", "Adult"))


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





  
  
## Comparing niche charactreistics to species success
  


med_area.1 = med_area %>% 
  filter(CODE %nin% c("SMB", "SMB1", "SMB2", "SMB3")) %>%
  left_join(suc) %>%
 filter(suc != "na")%>%
  unite("ID", c(suc, community), remove = F) %>%
  ungroup( ) %>%
  select(ID, post_n, med_area, CODE, color) %>%
  group_by(ID, CODE, color) %>%
  sample_n(1000) %>%
  mutate(post_index = c(1:1000))

m_species = med_area.1 %>%
  ungroup() %>%
  group_by( ID, CODE, color) %>%
  mutate(m = mean(med_area)) %>%

  left_join(legend %>%
              select(CODE, scientific, color) %>%
              unique(), by = c("CODE")) %>%
  #mutate(CODE = factor(CODE, levels = c("BB","CC","CS","RWF","ST","LT","MM","PS","RS","SS"))) %>%
  
  mutate(scientific = factor(scientific,levels = su_order))

m_quant = med_area.1 %>%
  ungroup() %>%
  group_by(ID, post_index) %>%
  summarize(m = mean(med_area)) 


## Final graph of niche area comparisons
ggplot(data = m_quant,aes( x = ID, y = m)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4), fill = "gray") +
  geom_point(data = m_species, aes(x= ID, y = m, col = scientific, shape = scientific), size = 5) +
  geom_line(data = m_species, aes(x = ID, y = m, col = scientific, group = scientific)) +
  theme_minimal(base_size = 18) + 
  scale_x_discrete("Response|Period",
                   labels = c("ne_1" = "Declined | Post",
                                                 "ne_2" = "Declined | Modern",
                                                 "p_1" = "Emerging | Post",
                                                 "p_2" = "Emerging | Modern")) + 
   
  scale_shape_manual(values = c(16,16,16,16,16,16,17,17,17,17,17,17), 
                     labels = italic_labels) +
  ylab("SEAc") + 
  scale_color_manual("Species", 
                     labels = italic_labels,
                     values = (legend %>% 
                                 mutate(scientific = factor(scientific, levels = su_order)) %>%
                                 arrange(scientific) %>%
                       select(scientific, color) %>%
                       unique() %>%
                       filter(scientific %in% c(m_species$scientific)))$color)  +  # or manually assign shapes to species for consistency
  labs(
    x = "Response | Period",
    y = "Niche Area (SEAc) 95% CI",
    color = "Species",
    shape = "Species"
  ) +

  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.major.y = element_line(color = "gray85"),
    legend.position = "right"
  )


## Bayes and ROPE test
library(BayesFactor)

# Calculate Probability of Direction (pd)

## comparing declining between periods
ne_1 = m_quant %>% filter(ID == "ne_1")
ne_2 = m_quant %>% filter(ID == "ne_2")
p_1 = m_quant %>% filter(ID == "p_1")
p_2 = m_quant %>% filter(ID == "p_2")
pd_value <- mean(ne_1 > ne_2)
pd_value = mean(p_1 > p_2)
pd_value = mean(ne_2 > p_1)
pd_value = mean(ne_2 > p_2)
pd_value


## SMB area before/after

m_smb.sp = med_area %>%

  filter(CODE %in% c("SMB", "SMB1", "SMB2", "SMB3")) %>%
  left_join(suc) %>%
  filter(suc != "na") %>%
  group_by( suc, CODE, community) %>%
  mutate(m = mean(med_area)) %>%
  unite("ID", c(suc, community), remove = F) %>%
  left_join(legend %>%
              select(CODE, scientific, color) %>%
              unique(), by = c("CODE")) %>%
  #mutate(CODE = factor(CODE, levels = c("BB","CC","CS","RWF","ST","LT","MM","PS","RS","SS"))) %>%
  
  mutate(scientific = factor(scientific,levels = c("A. nebulosus", "S. atromaculatus","L. cornutus",
                                                   "P. cylindraceum", "S. fontinalis",  "S. namaycush",
                                                   "U. limi" ,"L. gibbosus", "O. mordax","C. cognatus" 
                                                   )))

m_smb.qu = med_area %>%
  filter(CODE %in% c("SMB", "SMB1", "SMB2", "SMB3")) %>%
  left_join(suc) %>%
  filter(suc != "na") %>%
  ungroup() %>%
  group_by( suc, post_n, community) %>%
  summarize(m = mean(med_area)) %>%
  unite("ID", c(suc, community)) 


## Final graph of niche area comparisons
ggplot(data = m_smb.qu,aes( x = ID, y = m)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4), fill = "gray") +
  geom_point(data = m_smb.sp, aes(x= ID, y = m, col = scientific, shape = scientific), size = 5)+
  geom_line(data = m_smb.sp, aes(x = ID, y = m, col = scientific, group = scientific)) +
  theme_minimal(base_size = 18) + 
  scale_x_discrete("Response|Period", labels = c("ne_1" = "Declined | Post",
                                                 "ne_2" = "Declined | Modern",
                                                 "p_1" = "Recovered | Post",
                                                 "p_2" = "Recovered | Modern")) + 
   
  scale_shape_manual(values = c(16,16,16,16,16,17,17,17,17,17)) +
  ylab("SEAc") + 
  scale_color_manual("Species", values = (legend %>% 
                       select(scientific, color) %>%
                       unique() %>%
                       filter(scientific %in% c(m_species$scientific)))$color)  +  # or manually assign shapes to species for consistency
  labs(
    x = "Response | Period",
    y = "Niche Area (SEAc) 95% CI",
    color = "Species",
    shape = "Species"
  ) +

  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.major.y = element_line(color = "gray85"),
    legend.position = "right"
  )



### Overlap figures


exclude = legend$CODE %>% unique() 
exclude = exclude[which(exclude %in% c("PS", "MM","CS","CC","SS","RS","LT", "SMB") ==F)]

overlap.df %>%
  #filter(!grepl("SMB", sorted_comparison)) %>%
  
  mutate(period = case_when(Community == 1 ~ "early", Community == 2 ~ "late")) %>%
  #filter(!grepl(exclude, sorted_comparison)) %>%
  group_by(Community, post, period) %>%
  summarize(mean_value = mean(Values, na.rm = T)) %>%
  ggplot(aes(x = period, y = mean_value, fill = period)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 14) + 
  ylab("Average Pairwise Overlap") +
  scale_fill_manual("Perid", values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  scale_x_discrete(labels = c("Early", "Late")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())



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


ov.fa  = (legend_sci%>% 
  filter(code %nin% c("BND", "SMB","RT", "NRD")) %>%
  arrange(code))$scientific

suc_ov = rbind(early, late) %>%
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
  filter(s1 %in% c("SMB1", "SMB2","SMB3") | s2 %in% c("SMB1", "SMB2","SMB3")) %>%
  mutate(duplicate = suc1) %>%
  mutate(suc1 = case_when(duplicate %in% c("SMB1", "SMB2", "SMB3") ~ suc1, 
                          duplicate %nin% c("SMB1", "SMB2", "SMB3") ~ suc2)) %>%
  mutate(suc2 = case_when(duplicate %in% c("SMB1", "SMB2","SMB3") ~ suc2, 
                          duplicate %nin% c("SMB1", "SMB2" ,"SMB3") ~ duplicate )) %>%
  filter(!(suc1 == "SMB1" & suc2 == "SMB2"),
         !(suc1 == "SMB2" & suc2 == "SMB1"),
         !(suc1 == "SMB1" & suc2 == "SMB3"), 
         !(suc1 == "SMB2" & suc2 == "SMB3"),
         !(suc1 == "SMB3" & suc2 == "SMB1"),
         !(suc1 == "SMB3" & suc2 == "SMB2"),
         !(suc1 == "SMB3" & suc2 == "SMB3")) %>%
  na.omit()



suc_ov.quant = suc_ov %>%
  group_by(post, suc1, suc2, Community) %>%
  summarise(mean = mean(Values)) %>%
  filter(suc2 != "na") %>%
  filter(!(suc1 == "SMB1" & suc2 == "SMB2"),
         !(suc1 == "SMB2" & suc2 == "SMB1"),
         !(suc1 == "SMB1" & suc2 == "SMB3"), 
         !(suc1 == "SMB2" & suc2 == "SMB3"),
         !(suc1 == "SMB3" & suc2 == "SMB2"),
         !(suc1 == "SMB3" & suc2 == "SMB3")) 

suc_ov.species =  suc_ov %>%
  group_by(s1, s2, suc1, suc2, Community) %>%
  summarise(mean = mean(Values)) %>%
  filter(suc2 != "na") %>%
  filter(!(suc1 == "SMB1" & suc2 == "SMB2"),
         !(suc1 == "SMB2" & suc2 == "SMB1"),
         !(suc1 == "SMB1" & suc2 == "SMB3"), 
         !(suc1 == "SMB2" & suc2 == "SMB3"),
         !(suc1 == "SMB3" & suc2 == "SMB1"),
         !(suc1 == "SMB3" & suc2 == "SMB2"),
         !(suc1 == "SMB3" & suc2 == "SMB3")) %>%
  mutate(s1_dup = s1) %>%
  mutate(s1 = case_when(s1 %in% c("SMB1", "SMB2","SMB3") ~ s2,
                        s1 %nin% c("SMB1", "SMB2","SMB3") ~ s1 )) %>% 
  mutate(s2 = case_when(s2 %in% c("SMB1", "SMB2","SMB3") ~ s2, 
                        s2 %nin% c("SMB1", "SMB2","SMB3") ~ s1_dup)) %>%
  select(-s1_dup) %>%
  left_join(legend %>%
              mutate(community = as.character(community)), by = c("s1" = "CODE", "Community" = "community")) %>%
  mutate(scientific = factor(scientific, levels = su_order))


su.order = (suc_ov.species %>% 
  ungroup() %>%
  select(suc2, scientific) %>%
  unique() %>%
  arrange(suc2, scientific))$scientific

## Trying a different facet

ggplot(data = suc_ov.quant, aes(x = as.factor(Community), y = mean)) + 
  stat_summary(fun.data=quantiles_95, geom="boxplot", aes(width=0.4), alpha = .5, fill = "gray") + 
  theme_minimal(base_size = 14) +
  #theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete("Community", labels = c("1" = "Post Initiation", 
                                              "2" = "Modern Observation")) + 
  geom_line(data = suc_ov.species, aes(x = Community, y = mean, col = scientific, group = scientific)) +
  ylab("Pairwise Overlap") + 
  geom_point(data = suc_ov.species, aes(x = Community, y = mean, col = scientific , shape = scientific), size = 4) +
  facet_wrap(~interaction(suc1, suc2), 
             labeller = labeller("interaction(suc1, suc2)"=  c(
                                              "SMB1.ne" = "(J) | Declined", 
                                              "SMB2.ne" = "(M) | Declined",
                                              "SMB1.p" = "(J) | Recovered",
                                              "SMB2.p" = "(M) | Recovered",
                                              "SMB3.ne" = "(A) | Declined",
                                              "SMB3.p" = "(A) | Recovered")),
             ncol = 2, dir = "v") +
  #labs(col = "Species") + 
  scale_x_discrete("Period", labels = c("1" = "Post \n Initiation",
                                                 "2" = "Modern \n Observation" )) + 
   
  scale_shape_manual(values =c(16,16,16,16,16,17,17,17,17,17), # c(16, 17, 16, 17, 17, 16, 16, 16, 17, 17)
                     #labels = italic_labels
                     ) +
  ylab("SEAc") + 
  scale_color_manual("Species", 
                     #labels = italic_labels,
                     values = (legend %>%
                                 mutate(scientific = factor(scientific, levels = su.order)) %>%
                                 arrange(scientific) %>%
                       select(scientific, color) %>%
                       unique() %>%
                       filter(scientific %in% c(m_species$scientific)))$color)  +  # or manually assign shapes to species for consistency
  labs(
    x = "Response | Period",
    y = expression("Overlap w/"~ italic(M.~dolomieu)~"(%) 95% CI"),
    color = "Species",
    shape = "Species"
  ) +

  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.major.y = element_line(color = "gray85"),
    legend.position = "right"
  )




(legend %>%
                                 arrange(CODE) %>%
                       select(scientific, color) %>%
                       unique() %>%
                       filter(scientific %in% c(m_species$scientific)))


## incorporating success between periods
cred.ints = suc_ov.quant %>%
  unite("ID", suc1, suc2, Community) %>%
  ungroup() %>%
  rename(value = mean)

print(cred.ints)
  
  
  
## juvenile overlaps neg
smb1.ne.1 = suc_ov.quant %>% filter(suc1 == "SMB1", suc2 == "ne", Community == 1)
smb1.ne.2 = suc_ov.quant %>% filter(suc1 == "SMB1", suc2 == "ne", Community == 2)
pd.smb1.ne.comp = mean(smb1.ne.1$mean > smb1.ne.2$mean) 

## juvenile overlaps pos 
smb1.p.1 = suc_ov.quant %>% filter(suc1 == "SMB1", suc2 == "p", Community == 1)
smb1.p.2 = suc_ov.quant %>% filter(suc1 == "SMB1", suc2 == "p", Community == 2)
pd.smb1.p.comp = mean(smb1.p.1$mean > smb1.p.2$mean)

## juvenile smb, overlap comparisons between sensitive and recovered

pd.smb1.comp.pre = mean(smb1.ne.1$mean > smb1.p.1$mean) ## neg/pos between early
pd.smb.1.comp.post = mean(smb1.ne.2$mean > smb1.p.2$mean) ## neg/pos between post
pd.smb1.comp.nep.12 = mean(smb1.ne.2$mean > smb1.p.1$mean)
#pd.smb1.comp.prepost = mean(smb1.ne.1$mean > smb1.p.2$mean) 
pd.smb1.comp.prepost2 = mean(smb1.ne.2$mean > smb1.p.2$mean)



## now comparing with medium smb

smb2.ne.2 = suc_ov.quant%>% filter(suc1 == "SMB2", suc2 == "ne", Community == 2)
smb2.p.2 = suc_ov.quant %>% filter(suc1 == "SMB2", suc2 == "p", Community == 2)
pd.smb2 = mean(smb2.ne.2$mean > smb2.p.2$mean)

## Are smb2 similar to smb1 declined

pd.smb2.smb1.post.ne = mean(smb2.ne.2$mean > smb1.ne.2$mean)
pd.smb2.smb1.prepost.ne = mean(smb2.ne.2$mean > smb1.ne.1$mean)


## Are smb2 similar to smb1 recovered

pd.smb2.smb1.post.p = mean(smb2.ne.2$mean > smb1.p.2$mean)
pd.smb2.smb1.prepost.p = mean(smb2.ne.2$mean > smb1.p.1$mean)



## are smb3 similar?

smb3.ne.1 = suc_ov.quant %>% filter(suc1 == "SMB3", suc2 == "ne", Community == 1)
smb3.ne.2 = suc_ov.quant %>% filter(suc1 == "SMB3", suc2 == "ne", Community == 2)


mean(smb3.ne.1$mean > smb3.ne.2$mean) ## both declined smb3 are the same
mean(smb3.ne.1$mean > smb1.ne.2$mean) ## smb3 1 

smb3.ne.2 = suc_ov.quant %>% filter(suc1 == "SMB3", suc2 == "p", Community == 1)
smb3.ne.2 = suc_ov.quant %>% filter(suc1 == "SMB3", suc2 == "p", Community == 2)






## for loop
library(tidyverse)

# Get all unique IDs
ids <- unique(cred.ints$ID)

# Initialize an empty tibble to store results
pd_results <- tibble(ID1 = character(), ID2 = character(), PD = numeric())

# Loop over all unique ID combinations
combs <- combn(ids, 2, simplify = FALSE)

for (pair in combs) {
  id1 <- pair[1]
  id2 <- pair[2]
  
  # Extract posterior samples
  vec1 <- cred.ints %>% filter(ID == id1) %>% pull(value)
  vec2 <- cred.ints %>% filter(ID == id2) %>% pull(value)
  
  # Calculate PD (Probability that id1 > id2)
  pd_val <- mean(vec1 > vec2)
  
  # Store both directions
  pd_results <- pd_results %>%
    add_row(ID1 = id1, ID2 = id2, PD = pd_val) %>%
    add_row(ID1 = id2, ID2 = id1, PD = 1 - pd_val)
}

pd_results = pd_results %>%
  filter(ID1 > ID2,
         ID1 != ID2) 

pd_results %>%
  arrange(ID1) %>%
   filter(PD < .05 | PD > .95) %>%
  print(n=100) 










## Visualizing posterior mean ellipses

legend = legend.a %>% left_join(col_join)

legend.ellips = legend %>% arrange(CODE) %>% 
  filter(CODE != "NRD", 
         CODE != "SMB",
         CODE != "BND",
         !(CODE == "RWF" & community == 1)) 

plots = list()
for(h in 2){
  
  dat = data_setup(data_siber %>% filter(group != 12), h)
  
  
  spp=length(names(dat[[2]]))
  
 
  dat.corr = dat[[3]] %>%

    left_join(baseline_simmr %>%
                mutate(community = as.numeric(as.factor(community))) %>%
                filter(GROUP == "B")) %>%
    mutate(D13C_c = iso1 - mean_d13C,
           D15N_c = (iso2 - mean_d15N)/3.4 + 1) %>%
    left_join(legend)
  
  p2 = ggplot() + geom_point(dat.corr, mapping = aes(x = D13C_c, y = D15N_c, 
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
  
  p2= p2 +
    geom_path(aes(x = x.vec,
                  y = y.vec,
                  color = (rep(as.factor(unique((dat[[3]]$group))),each = 1000))),
              lwd = 1) +
    
    ylab("Trophic Position") + xlab("d13C (corrected)")   + 
    scale_color_manual(values = (legend.ellips %>%
                                   filter(CODE %in% dat.corr$CODE) %>%
                                   filter(community == h))$color, 
                       labels =(legend.ellips %>%
                                   filter(CODE %in% dat.corr$CODE) %>%
                                  filter(community == h))$scientific,
                       name = "Species")+ 
    theme(text = element_text(size = 13)) +
    ylim(1.75,4.25) + xlim(-11,5)
  
  

  
}


### Trophic Position ----------------------------




## Plot
dat_with_residuals %>% 
  left_join(suc, by = c("TAXON" = "CODE")) %>%
  mutate(
    sampled_both_periods = ifelse(
      n_distinct(period) == 2,  # 2 different periods exist
      "yes", 
      "no"
    )
  ) %>%
  left_join(legend_sci, by = c("TAXON" = "code")) %>%
  mutate(sampled_both_periods = case_when(sampled_both_periods == "yes" & TAXON %in% c("BND", "NRD", "LLS", "ST", "BB", "WS", "RWF", "SS") ~ "no",
                                          sampled_both_periods == "yes" & TAXON %nin% c("BND", "NRD", "LLS", "ST", "BB","WS","RWF", "SS") ~ "yes",
                                          sampled_both_periods == "no" ~ "no")) %>%
  ungroup()%>%
  ggplot(aes(x = scientific, y = D15N_c, color = period)) +
  geom_boxplot(position = position_dodge(width = 0.8), aes(fill = interaction(sampled_both_periods, period)), alpha = .25) +
  # Add asterisks
  geom_text(
    data = sig_results %>% filter(significance != ""),   # Only plot significant results
    aes(x = scientific, y = label_y, label = significance),
    inherit.aes = FALSE,
    size = 8,
    vjust = 0
  ) +
  facet_wrap(~suc, scales = "free_x", 
             labeller = labeller(suc = c("na" = "No Change",
                                         "ne" = "Declined", 
                                         "p" = "Recovered"))) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Trophic Position", color = "Period") +
  geom_jitter(width = .05, height = 0) + 
  scale_fill_manual(values = c("transparent", wes_palette("AsteroidCity1", n = 2)[1],
                               "transparent", wes_palette("AsteroidCity1", n = 2)[2])) +
  scale_color_manual("Period", values = c( wes_palette("AsteroidCity1", n = 2)[1],  wes_palette("AsteroidCity1", n = 2)[2])) + 
  guides(fill = "none") +
  theme(axis.text.x =element_text(face = "italic"), 
        axis.title.x = element_blank()) 

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

## Differences in body size between the two sampling periods

## Note that RWF lengths are messed up and I have to see why theyre not getting found in the fish_measurement csv
full %>%
  ggplot(aes(x = TAXON, y = log(LENGTH), fill = period,
             col = period)) + 
  geom_boxplot() + 
  theme_minimal() + geom_jitter()

full %>%
  filter(TAXON == "RWF")


