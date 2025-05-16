library(tidyverse)
library(wesanderson)

## Successful in LML SUC## Species success

#### Load in the legend
load("Data/IsotopeComparisons/simmr_legend.RData")

load("Data/Legend.col.RData") ## Has colors for plotting individuals species


legend_sci = read.csv("Data/IsotopeComparisons/legend.csv")




legend = legend %>% left_join(legend_sci, by = c("CODE" = "code"))



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

## Legend
## Legend
species.early = c("BB", "CC","CS","MM","PS","SMB","SS","WS") # 8 long

species.pre = c("BB", "CC","CS","PS","SMB","SS","WS") #7 long

species.late =  c( "CC","CS","LT","MM","PS","RS","SMB","SS","WS") # 9 long
## Legend for joining together the posteriors - includes group/community combinations
legend = data.frame(group = as.character(c(1:9)), 
                    A = c(species.pre, rep("NA", 2)),
                    B = c(species.early, rep("NA",1)),
                    C  = species.late
                     
                    )%>%
  pivot_longer(2:4, names_to = "community", values_to = "species") %>%
  arrange(community) %>%
  mutate(community = as.numeric(as.factor(community)))
  
legend

bayes_cri <- function(x) {
  data.frame(
    y = mean(x),  # Mean of posterior samples
    ymin = quantile(x, 0.025),  # 2.5% quantile (Lower bound)
    ymax = quantile(x, 0.975)   # 97.5% quantile (Upper bound)
  )
}



mean.legend = legend$species %>% unique() %>% as.data.frame() %>% 
  na.omit() %>% 
  filter(. != "NA") %>%
  mutate(group = as.character(as.numeric(as.factor(.)))) %>%
  arrange(group) %>%
  rename("species" = ".")
### Volume Figures
load("Data/VaryingIsotopesData/Final/total.vol.frame (1).RData")



total_vols = total.vol %>%
  group_by(group, post, community) %>%
  summarise(total_vol = sum(total_vol)) %>%
  left_join(mean.legend)


## Total volume of each species individually
species_vols= total_vols %>%
  as.data.frame() %>%
  complete(species, post, community) %>%
  replace_na(list(total_vol = 0)) 

## Contemplating the incorporation of sensitive vs. tolerant natives in this conversation
summary.vol = species_vols %>% 
  filter(species == "SMB") %>%
  left_join(suc, by = c("species" = "CODE")) %>%
  group_by(suc, post, community) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ungroup() %>%
  group_by(suc, community) %>%
  summarize(mean = quantiles_95(total_vol)[3], min = quantiles_95(total_vol)[2],
            max = quantiles_95(total_vol)[4])

species_spec.vol = species_vols %>% 
  filter(species %nin% c("SMB")) %>%
  left_join(suc, by = c("species" = "CODE")) %>%
  group_by(suc, species, community) %>%
  summarize(mean =  quantiles_95(total_vol)[3], min = quantiles_95(total_vol)[2],
            max = quantiles_95(total_vol)[4])


  #### Load in the legend
load("Data/IsotopeComparisons/simmr_legend.RData")

load("Data/Legend.col.RData") ## Has colors for plotting individuals species


legend_sci = read.csv("Data/IsotopeComparisons/legend.csv")



legend = legend %>% left_join(legend_sci, by = c("CODE" = "code"))


ggplot() +
  geom_bar(data = species_spec.vol, aes(x = suc, y = mean, fill = species), stat = "identity") +
  facet_wrap(~community, labeller =labeller(community = c("1" = "Pre-Removal", "2" = "Post-Initiation", "3" = "Modern Observation"))) +
  theme_minimal(base_size = 14) +
  scale_x_discrete("Response", labels = c("ne" = "Declined", 
                              "p" = "Recovered",
                              "na" = "None")) +
  scale_fill_manual("Species", values = (legend %>% filter(CODE %in% c("BB","CC","CS","LT","MM","PS","RS", "SS", "WS")))$color,
                    labels =(legend %>% filter(CODE %in% c("BB","CC","CS","LT","MM","PS","RS","SS", "WS")))$scientific ) +
 # geom_hline(data = summary.vol, aes( yintercept = mean), lty = 2) +
  ylab("Niche Volume") +
  theme(legend.position = "bottom")
  

### Rugosity Figures

load("Data/VaryingIsotopesData/Final/rugosity.summary.RData")

rugosity.summary$post %>% unique() %>% length()

rugosity.summary  %>%
  ggplot(aes(x = as.factor(community), y = rug, col = as.factor(community))) + 
  #geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("Rugosity") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  scale_color_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 3),labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  ) 


r1 = rugosity.summary %>% filter(community == 1)
r2 = rugosity.summary %>% filter(community == 2)
r3 = rugosity.summary %>% filter(community == 3)

pd1 = mean(r1$rug > r2$rug)
pd2 = mean(r1$rug > r3$rug)

### DNHP
nhsp$post %>% unique() %>% length()

load("Data/VaryingIsotopesData/Final/nhsp.RData")

nhsp %>%
  ggplot(aes(x = as.factor(community), y = values, col = as.factor(community))) + 
  #geom_boxplot() +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("DNHP") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  facet_wrap(~metric, scales = "free",
             labeller = labeller("metric" = c("mean" = "Mean DNHP", "sd" = "SD DNHP"))) +

 stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )   + 
  scale_color_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 3),
                     labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) 

## Any sig differences between 1st and ... for mean
nhsp.1 = nhsp %>% filter(metric == "mean", 
                community == 1)
nhsp.2 = nhsp %>% filter(metric == "mean", 
                         community == 2)
nhsp.3 = nhsp %>% filter(metric == "mean", 
                         community == 3)

mean(nhsp.1$values > nhsp.2$values)
mean(nhsp.1$values > nhsp.3$values)
## No. None of those comparisons are significant


## Any sig differences between 1st and ... for sd

nhsp.1 = nhsp %>% filter(metric == "sd", 
                community == 1)
nhsp.2 = nhsp %>% filter(metric == "sd", 
                         community == 2)
nhsp.3 = nhsp %>% filter(metric == "sd", 
                         community == 3)

mean(nhsp.1$values > nhsp.2$values)
mean(nhsp.1$values > nhsp.3$values)
## No. None of those comparisons are significant

