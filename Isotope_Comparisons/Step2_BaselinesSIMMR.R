## Baselines and SIMMR Script
# Libraries
library(rjags)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
`%nin%` = Negate(`%in%`)


## Script


baselines.all = baselines %>% 
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

load( file = "Data/IsotopeComparisons/baselines.RData")

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

load( file = "Data/IsotopeComparisons/baseline_simmr.RData")

#baseline_simmr[4,4] = 1.1 ## Artificially lowering the d15N of zooplankton to see how it influences the diet contributions


load(file = "Data/IsotopeComparisons/IsotopeData.RData")

data_simmr = data %>%
  select(D13C, D15N, period, TAXON)  %>%
  rename(iso1 = D13C, 
         iso2 = D15N, 
         community = period, 
         group = TAXON) 

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


#save(CI_list, file = "Data/IsotopeComparisons/SimmrOutput.RData")


load(file = "Data/IsotopeComparisons/SimmrOutput.RData")

load("Data/VaryingIsotopesData/filtered_ellipses.RData")
load(file = "Data/VaryingIsotopesData/MeanEllipses.RData")




## I want to go through and modify each simmer species output for each period by their estimated volume

load(file = "Data/VaryingIsotopesData/species_vols.RData")

## Need a simmr specific legend

mod.leg = data.frame(early = c("BB","BND","CC","CS","LT","MM","PS","RS","RWF","SMB","SS","WS","NA"),
           late = c("CC","CS","LLS","LT","MM","NRD","PS","RS","RWF","SMB","SS","ST","WS")) %>%
  mutate(index = c(1:13)) %>%
  pivot_longer(c(early, late), names_to = "community", values_to = "CODE") %>%
  arrange(community) %>% 
  mutate(simmr_join = case_when(community == "early" ~ 1, community == "late" ~ 2)) %>%
  group_by(community) %>% 
  mutate(group = (c(1:13))) %>%
  ungroup() %>%
  select(-community)

rownames3 = data.frame(group = rep(c(1:13), each = 25), 
                      quantile = rep(rep(c("2.5", "25","50","75","97.5"), each = 5), by = 13),
                      type = rep(rep(c("dev", "ben","pel","sdiso1","sdiso2"), rep = 5),by=13),
                      community =3, 
                      simmr_join = 2)

simmr_estimates3 = data.frame(value = unlist(CI_list[[2]]$quantiles)) %>%
  as.matrix() %>%
  cbind(rownames3) %>%
  left_join(mod.leg, by = c("group","simmr_join"))



rownames12 =  data.frame(group = rep(c(1:12), each = 25), 
                      quantile = rep(rep(c("2.5", "25","50","75","97.5"), each = 5), by = 12),
                      type = rep(rep(c("dev", "ben","pel","sdiso1","sdiso2"), rep = 5),by=12), 
                      simmr_join = 1)
simmr_estimates12 = data.frame(value = unlist(CI_list[[1]]$quantiles)) %>%
  as.matrix() %>%
  cbind(rownames12) %>%
  left_join(mod.leg)
simmr_estimates12 = rbind(simmr_estimates12, simmr_estimates12) %>%
  mutate(community = rep(c(1,2), each = 300))

simmr_estimates12 %>% filter(group == 10)
  
simmr_estimates = rbind(simmr_estimates12, simmr_estimates3)







simmr_vol = species_vols %>% 
  mutate(community = as.numeric(community)) %>%
  left_join(simmr_estimates %>%
                                mutate(group = as.character(group)) %>%
                             filter(type == "ben",
                                    quantile == "50"), by = c("species" = "CODE", "community"))



## Here I'm estimmating the weighted average of benthic contibution to the food web

## viewing averages of simmr_vol 
simmr_vol %>% 
  group_by(community, species) %>%
  summarize(mean = mean(value))



simmr_vol %>%
#  filter(species != "SMB") %>%
  mutate(weighted.avg = total_vol * value) %>%
  group_by(post, community) %>%
  mutate(total_web_vol = sum(total_vol)) %>%
  ungroup() %>%
  group_by(post, community, total_web_vol) %>%
  summarize(wa = sum(weighted.avg, na.rm = T)) %>%
  ungroup() %>%
  mutate(benthic.contr = wa / total_web_vol) %>%
  ggplot(aes(x = as.factor(community), y = 100*benthic.contr, col = as.factor(community))) +
  stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )  +
  theme_minimal(base_size = 12) +
  ylab("Benthic Contribution (%)") +
  scale_x_discrete("Period",
                   labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.") ) +

  scale_color_manual("Community",
                    values = wes_palette("Darjeeling1", 
                                         type = "discrete", n = 3),
  labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  theme(legend.position = "non")

simmr_vol_graph = simmr_vol %>%
#  filter(species != "SMB") %>%
  mutate(weighted.avg = total_vol * value) %>%
  group_by(post, community) %>%
  mutate(total_web_vol = sum(total_vol)) %>%
  ungroup() %>%
  group_by(post, community, total_web_vol) %>%
  summarize(wa = sum(weighted.avg, na.rm = T)) %>%
  ungroup() %>%
  mutate(benthic.contr = wa / total_web_vol) %>% 
  select(-total_web_vol, -wa) %>%
  pivot_wider(names_from = community, values_from = benthic.contr) %>% 
  mutate(diff13 = `3` - `2`,
         diff12 = `2` - `1`,
         diff23 = `3` - `2`) %>%
  select(diff13, diff12, diff23) %>%
  reframe(mean13 = mean(diff13), low13 = quantile(diff13, .025), up13 = quantile(diff13, .975),
            mean12 = mean(diff12), low12 = quantile(diff12, .025), up12 = quantile(diff12, .975),
            mean23 = mean(diff23), low23 = quantile(diff23, .025), up23 = quantile(diff23), .975)
save(simmr_vol_graph, file = "Data/VaryingIsotopesData/simmr_vol_graph.RData")  









