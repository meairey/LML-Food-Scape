## Varied and unvaried estimates



load(file = "Data/VaryingIsotopesData/MeanEllipses.RData")

load(file = "Data/UnvariedIsotopeData/MeanEllipses.RData")


e.un = ellip.mean.filtered.un %>% mutate(method = "un")
e.va = ellip.mean.filtered %>% mutate(method = "va")

dat = rbind(e.un, e.va)


dat %>% 
  
  mutate(spp = as.character(spp)) %>%
  left_join(legend,by = c("spp" = "group")) %>% 
  ggplot(aes(x= xax, y = yax, col =as.factor(method))) +
  geom_jitter(alpha = .15) + 
  facet_wrap(~species +community) 




## Filtered heatmap with landscape for all species
dat %>%
  #left_join(back.trace) %>%
  #mutate(xax = xax*sd_C + mean_C, yax = yax*sd_N + mean_N) %>%
  group_by(community, xax, yax, method) %>%
  summarize(vol = sum(string)) %>%
  ungroup() %>%
  ggplot(aes(x = xax, y = yax, col = (vol))) +
  geom_point(size = 5) +
  
  scale_color_viridis() +
  theme_minimal() + 
  labs(fill = "Z height") +
  xlab("d13C") +
  ylab("d15N") + 
  facet_wrap(~community + method)


## Full landscape 
load("Data/UnvariedIsotopeData/filtered_ellipses.RData")
load("Data/VaryingIsotopesData/filtered_ellipses.RData")

d.un = filtered_data.un %>% mutate(method = "un")
d.va = filtered_data %>% mutate(method = "va")
d = rbind(d.un, d.va)

# Volume metrics and boxplot ---------------------------------


ellip.filtered = d %>%
  mutate(spp = as.character(spp)) %>%
  left_join(legend, by = c("spp" = "group")) %>%
  rename("code" = "species")


# Total volume of community 
total_vols = ellip.filtered %>%
  group_by(code, post, community, method) %>%
  summarise(total_vol = sum(string)) 


## Total volume of each species individually
total_vols %>%
  as.data.frame() %>%
  filter(community ==3) %>%
  complete(code, post, community, method) %>%
  replace_na(list(total_vol = 0)) %>%
  ggplot(aes(x = code, y = total_vol +1, col = method)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Total Volume") +
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_log10() +
  theme_minimal(base_size = 14) +
  labs(col = "Period") +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  facet_wrap(~code, scales = "free")


# Total volume of the entire community


total_vols %>% 
  ungroup() %>%
  group_by(post, community, method) %>%
  summarize(total_vol = sum(total_vol)) %>%
  ggplot(aes(x = community %>% as.factor(), y = total_vol, col = method)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = c("1" = "Pre","2" = "Early", "3" = "Late")) +
  ylab("Total Volume") +
  xlab("Community") +
  theme(legend.position = "none") 

