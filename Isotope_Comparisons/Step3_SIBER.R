## SIBER
library(tidyverse)
library(SIBER)

## Needed Functions
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

`%nin%` = Negate(`%in%`)
#----------------------------- Params ---------------------------------------- 

x.p = seq(-38,-20,length.out=100); 
y.p = seq(0,9, length.out = 100) 
parms <- list();

parms$n.iter=2*10^4; 
parms$n.burnin=1*10^3; 
parms$n.thin=10
parms$n.chains=3     
priors=list();
priors$R=1*diag(2);
priors$k=2; 
priors$tau.mu=1.0E-3 # Vague priors 
Nsamples=1000
n.posts <- 1000;
p.ell <- 0.90 # How much data to include? Standard ellipses --> p.ell = .9
n.points = 1000

source("Isotope_Comparisons/IsotopeFunctions.R")
## Must have legend

load(file = "Data/IsotopeComparisons/IsotopeData.RData")

## Regular data formatting
data_simmr = data %>%
  select(D13C, D15N, period, TAXON)  %>%
  rename(iso1 = D13C, 
         iso2 = D15N, 
         community = period, 
         group = TAXON) 

## Trying to split out SMB1/2

data_simmr = data %>% 
  unique() %>%
  filter(TAXON == "SMB") %>%
  bind_rows(.,data) %>%
  arrange(ITEM_N) %>%
  group_by(ITEM_N, ISO_FISH_N) %>%
  mutate(count = c(1:length(TAXON))) %>%

  select(count, TAXON, LENGTH, everything()) %>%
  ungroup() %>%

  mutate(TAXON = case_when(TAXON == "SMB" & LENGTH < 200  & LENGTH > 100 & count == 2 ~ "SMB2" ,
                           TAXON == "SMB" & LENGTH > 0 & LENGTH < 100 & count == 2 ~ "SMB1",
                           TAXON == "SMB" & LENGTH >= 200 & count == 2 ~ "SMB3",
                           TAXON == "SMB" & count == 1 ~ "SMB",
                           TAXON != "SMB" ~ TAXON)) %>%
  group_by(TAXON, period) %>%
  filter(TAXON != "NA") %>%
  select(D13C, D15N, period, TAXON)  %>%
  rename(iso1 = D13C, 
         iso2 = D15N, 
         community = period, 
         group = TAXON) %>%
  ungroup() 


data_simmr %>% 
  group_by(group, community) %>%
  summarize(n = n())

data_simmr$community %>% unique()


data_simmr %>%
  filter(group %in% c("SMB1", "SMB2", "SMB3")) %>%
  ggplot(aes(x = iso1, y = iso2, col = group)) + 
  geom_point() + 
  stat_ellipse() +
  facet_wrap(~community)

#save(legend, file = "Data/IsotopeComparisons/simmr_legend.RData")

data_siber = data_simmr %>%
  mutate(group = as.numeric(as.factor(group)),
         community = as.numeric(as.factor(community)),
         iso1 = as.numeric(iso1), 
         iso2 = as.numeric(iso2)) %>%
  na.omit() %>%
  group_by(community, group) %>%
  mutate(total = n()) %>%
  filter(total > 3) %>%
  arrange(community, group) %>%
  select(iso1, iso2,group, community) ## Does still have common shiner



data_siber %>%
  left_join(legend) %>%
  filter(CODE == "CS")

## Create Legend
legend.a = data_simmr %>%

  rename(TAXON = group) %>%
  mutate(group = as.numeric(as.factor(TAXON)), 
         community  = as.numeric(as.factor(community))) %>%
  arrange(community, group) %>%
  select(community, TAXON, group) %>%
  unique() %>%
  rename(CODE = TAXON) %>%
  mutate(common = "n")

col_join = data.frame(color = c("#FF8C00", "#FFD700", "#87CEEB", "#4682B4", 
  "#6A5ACD", "#F4A460", "#32CD32", "#FF4500", 
  "#8B0000", "#FF69B4", "#FF6347", "#40E0D0", 
  "#2E8B57", "#D2691E", "#DAA520", "#FFB6C1", "#800080" , "black"),
  CODE = unique(legend.a$CODE))

legend = legend.a %>% left_join(col_join)
save(legend, file = "Data/Legend.col.RData")

## Siber Object
siber.example <- createSiberObject(as.data.frame(data_siber)) ## Creates data object 

posterior <- siberMVN(siber.example, parms, priors) ## Generates posterior distributions 


overlap_list = list() 

## 100 posteriors takes about 2hrs

for(h in 1:length(unique(data_siber$community))){
  
  overlap_list[[h]] = overlap(data_siber,h, 100, posterior) %>% as.data.frame()
}

# Code to reduce list to matrix 
overlap_data = Reduce(full_join, overlap_list) %>%
  select('SppPair', everything())


save(overlap_data, file = "Data/IsotopeComparisons/OverlapData.RData")

save(posterior, file = "Data/IsotopeComparisons/Posterior.RData")


load("Data/IsotopeComparisons/OverlapData.RData")
load("Data/IsotopeComparisons/Posterior.RData")


overlap.long = overlap_data %>% 
  pivot_longer(2:length(.[1,]),
               names_to = "Community", 
               values_to = "Values") %>% 
  #na.omit() %>%
  rename("Species_Pair" = `SppPair`) %>% 
  separate(Community, into = c("c", "Community", "post")) %>%
  select(-c) %>%
  mutate(Values = as.numeric(Values)) 


#library(tidyr)
#library(tidyverse)


## table 

df = overlap.long%>% ungroup() 
df$group1 = sapply(strsplit(df$Species_Pair, "v"), function(x) min(x))
df$group2 <- sapply(strsplit(df$Species_Pair, "v"), function(x) max(x))
df = df %>%
  mutate(group1 = parse_character(group1), group2 = parse_character(group2))
df$sorted_comparison <- apply(df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "_"))


df.long = df %>% 
  select(sorted_comparison, everything(), 
         -Species_Pair, -group1, -group2) %>%
  unique()


dec.comp = df.long %>% select(sorted_comparison, Community) %>% unique() %>% group_by(sorted_comparison) %>%
  summarize(count = n()) %>%
  filter(count > 1)




df.long %>% filter(!grepl("SMB", sorted_comparison)) %>%
  mutate(period = case_when(Community == 1 ~ "early", Community == 2 ~ "late")) %>%
  group_by(Community, post, period) %>%
  summarize(mean_value = mean(Values)) %>%
  ggplot(aes(x = period, y = mean_value, fill = period)) + 
  geom_boxplot() + 
  theme_minimal() + 
  ylab("Average Pairwise Overlap")


overlap.df = df.long
save(overlap.df, file = "Data/IsotopeComparisons/overlap.df.RData")


overlap_test = (df.long %>%
                  filter(!grepl("SMB", sorted_comparison)) %>%
                  mutate(period = case_when(Community == 1 ~ "early", Community == 2 ~ "late")) %>%
                  group_by(Community, post, period) %>%
                  summarize(mean_value = mean(Values)) %>%
                  ungroup() %>%
                  select(-Community) %>%
                  pivot_wider(names_from = period, values_from = mean_value) %>%
                  mutate(diff = early - late) %>%
                  ungroup() )$diff %>% as.numeric() 



## Pairwise matrix table of all overlap outputs
### This is sent into a csv which is then modified in excel as an .xlsx file
df = df.long %>%
  mutate(period = case_when(Community == 1 ~ "early", Community == 2 ~ "late")) %>%
  group_by(sorted_comparison, period) %>%
  summarize(mean = round(mean(Values), digits = 2), 
            q2.5 = round(quantiles_95(Values)[1], digits = 2),
            q97.5 = round(quantiles_95(Values)[5], digits = 2)) %>%
  mutate(mean = paste(mean, "[", q2.5,",",q97.5,"]")) %>%
  separate(sorted_comparison, into = c("Species1","Species2"))

df  


# Get unique species
species <- unique(c(df$Species1, df$Species2))

# Create an empty 13x13 matrix
matrix_13x13 <- matrix(NA, nrow = length(species), ncol = length(species),
                       dimnames = list(species, species))

# Populate the matrix
for (i in seq_len(nrow(df))) {
  row <- df$Species1[i]
  col <- df$Species2[i]
  value <- df$mean[i]
  period <- df$period[i]
  
  if (period == "early") {
    matrix_13x13[row, col] <- value
  } else if (period == "late") {
    matrix_13x13[col, row] <- value
  }
}

# Output the matrix
write.csv(matrix_13x13, file = "Data/OverlapComparisons.csv")
print(matrix_13x13)


## Ellipse area just use full posterior with full length of species/ellipses  
ellipse.area = siberEllipses(posterior) %>% 
  as.data.frame() %>%
  rename_with(~ names(posterior)) %>% 
  mutate(post_n = seq(1:length(.[,1]))) %>%
  pivot_longer(1:length(posterior), 
               names_to = "comm",
               values_to = "area") %>%
  separate(comm, into = c("community", "group")) %>%
  
  group_by(post_n, community) %>%
  mutate(total_area = sum(area)) %>%
  mutate(relative_area = area / total_area) %>% 
  group_by(community, group) 


save(ellipse.area, file = "Data/IsotopeComparisons/EllipseArea.RData")

load("Data/IsotopeComparisons/EllipseArea.RData")

med_area = ellipse.area %>% 
  filter(relative_area >= quantile(relative_area, probs = 0.025) & 
           relative_area <= quantile(relative_area, probs = 0.975)) %>% 
  ungroup() %>%
  group_by(group, post_n, community) %>% 
  #summarize(med_area = mean(relative_area)) %>%
  summarize(med_area = mean(area)) %>%
  mutate(group = as.numeric(group), 
         community = as.numeric(community)) %>%
  left_join(legend) %>%
  mutate(CODE = str_replace(CODE, "NPD", "PD")) %>% 
  mutate(CODE = str_replace(CODE, "BNM", "BM")) %>%
  mutate(period = case_when(community == 1 ~ "early", 
                            community == 2 ~ "late"))

save(med_area, file ="Data/IsotopeComparisons/med_area.RData")


### Plotting isotope webs

posterior

source("../thermal-guild-otolith/isotope_functions_update.R")

### Ontogenetic shifts in bass?
full_corrected%>% 
  unique() %>%
  filter(TAXON == "SMB") %>%
  bind_rows(.,data) %>%
  arrange(ITEM_N) %>%
  group_by(ITEM_N, ISO_FISH_N) %>%
  mutate(count = c(1:length(TAXON))) %>%

  select(count, period,TAXON, LENGTH, everything()) %>%
  select(TAXON, LENGTH, D13C, D13C_c, D15N, D15N_c) %>%
  ungroup() %>%
  na.omit() %>%

  mutate(TAXON = case_when(TAXON == "SMB" & LENGTH >= 100 & LENGTH < 200 & count == 2 ~ "SMB2" ,
                           TAXON == "SMB" & LENGTH < 100 & count == 2 ~ "SMB1" ,
                           TAXON == "SMB" & LENGTH >= 200 & count == 2 ~ "SMB3",
                           TAXON == "SMB" & count == 1 ~ "SMB",
                           TAXON != "SMB" ~ TAXON)) %>%
  group_by(TAXON, period) %>%
  filter(TAXON != "NA") %>%
  select(D13C_c, D15N_c, period, TAXON)  %>%
  rename(iso1 = D13C_c, 
         iso2 = D15N_c, 
         community = period, 
         group = TAXON) %>%
  ungroup() %>%
  filter(group %in% c("SMB1", "SMB2", "SMB3")) %>%
  ggplot(aes(x = iso1, y = iso2, col = group)) + 
  geom_point() + 
  stat_ellipse(level = .4) +
  facet_wrap(~community, scales = ) +
  theme_minimal(base_size = 13)

## Quantifying SMB overlap


df %>% 
  filter(grepl("SMB",group1), 
         grepl("SMB", group2)) %>%
  filter(group1 != "SMB",
         group2 != "SMB") %>%
  group_by(Community, sorted_comparison) %>%
  na.omit() %>%
  summarize(mean_overlap = quantiles_95(Values)[3],
            lower = quantiles_95(Values)[1],
            upper = quantiles_95(Values)[5])

