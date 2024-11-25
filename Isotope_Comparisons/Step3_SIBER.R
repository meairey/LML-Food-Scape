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

data_simmr = full %>%
  select(C, N, period, Species)  %>%
  rename(iso1 = C, 
         iso2 = N, 
         community = period, 
         group = Species) 

legend = data_simmr %>%
  rename(TAXON = group) %>%
  mutate(group = as.numeric(as.factor(TAXON)), 
         community = as.numeric(as.factor(community))) %>%
  arrange(community, group) %>%
  select(community, TAXON, group) %>%
  unique() %>%
  rename(CODE = TAXON) %>%
  mutate(color = "black") %>%
  mutate(common = "n")

save(legend, file = "Data/IsotopeComparisons/simmr_legend.RData")

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
  select(iso1, iso2,group, community) 

siber.example <- createSiberObject(as.data.frame(data_siber)) ## Creates data object 

posterior <- siberMVN(siber.example, parms, priors) ## Generates posterior distributions 


overlap_list = list() 



for(h in 1:length(unique(data_siber$community))){
  
  overlap_list[[h]] = overlap(data_siber,h,30, posterior) %>% as.data.frame()
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
  na.omit() %>%
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
            q2.5 = round(quantiles_95(Values), digits = 2),
            q97.5 = round(quantiles_95(Values), digits = 2)) %>%
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

ellipse.area %>% ggplot(aes(x = group, y = area, fill = community)) + 
  geom_boxplot()



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
