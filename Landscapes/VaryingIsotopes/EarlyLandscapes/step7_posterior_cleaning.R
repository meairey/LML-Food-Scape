# Libraries ----------------------
library(easypackages)
libraries("snow","plotrix", "SIBER","ellipse","mixtools",
          "mvtnorm","plot3D","scatterplot3d","scales","viridis","ggplot2",
          "gridExtra", "MASS","plotly", "knitr","tidyverse", "mvtnorm")
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML-Food-Scape/")
library(coda)
## Load in data -----------------


# Source--------------
# For functions
source("Data/Models/functions.R")
source("Landscapes/VaryingIsotopes/EarlyLandscapes/step1_LML_source.R")

n.posts = 3000
n_species = 9 ## Change back to 8 if you remove LT
n_years = year_max - year_min
n_sites = 32
## Posterior data frame --------------------
## Combining chains 

load("Data/VaryingIsotopesData/EarlyData/EarlyJagsOutput.RData")

chains <- as.mcmc.list(jags_output.early) ## set chains
chains_combined <- gtable_combine(chains) ## combine chains



rhat_values <- gelman.diag(chains, multivariate = FALSE)$psrf[, "Point est."]

rhat_values %>% as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "rowname") %>%
  separate(rowname,  into = c("parameter", "index"), remove = F) %>%
  #  select(parameter, V1) %>%
  filter(V1 > 1.1)

chains_combined <- gtable_combine(chains) ## combine chains

chain1 = as.matrix(chains[[1]]) %>% as.data.frame() %>%
  mutate(chain = 1)

chain2 = as.matrix(chains[[2]]) %>% as.data.frame() %>%
  mutate(chain = 2)

combined_density = rbind(chain1, chain2) 
for(i in c(1:n_species)){
  g = combined_density %>% ggplot(aes(x = combined_density[,i], col = as.factor(chain))) +
    geom_density(alpha = .05) +
    xlab(paste("early", i))
  print(g)
}



# Define names for the 
names = (data.frame(colnames = colnames(jags_output.early[[1]]))  %>%
           mutate(colnames = str_replace_all(colnames, "\\[", "_"), 
                  colnames = str_replace_all(colnames, "\\]", ""),
                  colnames = str_replace_all(colnames,",", "_")))$colnames

chain_dat = chains_combined %>%
  as.matrix() %>%
  as.data.frame() %>%
  set_names(names)

## Clean sigma estimates
sig.dat = chain_dat %>%
  select(contains("sigma", ignore.case = TRUE))  %>%
  as.data.frame() %>% 
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]), names_to = "metric", 
               values_to = "sig.value") %>%
  separate(metric,
           into = c("metric", "species", 
                    "sig.order1", "sig.order2")) %>% 
  mutate(metric = str_remove_all(metric, "[0-9]")) %>%
  unite("names", c(metric, sig.order1, sig.order2)) %>% 
  pivot_wider(names_from = names, values_from = sig.value)


## Clean length estimates
mass.dat = chain_dat %>%
  select(contains("avg_mass", ignore.case = TRUE))  %>%
  as.data.frame() %>%
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]),
               names_to = "metric", values_to = "mass.avg") %>%
  separate(metric, into = c("metric1","metric", "species"), sep = "_") %>%
  select(-metric1, -metric)

# Figure for looking at the mass data
mass.dat %>%
  rename("group" = "species") %>% 
  mutate(group = as.numeric(group)) %>%
  left_join(legend, by = c(group = "group")) %>%
  ggplot(aes(x = species, y = exp(mass.avg)))  +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_log10() +
  ylab("Avg. Mass (g)") +
  xlab("Species") +
  labs(col = "Species")


# Clean abund estimates -- this does it year by year
abund.dat =  chain_dat %>%
  as.data.frame() %>%
  select(contains("N", ignore.case = TRUE)) %>%
  
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything() ) %>% 
  #
  pivot_longer(2:length(.[1,]),
               names_to = "metric", 
               values_to = "value") %>%
  separate(metric, into = c("metric", "site",  "species")) %>%
  mutate(species = as.numeric(species),
         site = as.numeric(site),
         post = as.numeric(post)) %>%
  group_by(post, species) %>%
  summarize(tot_abund = sum(value)) %>%
  mutate(species = as.character(species))

# Figure for looking at the abundance data
abund.dat %>% 
  mutate(species = as.numeric(species)) %>%
  rename(group = species) %>%
  left_join(legend) %>%
  ggplot(aes(x = species, y = tot_abund)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab("Total Abundance") +
  xlab("Year Index") +
  labs(col = "Species")


## Table for looking at the abundance data
abund.sum.early = abund.dat %>% 
  mutate(species = as.numeric(species)) %>%
  rename(group = species) %>%
  left_join(legend)  %>%
  group_by(species) %>%
  summarize(tot_abund = mean(tot_abund)) 

save(abund.sum.early, file =  "Data/VaryingIsotopesData/early.abund.RData")
remove(abund.sum.early)
load(file =  "Data/VaryingIsotopesData/early.abund.RData")


# Clean mu estimates
mu.dat = chain_dat %>%
  as.data.frame() %>% 
  select(contains("mu", ignore.case = TRUE)) %>%
  mutate(post = 1:length(.[,1])) %>%
  filter(post <= n.posts) %>%
  select(post, everything()) %>%
  pivot_longer(2:length(.[1,]), names_to = "metric", values_to = "mu.value") %>% 
  separate(metric, into = c("metric", "species", "isotope")) %>% 
  unite("names", c(metric, isotope)) %>%
  pivot_wider(names_from = names, values_from = mu.value)



## Final posterior frame
posterior.early = left_join(sig.dat, mass.dat) %>% 
  left_join(abund.dat) %>%
  left_join(mu.dat) %>%
  group_by(post, species) 

posterior.early %>%
  group_by(species) %>%
  summarize(mean_N = mean(mu_1), mean_ab = mean(mass.avg))


species.early
save(file = "Data/VaryingIsotopesData/posterior_early.csv", posterior.early)
