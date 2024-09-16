library(tidyverse)
## This file is for working on the SIMMR/mixing model/baselines for the food web landscapes

# This iso data is cleaned up in a script and written into a csv from `AFRP_SIA.Rmd`
iso_data = read.csv("Data/LML_iso_corrected.csv")



## Baselines 

## All potential baseline items
baselines.full = full %>% 
  mutate(D13C = as.numeric(D13C),
         D15N = as.numeric(D15N)) %>%
  mutate(community = case_when(YEAR < 2005 ~ "early", 
                               YEAR > 2019 ~ "late")) %>%
  filter(GROUP %in% c("ZOOP","SNAIL","INSECT"))
  



## bind together all three potential baselines
baselines.summary = baselines.full %>% 
  group_by(community, GROUP)  %>%
  summarize(mean_d13C = mean(D13C), mean_d15N = mean(D15N), sd_C = sd(D13C), sd_N = sd(D15N))
baselines.summary



## Create zooplankton data for late period (waiting to run this)
fake_late.zoop = data.frame(GROUP = "ZOOP", community = "late", mean_d13C = baselines.summary$mean_d13C[3], mean_d15N = (baselines.summary$mean_d15N[3] - 2.001345),
                            sd_C = baselines.summary$sd_C[3], sd_N = baselines.summary$sd_N[3])
# Join real and fake data together
baselines = baselines.summary %>% rbind(fake_late.zoop) %>%
  ungroup() %>%
  mutate(S_N = c("benthic", "benthic", "pelagic", "benthic", "bethic", "pelagic"))
baselines

# Creating the baselines that go into SIMMR
baseline_simmr = baselines %>% filter(GROUP %in% c("SNAIL", "ZOOP"))
# Creating the data table that goes into SIMMR
data_simmr = full %>% filter(GROUP == "FISH") %>%
  mutate(community = case_when(YEAR < 2005 ~ "early", 
                               YEAR > 2019 ~ "late")) %>% 
  na.omit() %>% 
  select(TAXON, community, D13C, D15N) %>%
  mutate(D13C = as.numeric(D13C), D15N = as.numeric(D15N)) %>%
  unique()


library(simmr)
CI_list = list()
for(i in 1:length(data_simmr$community %>% unique())){
  
  data_s = data_simmr %>% filter(community == (data_simmr$community %>% unique())[i]) %>% ungroup() %>% as.data.frame()
  
  baseline_s = baseline_simmr %>% filter(community == (data_simmr$community %>% unique())[i]) %>% ungroup() %>% as.data.frame()
  
  S_N = c("Benthic", "Pelagic")
  
  
  simmr = simmr_load(
    mixtures = data_s[,3:4] %>% as.matrix(),
    source_names = S_N,
    source_means = baseline_s[,3:4],
    source_sds = baseline_s[,5:6], 
    group = data_s[,1]
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
  plot(simmr_out, type = "histogram")
  
  compare_groups(simmr_out,
                 groups = 1:length(unique((data_s$TAXON))),
                 
                 source_name = "Benthic")
  
  
  
  CI_list[[i]] = summary(simmr_out, 
                         group = 1:length(unique(data_s$TAXON)), 
                         type = "quantiles")
  
  
}

#


## This is how you get out the posterior
post = simmr_out$output$CC$BUGSoutput$sims.array  

# Get the sims.array from the simmr output
sims_array <- simmr_out$output$CC$BUGSoutput$sims.array

# Flatten the array by combining the chains and iterations
# The dimensions are usually [iterations, chains, parameters]
flattened_sims <- apply(sims_array, 3, function(x) as.vector(x))

# Convert the flattened array to a data frame
flattened_df <- as.data.frame(flattened_sims)

# Check the result
head(flattened_df)
# Check the result
head(combined_quantiles)

# Community 1 simmr data
community1<- do.call(rbind, CI_list[[1]]$quantiles) %>% as.data.frame() %>%
  mutate(ID = rep(names(CI_list[[1]]$quantiles), each = 5)) %>%
  mutate(community = 1)
# Community 2 simmr data
community2<- do.call(rbind, CI_list[[2]]$quantiles) %>% as.data.frame() %>%
  mutate(ID = rep(names(CI_list[[2]]$quantiles), each = 5)) %>%
  mutate(community = 2)

## Modify simmr resutls into table so we can index for trophic position
simmr_results = rbind(community1, community2) %>%
  rownames_to_column(var = "rowname") %>%
  mutate(rowname = gsub("[0-9]", "", rowname))%>%
  mutate(rowname = gsub("\\.", "", rowname)) %>%
  filter(rowname == "Benthic") %>%
  select(-`25%`, -`75%`) %>%
  mutate(community= case_when(community == 1 ~ "early", community == 2 ~ "late")) %>%
  rename(group = ID)


## Create look-up table to join the appropiate group#s to each taxon
join_names = data_simmr %>% select(TAXON, community) %>% unique() %>%
  arrange(community, TAXON) %>%
  group_by(community) %>%
  mutate(group = paste("group_", 1:length(TAXON), sep = "")) 



tp =  data_simmr %>% left_join(join_names) %>%
  left_join(simmr_results, by = c("community", "group")) %>%
  mutate(TP = 2 + ((D15N - (baseline_simmr %>% filter(GROUP == "SNAIL" & community == community))$mean_d15N)*`50%` + 
                     (D15N - (baseline_simmr %>% filter(GROUP == "ZOOP" & community == community))$mean_d15N)*(1-`50%`)) / 3.4)


tp %>% ggplot(aes(x = TAXON, y = TP, fill = community)) +
  geom_boxplot() + theme_minimal() + ylab("Trophic Position")
