# Introduction

This is an introduction to the components and the generation of the
posterior distributions for the foodweb landscapes of Little Moose.
Below you will find snippits of code and examples of the data frames
that are used during this process. Here we will look at the data frames
that are associated with the `step1_LML_source_L.R` file. However, the
graphs will display the results of the posterior generation for both the
early and late periods of this experiment.

------------------------------------------------------------------------

`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)`

\`\`\`{r, echo = F, include = F, warnings = F, message=F} \## Setup for
the markdown \# Libraries ———————- library(easypackages)
libraries(“snow”,“plotrix”, “SIBER”,“ellipse”,“mixtools”,
“mvtnorm”,“plot3D”,“scatterplot3d”,“scales”,“viridis”,“ggplot2”,
“gridExtra”, “MASS”,“plotly”, “knitr”,“tidyverse”, “mvtnorm”,“rjags”)



    ```{r, echo = F, message=F, warning=F}

    ## Load in data -----------------


    # Sources
    # For functions
    source("Data/Models/functions.R")
    # For starting conditions. This tutorial is for the late period.
    source("Landscapes/LateLandscapes/step1_LML_source_L.R")

    ## Posterior data frame --------------------
    ## Combining chains 

    load("Data/LateData/LateJagsOutput.RData")

    ## Chains generated from previous script 
    chains <- as.mcmc.list(jags_output.late) ## set chains
    chains_combined <- gtable_combine(chains) ## combine chains

## Define starting conditions for the model snippits.

For all periods of this experiment, we will be looking at 20 posterior
draws `n.posts`, 10 species `n_species`, 32 sites around the shore of
Little Moose Lake `n_sites`, and 4 years `n_years`. The early period
includes the years 2001-2004 and the late period includes the years
2019-2022.

\`\`\`{r} \# Model conditions n.posts = 20 n\_species = 10 n\_sites = 32
n\_years = year\_max - year\_min -1 \# seems missing 2023 data?


    We begin by setting names for the data frame that holds the posterior distributions. These names are composed of useful information that we can later filter by. For example, `N_1_1_1` is the posterior distribution for the estimated true abundance of species 1 at site 1 during the first year of that period. 
    ```{r, echo = F, message=F, warning=F}
    # Define names for the 
    names = (data.frame(colnames = colnames(jags_output.late[[1]]))  %>%
               mutate(colnames = str_replace_all(colnames, "\\[", "_"), 
                      colnames = str_replace_all(colnames, "\\]", ""),
                      colnames = str_replace_all(colnames,",", "_")))$colnames

    chain_dat = chains_combined %>%
      as.matrix() %>%
      as.data.frame() %>%
      set_names(names)


    head(chain_dat)[1:5]

But the entire posterior distribution is a large data frame. At this
point, we are not subsetting the posterior distribution and the
abundances are estimated at each site, at each year, for each species.
This results in a large number of columns.

`{r} dim(chain_dat)`

## Subsetting the compontents of the landscape generation.

We begin with sigma. Sigma represents the shape of the isotope ellipse
for each species.It is composed of a covariance matrix.

Below we see an example of the sigma data from the posterior
distribution.

\`\`\`{r, echo = F, messages = F, warning = F} \## Clean sigma estimates
sig.dat = chain\_dat %&gt;% select(contains(“sigma”, ignore.case =
TRUE)) %&gt;% as.data.frame() %&gt;% mutate(post = 1:length(.\[,1\]))
%&gt;% filter(post &lt;= n.posts) %&gt;% select(post, everything())
%&gt;% pivot\_longer(2:length(.\[1,\]), names\_to = “metric”, values\_to
= “sig.value”) %&gt;% separate(metric, into = c(“metric”, “species”,
“sig.order1”, “sig.order2”)) %&gt;% mutate(metric =
str\_remove\_all(metric, “\[0-9\]”)) %&gt;% unite(“names”, c(metric,
sig.order1, sig.order2)) %&gt;% pivot\_wider(names\_from = names,
values\_from = sig.value)

head(sig.dat)

    For `species == 1` and `post == 1`, see the below matrix as an example of what the covariance matrix looks like.

    ```{r, echo = F, message=F, warning = F}
    cov.mat = sig.dat %>% 
      filter(post == 1, species ==1)

    matrix(c(cov.mat$Sigma_1_1, cov.mat$Sigma_1_2, cov.mat$Sigma_2_1, cov.mat$Sigma_2_2),ncol = 2, nrow = 2)

The mass data is generated through the use of length and weight
relationships for species in Little Moose Lake. The one exception here
is brown bullhead. For brown bullhead, we use weight-length
relationships that are available on FishBase from previously published
research.Below, the estimated length-weight relationships are
represented by a colored curve. The observed length-weight data is
displayed.

\`\`\`{r, echo = FALSE, message = FALSE, warning = FALSE} weight.frame =
read.csv(“Data/weight\_frame.csv”)

weight.frame %&gt;% ggplot(aes(y = LENGTH, x = as.numeric(weight\_e),
col = SPECIES)) + geom\_point(aes(y = LENGTH, x = WEIGHT), col =
“black”, alpha = .5, key\_glyph = “rect”) + geom\_smooth() +
facet\_wrap(~SPECIES, scales = “free”) + theme\_minimal() +
theme(legend.position = “none”) + xlab(“Weight (g)”) + ylab(“Length
(mm)”)


    See below the structure of the mass data.

    ```{r, echo = F, messages = F, warning = F}
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

    head(mass.dat)

Comparing the masses of fishes bewteen the two periods. Data from the
posterior for both the late and the early period is loaded in here.

\`\`\`{r, message = F, echo = F, warning = F} \## Let’s compare the
early and late periods… load(“Data/posterior\_pre.csv”)
load(“Data/posterior\_late.csv”) load(“Data/posterior\_early.csv”) \##
Legend species.late = c(“BB”,
“CC”,“CS”,“LT”,“MM”,“PS”,“RS”,“SMB”,“SS”,“WS”) legend = data.frame(group
= c(1:10), species = species.late)

## Comparison of posterior parameter distributions

load(“Data/posterior\_early.csv”) load(“Data/posterior\_late.csv”)
load(“Data/posterior\_pre.csv”) \## Legend species= c(“BB”,
“CC”,“CS”,“LT”,“MM”,“PS”,“RS”,“SMB”,“SS”,“WS”) species.pre = c(“BB”,
“CC”,“CS”,“PS”,“SMB”,“SS”,“WS”) legend = data.frame(group =
as.character(c(1:10)), species = species, species.pre = c(species.pre,
rep(“NA”, 3)))

## Comparison of posterior parameter distributions

posterior.pre = posterior.pre %&gt;% rename(“group” = “species”) %&gt;%
left\_join(legend) %&gt;% mutate(community = 1) %&gt;% mutate(species =
species.pre) %&gt;% select(-species.pre) posterior.early =
posterior.early %&gt;% rename(“group” = “species”) %&gt;%
left\_join(legend) %&gt;% mutate(community = 2) %&gt;%
select(-species.pre) posterior.late = posterior.late %&gt;%
rename(“group” = “species”) %&gt;% left\_join(legend) %&gt;%
mutate(community = 3) %&gt;% select(-species.pre)

posterior = rbind(posterior.pre, posterior.early, posterior.late)

## 

posterior %&gt;% mutate(group =species) %&gt;%

left\_join(legend) %&gt;% ggplot(aes(x = as.factor(community), y =
mass.avg)) + geom\_boxplot() + \#scale\_y\_log10() +
facet\_wrap(~species, scales = “free”) + xlab(“Period”) + ylab(“Average
Mass (g)”) + theme\_minimal(base\_size = 12) + scale\_x\_discrete(labels
= c(“1” = “Pre”,“2” = “Early”, “3” = “Late”))


    See below for the snippit of the abundance data.
    ```{r, echo = F, messages = F, warning = F}
    # Clean abund estimates -- this does it year by year
    abund.dat =  chain_dat %>%
      as.data.frame() %>%
      select(contains("N", ignore.case = TRUE)) %>%
      
      mutate(post = 1:length(.[,1])) %>%
      filter(post <= n.posts) %>%
      select(post, everything()) %>% 
      #
      pivot_longer(2:((n_species * n_years * n_sites)+1),
                   names_to = "metric", 
                   values_to = "value") %>%
      separate(metric, into = c("metric", "site", "year", "species")) %>%
      mutate(species = as.numeric(species),
             site = as.numeric(site),
             post = as.numeric(post), 
             year = as.numeric(year)) %>%
      group_by(post, species, year) %>%
      summarize(tot_abund = sum(value)) %>%
      mutate(species = as.character(species))

    head(abund.dat)

See below for a comparison between the estimated total abundance for
both periods. \`\`\`{r, echo = F, messages = F, warning = F}

posterior %&gt;% mutate(group = (species)) %&gt;%

ggplot(aes(x = as.factor(community), y = tot\_abund)) +
geom\_boxplot() + scale\_y\_log10() + facet\_wrap(~species, scales =
“free”) + xlab(“Period”) + ylab(“Total Abundance”) +
theme\_minimal(base\_size = 12) + scale\_x\_discrete(labels = c(“1” =
“Pre”,“2” = “Early”, “3” = “Late”))


    ```{r, echo = F, messages = F, warning = F}
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

    head(mu.dat)

For this analysis, the mean and sigma of each community’s isotopic web
is centered around zero with a SD of 1. Because of this, the webs appear
in an transformed space. This is to assist in making comparisons across
communities so that the variation within isotopic landscapes of each
system don’t largely influence the webs. Below, you can see a
visualization of the mean centroids and ellipse shape for each species
during each period. Note - that the Pre and Early communities use the
same isotopic data. However, several species (i.e. lake trout, central
mudminnow, and rainbow smelt) are not observed during that period, so
they are excluded from that isotopic web.

\`\`\`{r, echo = F, messages = FALSE, warning = FALSE}

load(“Data/MeanEllipses.RData”) ellip.mean.filtered %&gt;% mutate(spp =
as.character(spp)) %&gt;% left\_join(legend,by = c(“spp” = “group”))
%&gt;% mutate(community =as.character(community)) %&gt;%
\#mutate(community = case\_when(community == 1 ~ “Pre”, community == 2 ~
“Early”, community == 3 ~ “Post”)) %&gt;% ggplot(aes(x= xax, y = yax,
col =as.factor(species))) + theme\_minimal(base\_size = 14) +
stat\_ellipse(lwd = 1) + facet\_wrap(~community, labeller =
labeller(community = c(“1” = “Pre”, “2” = “Early”, “3” = “Late”))) +
labs(col = “Species”)



    Finally, all the previous steps are combined. When applicable (i.e. where abundances are estimated across multiple years) the estimates are summarized by mean.

    See below for the structure of the final posterior data frame is.
    ```{r, echo = FALSE, messages = FALSE, warning = FALSE, include = F}


    ## Final posterior frame
    posterior.late = left_join(sig.dat, mass.dat) %>% 
      left_join(abund.dat) %>%
      left_join(mu.dat) %>%
      group_by(post, species) %>%
      select(year, everything()) %>%
      summarize(across(2:9, mean, na.rm = TRUE))

`{r, echo = FALSE, warning = FALSE, message=FALSE} head(posterior.late)`
